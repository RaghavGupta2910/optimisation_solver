// QP engine benchmark harness.
//
// Builds a few representative test problems (small unconstrained, medium
// box-constrained, ill-conditioned with Ruiz scaling, large sparse) and
// reports solve time, iteration count, and primal/dual residuals.
//
// Run: ./qp_bench [cpu|cuda|auto]
//
// The optional backend argument applies to every case (default cpu). With the
// hybrid CUDA backend each case also reports where the time went: the CPU KKT
// solves, device setup, and the bytes moved each way.
//
// Run: ./qp_bench sweep <cpu|cuda|auto> <n> [nnzPerRow] [repeats] [threads]
//
// Times one generated QP of size n end to end for backend comparisons; see
// sweep() below.
#include "qp/admm_solver.h"
#include "qp/qp_model.h"
#include "qp/qp_solver.h"
#include "qp/qp_types.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>

// std::chrono::steady_clock and QueryPerformanceCounter on this MinGW host
// both produce out-of-scale or negative durations. std::clock() is portable
// and good enough for solve-time reporting on sub-second problems, which
// covers everything in the harness. Larger problems should use a real
// wall-clock profile elsewhere.
namespace bench {
double nowSeconds() {
    return static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
}
}  // namespace bench

namespace {

qp::SparseMatrix mkMat(int rows, int cols,
    const std::vector<double>& rv,
    const std::vector<double>& cv,
    const std::vector<double>& vv) {
    return qp::SparseMatrix::fromTriplets(rows, cols, rv, cv, vv);
}

void section(const char* name) {
    std::printf("\n=== %s ===\n", name);
}

struct BenchCase {
    const char* name;
    qp::QpModel model;
    qp::AdmmOptions options;
};

qp::QpModel smallUnconstrained() {
    // min 0.5 * ||x||^2 + 1^T x, n=20, no constraints.
    qp::QpModel m;
    const int n = 20;
    std::vector<double> rv(n), cv(n), vv(n, 1.0);
    for (int j = 0; j < n; ++j) { rv[j] = j; cv[j] = j; }
    m.P = mkMat(n, n, rv, cv, vv);
    m.A = mkMat(0, n, {}, {}, {});
    m.q.assign(static_cast<std::size_t>(n), 1.0);
    return m;
}

qp::QpModel mediumBoxed(int n) {
    // min 0.5 * x^T diag(1..n) x - sum(x_i), subject to 0 <= x <= 1.
    qp::QpModel m;
    std::vector<double> rv, cv, vv;
    rv.reserve(static_cast<std::size_t>(n));
    cv.reserve(static_cast<std::size_t>(n));
    vv.reserve(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) { rv.push_back(j); cv.push_back(j); vv.push_back(j + 1.0); }
    m.P = mkMat(n, n, rv, cv, vv);
    m.A = mkMat(n, n, rv, cv, std::vector<double>(static_cast<std::size_t>(n), 1.0));
    m.q.assign(static_cast<std::size_t>(n), -1.0);
    m.l.assign(static_cast<std::size_t>(n), 0.0);
    m.u.assign(static_cast<std::size_t>(n), 1.0);
    return m;
}

qp::QpModel illConditionedRanged(int n) {
    // 1e8 * I Hessian, badly scaled constraint matrix, box constraints.
    qp::QpModel m;
    std::vector<double> rv, cv, vv;
    rv.reserve(static_cast<std::size_t>(n));
    cv.reserve(static_cast<std::size_t>(n));
    vv.reserve(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) { rv.push_back(j); cv.push_back(j); vv.push_back(1e8); }
    m.P = mkMat(n, n, rv, cv, vv);

    // A = identity scaled by wildly different factors per row.
    std::vector<double> aR, aC, aV;
    for (int i = 0; i < n; ++i) {
        aR.push_back(i); aC.push_back(i);
        aV.push_back(std::pow(10.0, static_cast<double>((i % 7) - 3)));
    }
    m.A = mkMat(n, n, aR, aC, aV);

    m.q.assign(static_cast<std::size_t>(n), 0.5);
    m.l.assign(static_cast<std::size_t>(n), -1.0);
    m.u.assign(static_cast<std::size_t>(n), 1.0);
    return m;
}

qp::QpModel largeSparse(int n, double density) {
    // Random sparse P and A.  Tests the sparse Cholesky path when n > 500.
    qp::QpModel m;
    std::mt19937 g(42);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    std::uniform_int_distribution<int> colPick(0, n - 1);
    std::uniform_real_distribution<double> densityDist(0.0, 1.0);

    std::vector<double> pR, pC, pV;
    for (int i = 0; i < n; ++i) {
        // Diagonal entry to keep P positive definite.
        pR.push_back(i); pC.push_back(i); pV.push_back(1.0 + std::abs(val(g)));
        for (int j = 0; j < n; ++j) {
            if (i != j && densityDist(g) < density) {
                pR.push_back(i); pC.push_back(j); pV.push_back(val(g));
            }
        }
    }
    m.P = mkMat(n, n, pR, pC, pV);

    // m = n/2 equality-ish rows (boxed).
    const int mrows = n / 2;
    std::vector<double> aR, aC, aV;
    for (int i = 0; i < mrows; ++i) {
        int nnz = 1 + (colPick(g) % 5);
        for (int k = 0; k < nnz; ++k) {
            aR.push_back(i); aC.push_back(colPick(g)); aV.push_back(val(g));
        }
    }
    m.A = mkMat(mrows, n, aR, aC, aV);
    m.q.assign(static_cast<std::size_t>(n), 0.0);
    m.l.assign(static_cast<std::size_t>(mrows), -0.5);
    m.u.assign(static_cast<std::size_t>(mrows), 0.5);
    return m;
}

qp::ComputeBackend gBackend = qp::ComputeBackend::Cpu;

void run(const char* label, const BenchCase& bc) {
    qp::AdmmOptions options = bc.options;
    options.backend = gBackend;
    const double t0 = bench::nowSeconds();
    auto r = qp::AdmmSolver(bc.model, options).solve();
    double ms = 1000.0 * (bench::nowSeconds() - t0);
    // Sanitise: clock() can return -1 (failure sentinel) or negative values
    // on some MinGW builds. Clamp to 0 if the result is non-finite or negative.
    if (!(ms >= 0.0)) ms = 0.0;
    // Clamp negative or NaN times (can occur due to clock() quantization
    // on some hosts) to 0. All solve times here are well below 1 ms anyway.
    std::printf("  %-30s  status=%-9s  iters=%4d  obj=%12.6g  "
                "rPri=%10.3e  rDual=%10.3e  time=%7.2f ms\n",
                label, qp::toString(r.status), static_cast<int>(r.iterations),
                r.primalObjective, r.primalResidual, r.dualResidual, ms);
    if (gBackend != qp::ComputeBackend::Cpu) {
        std::printf("  %-30s  backend=%s  kkt solve=%.3g s  kkt factor=%.3g s  setup=%.3g s  "
                    "H2D=%lld B  D2H=%lld B  syncs=%lld\n",
                    "", qp::toString(r.executedBackend), r.kktSolveSeconds, r.kktFactorSeconds,
                    r.backendProfile.setupSeconds,
                    static_cast<long long>(r.backendProfile.hostToDeviceBytes),
                    static_cast<long long>(r.backendProfile.deviceToHostBytes),
                    static_cast<long long>(r.backendProfile.synchronisations));
    }
}

// Scalable convex QP for backend sweeps: n variables, n constraint rows with
// nnzPerRow entries each inside a band of width 2*nnzPerRow around the
// diagonal, so the CPU KKT factor stays banded and its cost grows linearly
// with n. P is a diagonally dominant tridiagonal (positive definite), and the
// row bounds are built around a reference point, so every instance has an
// optimum. Rows cycle through equalities, one-sided and ranged constraints.
qp::QpModel bandedSweepModel(int n, int nnzPerRow, unsigned seed) {
    std::mt19937 g(seed);
    std::uniform_real_distribution<double> val(-1.0, 1.0);
    const int band = std::min(n, 2 * nnzPerRow);
    std::uniform_int_distribution<int> offset(0, band - 1);

    std::vector<double> pR, pC, pV;
    for (int j = 0; j < n; ++j) {
        pR.push_back(j); pC.push_back(j); pV.push_back(1.0 + std::abs(val(g)));
        if (j + 1 < n) {
            const double v = 0.25 * val(g);
            pR.push_back(j); pC.push_back(j + 1); pV.push_back(v);
            pR.push_back(j + 1); pC.push_back(j); pV.push_back(v);
        }
    }

    std::vector<double> aR, aC, aV;
    aR.reserve(static_cast<std::size_t>(n) * nnzPerRow);
    aC.reserve(static_cast<std::size_t>(n) * nnzPerRow);
    aV.reserve(static_cast<std::size_t>(n) * nnzPerRow);
    for (int i = 0; i < n; ++i) {
        const int first = std::min(std::max(i - band / 2, 0), n - band);
        for (int k = 0; k < nnzPerRow; ++k) {
            double v = val(g);
            if (v == 0.0) v = 1.0;
            aR.push_back(i); aC.push_back(first + offset(g)); aV.push_back(v);
        }
    }

    qp::QpModel m;
    m.P = mkMat(n, n, pR, pC, pV);
    m.A = mkMat(n, n, aR, aC, aV);
    m.q.resize(static_cast<std::size_t>(n));
    std::vector<double> reference(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) {
        m.q[j] = val(g);
        reference[j] = val(g);
    }
    std::vector<double> activity(static_cast<std::size_t>(n), 0.0);
    for (std::size_t k = 0; k < aR.size(); ++k) {
        activity[static_cast<std::size_t>(aR[k])] += aV[k] * reference[static_cast<std::size_t>(aC[k])];
    }
    const double inf = std::numeric_limits<double>::infinity();
    m.l.resize(static_cast<std::size_t>(n));
    m.u.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        switch (i % 3) {
            case 0: m.l[i] = activity[i]; m.u[i] = activity[i]; break;
            case 1: m.l[i] = -inf; m.u[i] = activity[i] + 1.0; break;
            default: m.l[i] = activity[i] - 1.0; m.u[i] = activity[i] + 1.0; break;
        }
    }
    return m;
}

std::string jsonNumber(double value) {
    if (!std::isfinite(value)) {
        return "null";
    }
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    return buffer;
}

// qp_bench sweep <cpu|cuda|auto> <n> [nnzPerRow] [repeats] [threads]
//
// One untimed warm-up solve, then `repeats` timed end-to-end QpSolver::solve
// calls on bandedSweepModel; the median is reported. Wall time is measured
// with steady_clock (the MSVC and Linux builds this mode is meant for have a
// reliable one). The last line is "json " plus one JSON object.
int sweep(int argc, char** argv) {
    if (argc < 4) {
        std::printf("usage: qp_bench sweep <cpu|cuda|auto> <n> [nnzPerRow] [repeats] [threads]\n");
        return 2;
    }
    qp::ComputeBackend backend = qp::ComputeBackend::Cpu;
    if (std::strcmp(argv[2], "cuda") == 0) {
        backend = qp::ComputeBackend::Cuda;
    } else if (std::strcmp(argv[2], "auto") == 0) {
        backend = qp::ComputeBackend::Auto;
    } else if (std::strcmp(argv[2], "cpu") != 0) {
        std::printf("unknown backend '%s' (cpu|cuda|auto)\n", argv[2]);
        return 2;
    }
    const int n = std::atoi(argv[3]);
    const int nnzPerRow = argc > 4 ? std::atoi(argv[4]) : 8;
    const int repeats = std::max(argc > 5 ? std::atoi(argv[5]) : 5, 1);
    const int threads = argc > 6 ? std::atoi(argv[6]) : 0;
    if (n < 2 || nnzPerRow < 1) {
        std::printf("n must be >= 2 and nnzPerRow >= 1\n");
        return 2;
    }

    const qp::QpModel model = bandedSweepModel(n, nnzPerRow, 2024u);
    qp::AdmmOptions options;
    options.backend = backend;
    options.threadCount = threads;

    struct Timed {
        qp::AdmmResult result;
        double seconds = 0.0;
    };
    auto timedSolve = [&] {
        const auto start = std::chrono::steady_clock::now();
        Timed t;
        t.result = qp::QpSolver{}.solve(model, options);
        t.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        return t;
    };

    const Timed warmUp = timedSolve();
    if (warmUp.result.status == qp::QpStatus::InvalidProblem) {
        std::printf("solve refused: %s\n", warmUp.result.backendMessage.c_str());
        return 1;
    }
    std::vector<Timed> runs;
    for (int r = 0; r < repeats; ++r) {
        runs.push_back(timedSolve());
    }
    std::vector<double> runSeconds;
    for (const Timed& t : runs) runSeconds.push_back(t.seconds);
    std::sort(runs.begin(), runs.end(), [](const Timed& a, const Timed& b) { return a.seconds < b.seconds; });
    const Timed& median = runs[runs.size() / 2];
    const qp::AdmmResult& r = median.result;
    const qp::CudaAvailability cuda = qp::cudaAvailability(0);
    const long long nnz = static_cast<long long>(model.P.nonzeros() + model.A.nonzeros());

    std::printf("n=%d m=%d nnz(P)=%zu nnz(A)=%zu backend=%s (%s)\n", n, model.numConstraints(),
                model.P.nonzeros(), model.A.nonzeros(), qp::toString(r.executedBackend),
                r.backendMessage.c_str());
    std::printf("status=%s iters=%lld obj=%.12g rPri=%.3e rDual=%.3e warm-up=%.4f s median=%.4f s "
                "(%d runs) kkt solve=%.4f s kkt factor=%.4f s setup=%.4f s H2D=%lld B D2H=%lld B\n",
                qp::toString(r.status), static_cast<long long>(r.iterations), r.primalObjective,
                r.primalResidual, r.dualResidual, warmUp.seconds, median.seconds, repeats,
                r.kktSolveSeconds, r.kktFactorSeconds, r.backendProfile.setupSeconds,
                static_cast<long long>(r.backendProfile.hostToDeviceBytes),
                static_cast<long long>(r.backendProfile.deviceToHostBytes));

    std::string list;
    for (std::size_t k = 0; k < runSeconds.size(); ++k) {
        list += (k == 0 ? "" : ",") + jsonNumber(runSeconds[k]);
    }
    std::printf("json {\"engine\":\"qp\",\"rows\":%d,\"columns\":%d,\"nonzeros\":%lld,"
                "\"nnz_p\":%zu,\"nnz_a\":%zu,\"threads\":%d,"
                "\"requested_backend\":\"%s\",\"executed_backend\":\"%s\",\"cuda_device\":\"%s\","
                "\"cuda_runtime\":%d,\"status\":\"%s\",\"iterations\":%lld,"
                "\"warmup_seconds\":%s,\"run_seconds\":[%s],\"median_seconds\":%s,"
                "\"kkt_solve_seconds\":%s,\"kkt_factor_seconds\":%s,\"setup_seconds\":%s,"
                "\"snapshot_seconds\":%s,\"h2d_bytes\":%lld,\"d2h_bytes\":%lld,\"host_syncs\":%lld,"
                "\"objective\":%s,\"primal_residual\":%s,\"dual_residual\":%s}\n",
                model.numConstraints(), n, nnz, model.P.nonzeros(), model.A.nonzeros(), threads,
                qp::toString(backend), qp::toString(r.executedBackend),
                cuda.usable ? cuda.deviceName.c_str() : "", cuda.runtimeVersion,
                qp::toString(r.status), static_cast<long long>(r.iterations),
                jsonNumber(warmUp.seconds).c_str(), list.c_str(), jsonNumber(median.seconds).c_str(),
                jsonNumber(r.kktSolveSeconds).c_str(), jsonNumber(r.kktFactorSeconds).c_str(),
                jsonNumber(r.backendProfile.setupSeconds).c_str(),
                jsonNumber(r.backendProfile.snapshotSeconds).c_str(),
                static_cast<long long>(r.backendProfile.hostToDeviceBytes),
                static_cast<long long>(r.backendProfile.deviceToHostBytes),
                static_cast<long long>(r.backendProfile.synchronisations),
                jsonNumber(r.primalObjective).c_str(), jsonNumber(r.primalResidual).c_str(),
                jsonNumber(r.dualResidual).c_str());
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc > 1 && std::strcmp(argv[1], "sweep") == 0) {
        return sweep(argc, argv);
    }
    if (argc > 1) {
        if (std::strcmp(argv[1], "cuda") == 0) {
            gBackend = qp::ComputeBackend::Cuda;
        } else if (std::strcmp(argv[1], "auto") == 0) {
            gBackend = qp::ComputeBackend::Auto;
        } else if (std::strcmp(argv[1], "cpu") != 0) {
            std::printf("unknown backend '%s' (cpu|cuda|auto)\n", argv[1]);
            return 2;
        }
    }
    const qp::CudaAvailability cuda = qp::cudaAvailability(0);
    std::printf("backend=%s  cuda compiled=%s  device=%s  runtime=%d driver=%d\n",
                qp::toString(gBackend), cuda.compiled ? "yes" : "no",
                cuda.usable ? cuda.deviceName.c_str() : "none", cuda.runtimeVersion, cuda.driverVersion);

    section("small unconstrained (n=20)");
    {
        BenchCase bc{"", smallUnconstrained(), qp::AdmmOptions{}};
        bc.options.rho = 1.0;
        bc.options.useRuizScaling = false;
        run("n=20, rho=1.0", bc);
    }

    section("medium box-constrained");
    {
        for (int n : {50, 200}) {
            BenchCase bc{"", mediumBoxed(n), qp::AdmmOptions{}};
            bc.options.rho = 5.0;
            bc.options.useRuizScaling = true;
            bc.options.iterationLimit = 5000;
            char label[64];
            std::snprintf(label, sizeof(label), "n=%d, rho=5.0, Ruiz", n);
            run(label, bc);
        }
    }

    section("ill-conditioned ranged (Ruiz on)");
    {
        for (int n : {50, 200}) {
            BenchCase bc{"", illConditionedRanged(n), qp::AdmmOptions{}};
            bc.options.rho = 1.0;
            bc.options.useRuizScaling = true;
            bc.options.ruizIterations = 10;
            bc.options.iterationLimit = 5000;
            char label[64];
            std::snprintf(label, sizeof(label), "n=%d, Ruiz=10", n);
            run(label, bc);
        }
    }

    section("large sparse (sparse Cholesky path)");
    {
        for (int n : {600, 1000}) {
            BenchCase bc{"", largeSparse(n, 0.01), qp::AdmmOptions{}};
            bc.options.rho = 1.0;
            bc.options.useRuizScaling = true;
            bc.options.iterationLimit = 3000;
            char label[64];
            std::snprintf(label, sizeof(label), "n=%d, density=0.01", n);
            run(label, bc);
        }
    }

    std::printf("\n");
    return 0;
}
