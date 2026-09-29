#include "qp/super_admm_solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qp {
namespace {

/*
 * Validate that P is symmetric positive semidefinite.
 *
 * SuperADMM is intended for convex QPs. The repository's generic
 * QpModel::validate() checks dimensions and bounds, but does not
 * currently enforce symmetry or positive semidefiniteness of P.
 *
 * This implementation therefore performs that validation locally.
 *
 * A dense Jacobi eigenvalue computation is used because the current
 * SuperADMM KKT backend is itself a dense reference implementation.
 */
bool validatePositiveSemidefiniteP(
    const QpModel& model,
    std::string& message
) {
    const int n = model.numVariables();

    if (n == 0) {
        return true;
    }

    std::vector<double> dense(
        static_cast<std::size_t>(n) * n,
        0.0
    );

    const auto& rowStart = model.P.csrRowStart();
    const auto& columnIndex = model.P.csrColumnIndex();
    const auto& values = model.P.csrValues();

    for (int row = 0; row < n; ++row) {
        const int begin =
            rowStart[static_cast<std::size_t>(row)];

        const int end =
            rowStart[static_cast<std::size_t>(row + 1)];

        for (int k = begin; k < end; ++k) {
            const int col =
                columnIndex[static_cast<std::size_t>(k)];

            if (col < 0 || col >= n) {
                message =
                    "SuperADMM quadratic matrix P has an invalid index";
                return false;
            }

            const double value =
                values[static_cast<std::size_t>(k)];

            if (!std::isfinite(value)) {
                message =
                    "SuperADMM quadratic matrix P contains a non-finite value";
                return false;
            }

            dense[
                static_cast<std::size_t>(row) * n + col
            ] += value;
        }
    }

    /*
     * Check symmetry.
     *
     * The problem is represented using a full symmetric P matrix.
     * We do not silently symmetrize P because that would change the
     * optimization problem supplied by the user.
     */
    double scale = 1.0;

    for (double value : dense) {
        scale = std::max(scale, std::abs(value));
    }

    const double symmetryTolerance =
        1e-10 * scale;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double pij =
                dense[
                    static_cast<std::size_t>(i) * n + j
                ];

            const double pji =
                dense[
                    static_cast<std::size_t>(j) * n + i
                ];

            if (std::abs(pij - pji) >
                symmetryTolerance) {

                message =
                    "SuperADMM requires a symmetric quadratic matrix P";

                return false;
            }
        }
    }

    /*
     * Jacobi diagonalization for a real symmetric matrix.
     *
     * After convergence, the diagonal entries approximate the
     * eigenvalues of P.
     */
    const int maxSweeps =
        std::max(50, 5 * n * n);

    const double eigenTolerance =
        1e-10 * scale;

    for (int sweep = 0; sweep < maxSweeps; ++sweep) {
        double maxOffDiagonal = 0.0;
        int p = -1;
        int q = -1;

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double value =
                    std::abs(
                        dense[
                            static_cast<std::size_t>(i) * n + j
                        ]
                    );

                if (value > maxOffDiagonal) {
                    maxOffDiagonal = value;
                    p = i;
                    q = j;
                }
            }
        }

        if (maxOffDiagonal <= eigenTolerance) {
            break;
        }

        if (p < 0 || q < 0) {
            break;
        }

        const double app =
            dense[
                static_cast<std::size_t>(p) * n + p
            ];

        const double aqq =
            dense[
                static_cast<std::size_t>(q) * n + q
            ];

        const double apq =
            dense[
                static_cast<std::size_t>(p) * n + q
            ];

        if (std::abs(apq) <= eigenTolerance) {
            continue;
        }

        /*
         * Stable Jacobi rotation.
         */
        const double theta =
            0.5 * std::atan2(
                2.0 * apq,
                aqq - app
            );

        const double c = std::cos(theta);
        const double s = std::sin(theta);

        /*
         * Rotate rows and columns p and q.
         */
        for (int k = 0; k < n; ++k) {
            if (k == p || k == q) {
                continue;
            }

            const double akp =
                dense[
                    static_cast<std::size_t>(k) * n + p
                ];

            const double akq =
                dense[
                    static_cast<std::size_t>(k) * n + q
                ];

            const double newKp =
                c * akp - s * akq;

            const double newKq =
                s * akp + c * akq;

            dense[
                static_cast<std::size_t>(k) * n + p
            ] = newKp;

            dense[
                static_cast<std::size_t>(p) * n + k
            ] = newKp;

            dense[
                static_cast<std::size_t>(k) * n + q
            ] = newKq;

            dense[
                static_cast<std::size_t>(q) * n + k
            ] = newKq;
        }

        const double newApp =
            c * c * app -
            2.0 * s * c * apq +
            s * s * aqq;

        const double newAqq =
            s * s * app +
            2.0 * s * c * apq +
            c * c * aqq;

        dense[
            static_cast<std::size_t>(p) * n + p
        ] = newApp;

        dense[
            static_cast<std::size_t>(q) * n + q
        ] = newAqq;

        dense[
            static_cast<std::size_t>(p) * n + q
        ] = 0.0;

        dense[
            static_cast<std::size_t>(q) * n + p
        ] = 0.0;
    }

    /*
     * Check whether the Jacobi iteration converged sufficiently.
     */
    double finalOffDiagonal = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            finalOffDiagonal =
                std::max(
                    finalOffDiagonal,
                    std::abs(
                        dense[
                            static_cast<std::size_t>(i) * n + j
                        ]
                    )
                );
        }
    }

    if (finalOffDiagonal >
        100.0 * eigenTolerance) {

        message =
            "SuperADMM could not reliably determine whether P is PSD";

        return false;
    }

    /*
     * A symmetric matrix is PSD iff every eigenvalue is >= 0.
     *
     * A small negative tolerance is allowed for floating-point
     * roundoff.
     */
    const double psdTolerance =
        1e-9 * scale;

    for (int i = 0; i < n; ++i) {
        const double eigenvalue =
            dense[
                static_cast<std::size_t>(i) * n + i
            ];

        if (!std::isfinite(eigenvalue)) {
            message =
                "SuperADMM encountered a non-finite eigenvalue of P";
            return false;
        }

        if (eigenvalue < -psdTolerance) {
            message =
                "SuperADMM requires P to be positive semidefinite";
            return false;
        }
    }

    return true;
}

}  // namespace

SuperAdmmSolver::SuperAdmmSolver(
    const QpModel& model,
    const SuperAdmmOptions& options
)
    : model_(model),
      options_(options),
      kktSolver_(model) {

    model_.validate();

    if (options_.iterationLimit < 0) {
        throw std::invalid_argument(
            "SuperADMM iterationLimit must be non-negative"
        );
    }

    if (options_.timeLimitSeconds < 0.0 ||
        !std::isfinite(options_.timeLimitSeconds)) {
        throw std::invalid_argument(
            "SuperADMM timeLimitSeconds must be non-negative and finite"
        );
    }

    if (options_.primalTolerance <= 0.0 ||
        !std::isfinite(options_.primalTolerance)) {
        throw std::invalid_argument(
            "SuperADMM primalTolerance must be positive"
        );
    }

    if (options_.dualTolerance <= 0.0 ||
        !std::isfinite(options_.dualTolerance)) {
        throw std::invalid_argument(
            "SuperADMM dualTolerance must be positive"
        );
    }

    if (options_.alpha <= 1.0 ||
        !std::isfinite(options_.alpha)) {
        throw std::invalid_argument(
            "SuperADMM alpha must be greater than 1"
        );
    }

    if (options_.sigma <= 0.0 ||
        !std::isfinite(options_.sigma)) {
        throw std::invalid_argument(
            "SuperADMM sigma must be positive"
        );
    }

    if (options_.b0 < 1.0 ||
        !std::isfinite(options_.b0)) {
        throw std::invalid_argument(
            "SuperADMM b0 must be at least 1"
        );
    }

    if (options_.tau <= 0.0 ||
        options_.tau >= 1.0 ||
        !std::isfinite(options_.tau)) {
        throw std::invalid_argument(
            "SuperADMM tau must be in (0,1)"
        );
    }

    if (options_.rho0 <= 0.0 ||
        !std::isfinite(options_.rho0)) {
        throw std::invalid_argument(
            "SuperADMM rho0 must be positive"
        );
    }

    if (options_.infeasibilityCheckInterval <= 0) {
        throw std::invalid_argument(
            "SuperADMM infeasibilityCheckInterval must be positive"
        );
    }

    if (options_.infeasibilityTolerance <= 0.0 ||
        !std::isfinite(options_.infeasibilityTolerance)) {
        throw std::invalid_argument(
            "SuperADMM infeasibilityTolerance must be positive"
        );
    }

    result_ = AdmmResult{};
    result_.status = QpStatus::IterationLimit;
    result_.statusMessage =
        "SuperADMM iteration limit reached";
}

bool SuperAdmmSolver::initialize() {
    const int n = model_.numVariables();
    const int m = model_.numConstraints();

    x_.assign(
        static_cast<std::size_t>(n),
        0.0
    );

    z_.assign(
        static_cast<std::size_t>(m),
        0.0
    );

    y_.assign(
        static_cast<std::size_t>(m),
        0.0
    );

    nu_.assign(
        static_cast<std::size_t>(m),
        0.0
    );

    zTilde_.assign(
        static_cast<std::size_t>(m),
        0.0
    );

    rho_.assign(
        static_cast<std::size_t>(m),
        options_.rho0
    );

    b_ = options_.b0;

    previousX_.clear();
    previousY_.clear();

    havePreviousIterate_ = false;

    result_.primal.assign(
        static_cast<std::size_t>(n),
        0.0
    );

    result_.constraintDual.assign(
        static_cast<std::size_t>(m),
        0.0
    );

    /*
     * z^0 = A x^0.
     *
     * Since x^0 = 0, this is the zero vector.
     */
    if (m > 0) {
        model_.A.multiply(
            x_,
            z_
        );
    }

    result_.primalResidual =
        primalResidual();

    result_.dualResidual =
        dualResidual();

    result_.primalObjective =
        objective(x_);

    result_.bestObjective =
        result_.primalObjective;

    result_.finalRho =
        rho_.empty()
            ? 0.0
            : rho_.front();

    result_.iterations = 0;

    return true;
}

AdmmResult SuperAdmmSolver::solve() {
    const auto start =
        std::chrono::steady_clock::now();

    try {
        model_.validate();
    } catch (const std::exception& e) {
        result_.status =
            QpStatus::InvalidProblem;

        result_.statusMessage =
            e.what();

        return result_;
    }

    /*
     * SuperADMM is a convex-QP algorithm.
     *
     * Reject an asymmetric or indefinite P instead of silently
     * solving a different problem.
     */
    {
        std::string convexityMessage;

        if (!validatePositiveSemidefiniteP(
                model_,
                convexityMessage)) {

            result_.status =
                QpStatus::InvalidProblem;

            result_.statusMessage =
                convexityMessage;

            return result_;
        }
    }

    if (!initialize()) {
        result_.status =
            QpStatus::NumericalFailure;

        result_.statusMessage =
            "SuperADMM initialization failed";

        return result_;
    }

    const int n =
        model_.numVariables();

    const int m =
        model_.numConstraints();

    /*
     * ------------------------------------------------------------
     * Zero-variable QP
     * ------------------------------------------------------------
     *
     * If n = 0, then Ax = 0 for every possible x.
     *
     * Therefore the problem is feasible iff
     *
     *     l_i <= 0 <= u_i
     *
     * for every constraint.
     */
    if (n == 0) {
        for (int i = 0; i < m; ++i) {
            const std::size_t index =
                static_cast<std::size_t>(i);

            const double lower =
                model_.l[index];

            const double upper =
                model_.u[index];

            if ((std::isfinite(lower) &&
                 lower > 0.0) ||
                (std::isfinite(upper) &&
                 upper < 0.0)) {

                result_.status =
                    QpStatus::Infeasible;

                result_.statusMessage =
                    "SuperADMM zero-variable problem is infeasible";

                result_.iterations = 0;

                result_.primal.clear();

                result_.constraintDual.assign(
                    static_cast<std::size_t>(m),
                    0.0
                );

                result_.primalObjective = 0.0;
                result_.bestObjective = 0.0;
                result_.primalResidual = 0.0;
                result_.dualResidual = 0.0;

                result_.dualObjective =
                    -std::numeric_limits<double>::infinity();

                result_.solveTimeSeconds =
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() -
                        start
                    ).count();

                return result_;
            }
        }

        result_.status =
            QpStatus::Optimal;

        result_.statusMessage =
            "SuperADMM solved zero-variable problem";

        result_.iterations = 0;

        result_.primal.clear();

        result_.constraintDual.assign(
            static_cast<std::size_t>(m),
            0.0
        );

        result_.primalObjective = 0.0;
        result_.bestObjective = 0.0;
        result_.primalResidual = 0.0;
        result_.dualResidual = 0.0;

        result_.dualObjective =
            -std::numeric_limits<double>::infinity();

        result_.solveTimeSeconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() -
                start
            ).count();

        return result_;
    }

    /*
     * If no iterations are allowed, return immediately.
     */
    if (options_.iterationLimit == 0) {
        result_.status =
            QpStatus::IterationLimit;

        result_.statusMessage =
            "SuperADMM iteration limit reached";

        result_.solveTimeSeconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() -
                start
            ).count();

        return result_;
    }

    /*
     * ------------------------------------------------------------
     * Main SuperADMM loop
     * ------------------------------------------------------------
     */
    for (std::int64_t iterationCount = 1;
         iterationCount <= options_.iterationLimit;
         ++iterationCount) {

        /*
         * Time-limit check before starting the next KKT solve.
         */
        if (options_.timeLimitSeconds > 0.0) {
            const double elapsed =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() -
                    start
                ).count();

            if (elapsed >= options_.timeLimitSeconds) {
                result_.status =
                    QpStatus::TimeLimit;

                result_.statusMessage =
                    "SuperADMM time limit reached";

                break;
            }
        }

        /*
         * Perform one complete SuperADMM iteration.
         */
        if (!iteration()) {
            result_.status =
                QpStatus::NumericalFailure;

            result_.statusMessage =
                "SuperADMM KKT solve or iteration failed";

            break;
        }

        result_.iterations =
            iterationCount;

        result_.primal =
            x_;

        result_.constraintDual =
            y_;

        result_.primalResidual =
            primalResidual();

        result_.dualResidual =
            dualResidual();

        result_.primalObjective =
            objective(x_);

        result_.bestObjective =
            result_.primalObjective;

        if (!std::isfinite(result_.primalObjective) ||
            !std::isfinite(result_.primalResidual) ||
            !std::isfinite(result_.dualResidual)) {

            result_.status =
                QpStatus::NumericalFailure;

            result_.statusMessage =
                "SuperADMM produced a non-finite iterate";

            break;
        }

        result_.finalRho =
            rho_.empty()
                ? 0.0
                : *std::max_element(
                      rho_.begin(),
                      rho_.end()
                  );

        /*
         * Optimality test.
         *
         * r_prim = ||Ax - z||_inf
         *
         * r_dual = ||Px + q + A^T y||_inf
         */
        if (result_.primalResidual <=
                options_.primalTolerance &&
            result_.dualResidual <=
                options_.dualTolerance) {

            result_.status =
                QpStatus::Optimal;

            result_.statusMessage =
                "SuperADMM converged";

            break;
        }

        /*
         * Infeasibility checks.
         */
        if (iterationCount %
                options_.infeasibilityCheckInterval ==
            0) {

            if (checkPrimalInfeasibility()) {
                result_.status =
                    QpStatus::Infeasible;

                result_.statusMessage =
                    "SuperADMM primal infeasibility detected";

                break;
            }

            if (checkDualInfeasibility()) {
                result_.status =
                    QpStatus::Unbounded;

                result_.statusMessage =
                    "SuperADMM dual infeasibility detected";

                break;
            }
        }

        /*
         * The numerical-stability requirement is
         *
         *     1 / b <= rho_i <= b.
         *
         * Once b < 1 this condition cannot be satisfied.
         */
        if (b_ < 1.0) {
            result_.status =
                QpStatus::NumericalFailure;

            result_.statusMessage =
                "SuperADMM stability bound fell below 1";

            break;
        }
    }

    /*
     * Store the final state.
     */
    result_.primal =
        x_;

    result_.constraintDual =
        y_;

    result_.primalResidual =
        primalResidual();

    result_.dualResidual =
        dualResidual();

    result_.primalObjective =
        objective(x_);

    result_.bestObjective =
        result_.primalObjective;

    result_.finalRho =
        rho_.empty()
            ? 0.0
            : *std::max_element(
                  rho_.begin(),
                  rho_.end()
              );

    /*
     * This implementation does not construct a separate dual
     * objective. Do not fabricate one.
     */
    result_.dualObjective =
        -std::numeric_limits<double>::infinity();

    result_.solveTimeSeconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() -
            start
        ).count();

    return result_;
}

bool SuperAdmmSolver::iteration() {
    const int n =
        model_.numVariables();

    const int m =
        model_.numConstraints();

    /*
     * Save the current k-th iterate.
     */
    const std::vector<double> oldX =
        x_;

    const std::vector<double> oldZ =
        z_;

    const std::vector<double> oldY =
        y_;

    /*
     * Save the previous iterate for the infeasibility
     * certificates.
     */
    if (havePreviousIterate_) {
        previousX_ =
            oldX;

        previousY_ =
            oldY;
    }

    /*
     * ------------------------------------------------------------
     * KKT system
     * ------------------------------------------------------------
     *
     * [ P + sigma I     A^T      ] [x^{k+1} ] =
     * [ A              -R^{-1}  ] [nu^{k+1}]
     *
     * RHS:
     *
     * [ sigma x^k - q              ]
     * [ z^k - R^{-1} y^k          ]
     */
    std::vector<double> rhsX(
        static_cast<std::size_t>(n),
        0.0
    );

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        rhsX[index] =
            options_.sigma *
                oldX[index] -
            model_.q[index];
    }

    std::vector<double> rhsNu(
        static_cast<std::size_t>(m),
        0.0
    );

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double ri =
            rho_[index];

        if (!std::isfinite(ri) ||
            ri <= 0.0) {
            return false;
        }

        rhsNu[index] =
            oldZ[index] -
            oldY[index] / ri;
    }

    std::vector<double> newX;
    std::vector<double> newNu;

    if (!kktSolver_.solve(
            rho_,
            options_.sigma,
            rhsX,
            rhsNu,
            newX,
            newNu)) {
        return false;
    }

    if (!finiteVector(newX) ||
        !finiteVector(newNu)) {
        return false;
    }

    x_ =
        std::move(newX);

    nu_ =
        std::move(newNu);

    /*
     * ------------------------------------------------------------
     * z-tilde update
     * ------------------------------------------------------------
     *
     * ztilde^{k+1}
     * =
     * z^k + R^{-1}(nu^{k+1} - y^k)
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        zTilde_[index] =
            oldZ[index] +
            (
                nu_[index] -
                oldY[index]
            ) / rho_[index];
    }

    /*
     * ------------------------------------------------------------
     * Projection
     * ------------------------------------------------------------
     *
     * z^{k+1}
     * =
     * Pi_[l,u](
     *     ztilde^{k+1} + R^{-1}y^k
     * )
     */
    project(
        zTilde_,
        oldY,
        z_
    );

    /*
     * ------------------------------------------------------------
     * Dual update
     * ------------------------------------------------------------
     *
     * y^{k+1}
     * =
     * y^k + R(
     *     ztilde^{k+1} - z^{k+1}
     * )
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        y_[index] =
            oldY[index] +
            rho_[index] *
                (
                    zTilde_[index] -
                    z_[index]
                );
    }

    if (!finiteVector(x_) ||
        !finiteVector(z_) ||
        !finiteVector(y_) ||
        !finiteVector(nu_)) {
        return false;
    }

    /*
     * KKT residual epsilon.
     */
    const double epsilon =
        kktResidual(
            oldX,
            oldZ,
            oldY
        );

    if (!std::isfinite(epsilon)) {
        return false;
    }

    const double rPrim =
        primalResidual();

    if (!std::isfinite(rPrim)) {
        return false;
    }

    /*
     * ------------------------------------------------------------
     * Stability-bound update
     * ------------------------------------------------------------
     *
     * b^{k+1} =
     *
     *     tau b^k,   if epsilon >= r_prim
     *
     *     b^k,       otherwise.
     */
    updateBound(
        epsilon,
        rPrim
    );

    /*
     * IMPORTANT:
     *
     * The SuperADMM algorithm updates R using z^k,
     * i.e. the OLD projected variable.
     */
    if (b_ < 1.0) {
        /*
         * Do not construct a new R when the stability condition
         *
         *     1/b <= R_i <= b
         *
         * has become impossible.
         *
         * solve() will report NumericalFailure after this
         * iteration.
         */
        havePreviousIterate_ =
            true;

        return true;
    }

    updateWeights(
        oldZ
    );

    havePreviousIterate_ =
        true;

    return true;
}

void SuperAdmmSolver::project(
    const std::vector<double>& zTilde,
    const std::vector<double>& yOld,
    std::vector<double>& zNew
) const {
    const int m =
        model_.numConstraints();

    zNew.resize(
        static_cast<std::size_t>(m)
    );

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double value =
            zTilde[index] +
            yOld[index] / rho_[index];

        const double lower =
            model_.l[index];

        const double upper =
            model_.u[index];

        double projected =
            value;

        if (std::isfinite(lower) &&
            projected < lower) {
            projected =
                lower;
        }

        if (std::isfinite(upper) &&
            projected > upper) {
            projected =
                upper;
        }

        zNew[index] =
            projected;
    }
}

void SuperAdmmSolver::updateBound(
    double epsilon,
    double primalResidualValue
) {
    if (!std::isfinite(epsilon) ||
        !std::isfinite(primalResidualValue)) {
        return;
    }

    if (epsilon >= primalResidualValue) {
        b_ *= options_.tau;
    }
}

void SuperAdmmSolver::updateWeights(
    const std::vector<double>& oldZ
) {
    const int m =
        model_.numConstraints();

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double lower =
            model_.l[index];

        const double upper =
            model_.u[index];

        const double zi =
            oldZ[index];

        /*
         * A constraint is considered active when z_i^k is
         * exactly on a finite lower or upper bound.
         */
        const bool atLower =
            std::isfinite(lower) &&
            zi == lower;

        const bool atUpper =
            std::isfinite(upper) &&
            zi == upper;

        const bool active =
            atLower ||
            atUpper;

        double newRho =
            0.0;

        if (active) {
            /*
             * Active constraint:
             *
             * R_i^{k+1}
             * =
             * min(b^{k+1}, alpha R_i^k)
             */
            newRho =
                std::min(
                    b_,
                    options_.alpha *
                        rho_[index]
                );
        } else {
            /*
             * Inactive constraint:
             *
             * R_i^{k+1}
             * =
             * max(1/b^{k+1}, R_i^k / alpha)
             */
            newRho =
                std::max(
                    1.0 / b_,
                    rho_[index] /
                        options_.alpha
                );
        }

        if (!std::isfinite(newRho) ||
            newRho <= 0.0) {
            return;
        }

        /*
         * Explicitly enforce the numerical stability interval.
         */
        newRho =
            std::max(
                newRho,
                1.0 / b_
            );

        newRho =
            std::min(
                newRho,
                b_
            );

        rho_[index] =
            newRho;
    }
}

double SuperAdmmSolver::primalResidual() const {
    const int m =
        model_.numConstraints();

    if (m == 0) {
        return 0.0;
    }

    std::vector<double> Ax;

    model_.A.multiply(
        x_,
        Ax
    );

    double residual =
        0.0;

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        residual =
            std::max(
                residual,
                std::abs(
                    Ax[index] -
                    z_[index]
                )
            );
    }

    return residual;
}

double SuperAdmmSolver::dualResidual() const {
    const int n =
        model_.numVariables();

    std::vector<double> Px;
    std::vector<double> Aty;

    model_.P.multiply(
        x_,
        Px
    );

    model_.A.transposeMultiply(
        y_,
        Aty
    );

    double residual =
        0.0;

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        const double value =
            Px[index] +
            model_.q[index] +
            Aty[index];

        residual =
            std::max(
                residual,
                std::abs(value)
            );
    }

    return residual;
}

double SuperAdmmSolver::kktResidual(
    const std::vector<double>& oldX,
    const std::vector<double>& oldZ,
    const std::vector<double>& oldY
) const {
    const int n =
        model_.numVariables();

    const int m =
        model_.numConstraints();

    std::vector<double> Px;
    std::vector<double> Atnu;
    std::vector<double> Ax;

    model_.P.multiply(
        x_,
        Px
    );

    model_.A.transposeMultiply(
        nu_,
        Atnu
    );

    model_.A.multiply(
        x_,
        Ax
    );

    double epsilon =
        0.0;

    /*
     * Top KKT block:
     *
     * (P + sigma I)x^{k+1}
     * + A^T nu^{k+1}
     *
     * =
     *
     * sigma x^k - q
     */
    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        const double lhs =
            Px[index] +
            options_.sigma *
                x_[index] +
            Atnu[index];

        const double rhs =
            options_.sigma *
                oldX[index] -
            model_.q[index];

        epsilon =
            std::max(
                epsilon,
                std::abs(lhs - rhs)
            );
    }

    /*
     * Bottom KKT block:
     *
     * A x^{k+1}
     * - R^{-1} nu^{k+1}
     *
     * =
     *
     * z^k - R^{-1} y^k
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double ri =
            rho_[index];

        if (!std::isfinite(ri) ||
            ri <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }

        const double lhs =
            Ax[index] -
            nu_[index] / ri;

        const double rhs =
            oldZ[index] -
            oldY[index] / ri;

        epsilon =
            std::max(
                epsilon,
                std::abs(lhs - rhs)
            );
    }

    return epsilon;
}

double SuperAdmmSolver::objective(
    const std::vector<double>& x
) const {
    const int n =
        model_.numVariables();

    std::vector<double> Px;

    model_.P.multiply(
        x,
        Px
    );

    double value =
        0.0;

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        value +=
            0.5 *
            x[index] *
            Px[index];

        value +=
            model_.q[index] *
            x[index];
    }

    return value;
}

bool SuperAdmmSolver::checkPrimalInfeasibility() {
    const int m =
        model_.numConstraints();

    if (!havePreviousIterate_ ||
        m == 0) {
        return false;
    }

    if (previousY_.size() !=
        y_.size()) {
        return false;
    }

    std::vector<double> deltaY(
        static_cast<std::size_t>(m),
        0.0
    );

    bool hasNonzeroDirection =
        false;

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        deltaY[index] =
            y_[index] -
            previousY_[index];

        if (std::abs(deltaY[index]) >
            options_.infeasibilityTolerance) {
            hasNonzeroDirection =
                true;
        }
    }

    if (!hasNonzeroDirection) {
        return false;
    }

    /*
     * Primal infeasibility certificate:
     *
     * A^T delta_y = 0
     *
     * and
     *
     * u^T delta_y_+
     * +
     * l^T delta_y_-
     * < 0
     */
    std::vector<double> ATdeltaY;

    model_.A.transposeMultiply(
        deltaY,
        ATdeltaY
    );

    for (double value : ATdeltaY) {
        if (std::abs(value) >
            options_.infeasibilityTolerance) {
            return false;
        }
    }

    double certificate =
        0.0;

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double dy =
            deltaY[index];

        if (dy >
            options_.infeasibilityTolerance) {

            const double upper =
                model_.u[index];

            if (!std::isfinite(upper)) {
                return false;
            }

            certificate +=
                upper * dy;

        } else if (
            dy <
            -options_.infeasibilityTolerance
        ) {

            const double lower =
                model_.l[index];

            if (!std::isfinite(lower)) {
                return false;
            }

            certificate +=
                lower * dy;
        }
    }

    return certificate <
           -options_.infeasibilityTolerance;
}

bool SuperAdmmSolver::checkDualInfeasibility() {
    const int n =
        model_.numVariables();

    const int m =
        model_.numConstraints();

    if (!havePreviousIterate_ ||
        n == 0) {
        return false;
    }

    if (previousX_.size() !=
        x_.size()) {
        return false;
    }

    std::vector<double> deltaX(
        static_cast<std::size_t>(n),
        0.0
    );

    bool hasNonzeroDirection =
        false;

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        deltaX[index] =
            x_[index] -
            previousX_[index];

        if (std::abs(deltaX[index]) >
            options_.infeasibilityTolerance) {
            hasNonzeroDirection =
                true;
        }
    }

    if (!hasNonzeroDirection) {
        return false;
    }

    /*
     * q^T delta_x < 0
     */
    double qDeltaX =
        0.0;

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        qDeltaX +=
            model_.q[index] *
            deltaX[index];
    }

    if (qDeltaX >=
        -options_.infeasibilityTolerance) {
        return false;
    }

    /*
     * P delta_x = 0
     */
    std::vector<double> PdeltaX;

    model_.P.multiply(
        deltaX,
        PdeltaX
    );

    for (double value : PdeltaX) {
        if (std::abs(value) >
            options_.infeasibilityTolerance) {
            return false;
        }
    }

    /*
     * Constraint-direction conditions:
     *
     * finite lower + finite upper:
     *
     *     A_i delta_x = 0
     *
     * finite lower + infinite upper:
     *
     *     A_i delta_x >= 0
     *
     * infinite lower + finite upper:
     *
     *     A_i delta_x <= 0
     *
     * infinite lower + infinite upper:
     *
     *     no restriction
     */
    if (m > 0) {
        std::vector<double> AdeltaX;

        model_.A.multiply(
            deltaX,
            AdeltaX
        );

        for (int i = 0; i < m; ++i) {
            const std::size_t index =
                static_cast<std::size_t>(i);

            const double a =
                AdeltaX[index];

            const bool lowerFinite =
                std::isfinite(
                    model_.l[index]
                );

            const bool upperFinite =
                std::isfinite(
                    model_.u[index]
                );

            if (lowerFinite &&
                upperFinite) {

                if (std::abs(a) >
                    options_.infeasibilityTolerance) {
                    return false;
                }

            } else if (
                lowerFinite &&
                !upperFinite
            ) {

                if (a <
                    -options_.infeasibilityTolerance) {
                    return false;
                }

            } else if (
                !lowerFinite &&
                upperFinite
            ) {

                if (a >
                    options_.infeasibilityTolerance) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool SuperAdmmSolver::finiteVector(
    const std::vector<double>& values
) const {
    for (double value : values) {
        if (!std::isfinite(value)) {
            return false;
        }
    }

    return true;
}

}  // namespace qp