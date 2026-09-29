#include "qp/super_admm_solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace qp {
namespace {

/*
 * Dense PSD validation.
 *
 * SuperADMM is defined for convex QPs, so P must be
 * symmetric positive semidefinite.
 *
 * This is deliberately a correctness-first implementation.
 * The production sparse KKT backend can replace this later.
 */
bool isPositiveSemidefinite(
    const SparseMatrix& matrix
) {
    const int rows = matrix.rows();
    const int columns = matrix.columns();

    if (rows != columns) {
        return false;
    }

    const int n = rows;

    if (n == 0) {
        return true;
    }

    std::vector<double> dense(
        static_cast<std::size_t>(n) *
            static_cast<std::size_t>(n),
        0.0
    );

    const auto& rowStart =
        matrix.csrRowStart();

    const auto& columnIndex =
        matrix.csrColumnIndex();

    const auto& values =
        matrix.csrValues();

    for (int i = 0; i < n; ++i) {
        const Offset begin =
            rowStart[
                static_cast<std::size_t>(i)
            ];

        const Offset end =
            rowStart[
                static_cast<std::size_t>(i + 1)
            ];

        for (Offset k = begin; k < end; ++k) {
            const int j =
                columnIndex[
                    static_cast<std::size_t>(k)
                ];

            if (j < 0 || j >= n) {
                return false;
            }

            const double value =
                values[
                    static_cast<std::size_t>(k)
                ];

            if (!std::isfinite(value)) {
                return false;
            }

            dense[
                static_cast<std::size_t>(i) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(j)
            ] += value;
        }
    }

    /*
     * Symmetry check.
     */
    double matrixScale = 0.0;

    for (double value : dense) {
        matrixScale =
            std::max(
                matrixScale,
                std::abs(value)
            );
    }

    const double symmetryTolerance =
        1e-10 *
        std::max(1.0, matrixScale);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double a =
                dense[
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j)
                ];

            const double b =
                dense[
                    static_cast<std::size_t>(j) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(i)
                ];

            if (std::abs(a - b) >
                symmetryTolerance) {
                return false;
            }
        }
    }

    /*
     * Jacobi eigenvalue iteration.
     *
     * A symmetric matrix is PSD iff all eigenvalues
     * are non-negative.
     */
    std::vector<double> a =
        std::move(dense);

    const int maxSweeps =
        std::max(
            10,
            5 * n * n
        );

    const double eigenTolerance =
        1e-10 *
        std::max(1.0, matrixScale);

    for (int sweep = 0;
         sweep < maxSweeps;
         ++sweep) {

        double largestOffDiagonal = 0.0;
        int p = 0;
        int q = 0;

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double value =
                    std::abs(
                        a[
                            static_cast<std::size_t>(i) *
                                static_cast<std::size_t>(n) +
                            static_cast<std::size_t>(j)
                        ]
                    );

                if (value > largestOffDiagonal) {
                    largestOffDiagonal = value;
                    p = i;
                    q = j;
                }
            }
        }

        if (largestOffDiagonal <= eigenTolerance) {
            break;
        }

        const double app =
            a[
                static_cast<std::size_t>(p) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(p)
            ];

        const double aqq =
            a[
                static_cast<std::size_t>(q) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(q)
            ];

        const double apq =
            a[
                static_cast<std::size_t>(p) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(q)
            ];

        if (apq == 0.0) {
            continue;
        }

        const double theta =
            0.5 *
            std::atan2(
                2.0 * apq,
                aqq - app
            );

        const double c =
            std::cos(theta);

        const double s =
            std::sin(theta);

        /*
         * Apply the Jacobi rotation.
         */
        for (int k = 0; k < n; ++k) {
            if (k == p || k == q) {
                continue;
            }

            const double akp =
                a[
                    static_cast<std::size_t>(k) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(p)
                ];

            const double akq =
                a[
                    static_cast<std::size_t>(k) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(q)
                ];

            const double newKp =
                c * akp - s * akq;

            const double newKq =
                s * akp + c * akq;

            a[
                static_cast<std::size_t>(k) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(p)
            ] = newKp;

            a[
                static_cast<std::size_t>(p) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(k)
            ] = newKp;

            a[
                static_cast<std::size_t>(k) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(q)
            ] = newKq;

            a[
                static_cast<std::size_t>(q) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(k)
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

        a[
            static_cast<std::size_t>(p) *
                static_cast<std::size_t>(n) +
            static_cast<std::size_t>(p)
        ] = newApp;

        a[
            static_cast<std::size_t>(q) *
                static_cast<std::size_t>(n) +
            static_cast<std::size_t>(q)
        ] = newAqq;

        a[
            static_cast<std::size_t>(p) *
                static_cast<std::size_t>(n) +
            static_cast<std::size_t>(q)
        ] = 0.0;

        a[
            static_cast<std::size_t>(q) *
                static_cast<std::size_t>(n) +
            static_cast<std::size_t>(p)
        ] = 0.0;
    }

    /*
     * Check the diagonalized matrix.
     */
    const double psdTolerance =
        1e-8 *
        std::max(1.0, matrixScale);

    for (int i = 0; i < n; ++i) {
        const double eigenvalue =
            a[
                static_cast<std::size_t>(i) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(i)
            ];

        if (!std::isfinite(eigenvalue)) {
            return false;
        }

        if (eigenvalue < -psdTolerance) {
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

    if (!std::isfinite(options_.timeLimitSeconds) ||
        options_.timeLimitSeconds < 0.0) {
        throw std::invalid_argument(
            "SuperADMM timeLimitSeconds must be non-negative"
        );
    }

    if (!std::isfinite(options_.primalTolerance) ||
        options_.primalTolerance <= 0.0) {
        throw std::invalid_argument(
            "SuperADMM primalTolerance must be positive"
        );
    }

    if (!std::isfinite(options_.dualTolerance) ||
        options_.dualTolerance <= 0.0) {
        throw std::invalid_argument(
            "SuperADMM dualTolerance must be positive"
        );
    }

    if (!std::isfinite(options_.alpha) ||
        options_.alpha <= 1.0) {
        throw std::invalid_argument(
            "SuperADMM alpha must be greater than 1"
        );
    }

    if (!std::isfinite(options_.sigma) ||
        options_.sigma <= 0.0) {
        throw std::invalid_argument(
            "SuperADMM sigma must be positive"
        );
    }

    if (!std::isfinite(options_.b0) ||
        options_.b0 < 1.0) {
        throw std::invalid_argument(
            "SuperADMM b0 must be at least 1"
        );
    }

    if (!std::isfinite(options_.tau) ||
        options_.tau <= 0.0 ||
        options_.tau >= 1.0) {
        throw std::invalid_argument(
            "SuperADMM tau must be in (0,1)"
        );
    }

    if (!std::isfinite(options_.rho0) ||
        options_.rho0 <= 0.0) {
        throw std::invalid_argument(
            "SuperADMM rho0 must be positive"
        );
    }

    /*
     * The numerical stability condition is
     *
     *     1 / b <= rho_i <= b.
     *
     * Therefore the initial rho must also satisfy it.
     */
    if (options_.rho0 < 1.0 / options_.b0 ||
        options_.rho0 > options_.b0) {
        throw std::invalid_argument(
            "SuperADMM rho0 must satisfy "
            "1/b0 <= rho0 <= b0"
        );
    }

    if (options_.infeasibilityCheckInterval <= 0) {
        throw std::invalid_argument(
            "SuperADMM infeasibilityCheckInterval "
            "must be positive"
        );
    }

    if (!std::isfinite(options_.infeasibilityTolerance) ||
        options_.infeasibilityTolerance <= 0.0) {
        throw std::invalid_argument(
            "SuperADMM infeasibilityTolerance "
            "must be positive"
        );
    }

    result_ = AdmmResult{};

    result_.status =
        QpStatus::IterationLimit;

    result_.statusMessage =
        "SuperADMM iteration limit reached";
}

bool SuperAdmmSolver::initialize() {
    const int n =
        model_.numVariables();

    const int m =
        model_.numConstraints();

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
     * Algorithm 1 initializes z^0 = A x^0.
     * Here x^0 = 0.
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

    return true;
}

bool SuperAdmmSolver::validateConvexObjective() const {
    return isPositiveSemidefinite(
        model_.P
    );
}

bool SuperAdmmSolver::validateZeroVariableProblem() const {
    const int m =
        model_.numConstraints();

    /*
     * With no optimization variables, Ax = 0.
     * Therefore feasibility is exactly
     *
     *     l_i <= 0 <= u_i.
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double lower =
            model_.l[index];

        const double upper =
            model_.u[index];

        if (0.0 < lower ||
            0.0 > upper) {
            return false;
        }
    }

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

    if (!validateConvexObjective()) {
        result_.status =
            QpStatus::InvalidProblem;

        result_.statusMessage =
            "SuperADMM requires P to be "
            "symmetric positive semidefinite";

        return result_;
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
     * Special case: zero-variable QP.
     */
    if (n == 0) {
        result_.solveTimeSeconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start
            ).count();

        if (!validateZeroVariableProblem()) {
            result_.status =
                QpStatus::Infeasible;

            result_.statusMessage =
                "Zero-variable problem is infeasible";

            result_.primal.clear();

            result_.constraintDual.assign(
                static_cast<std::size_t>(m),
                0.0
            );

            result_.primalObjective = 0.0;
            result_.bestObjective = 0.0;
            result_.primalResidual = 0.0;
            result_.dualResidual = 0.0;

            return result_;
        }

        result_.status =
            QpStatus::Optimal;

        result_.statusMessage =
            "SuperADMM solved zero-variable problem";

        result_.primal.clear();

        result_.constraintDual.assign(
            static_cast<std::size_t>(m),
            0.0
        );

        result_.primalObjective = 0.0;
        result_.bestObjective = 0.0;
        result_.primalResidual = 0.0;
        result_.dualResidual = 0.0;

        return result_;
    }

    if (options_.iterationLimit == 0) {
        result_.status =
            QpStatus::IterationLimit;

        result_.statusMessage =
            "SuperADMM iteration limit reached";

        result_.solveTimeSeconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start
            ).count();

        return result_;
    }

    for (std::int64_t iterationCount = 1;
         iterationCount <= options_.iterationLimit;
         ++iterationCount) {

        if (options_.timeLimitSeconds > 0.0) {
            const double elapsed =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start
                ).count();

            if (elapsed >=
                options_.timeLimitSeconds) {

                result_.status =
                    QpStatus::TimeLimit;

                result_.statusMessage =
                    "SuperADMM time limit reached";

                break;
            }
        }

        if (!iteration()) {
            result_.status =
                QpStatus::NumericalFailure;

            result_.statusMessage =
                "SuperADMM KKT solve failed";

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

        if (!std::isfinite(
                result_.primalObjective
            ) ||
            !std::isfinite(
                result_.primalResidual
            ) ||
            !std::isfinite(
                result_.dualResidual
            )) {

            result_.status =
                QpStatus::NumericalFailure;

            result_.statusMessage =
                "SuperADMM produced a "
                "non-finite iterate";

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
         * Standard SuperADMM stopping criterion:
         *
         *     r_prim <= tolerance
         *     r_dual <= tolerance
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
         * Infeasibility checks are performed every
         * 10 iterations by default, as specified in
         * the paper.
         */
        if (iterationCount %
                options_.infeasibilityCheckInterval ==
            0) {

            if (checkPrimalInfeasibility()) {
                result_.status =
                    QpStatus::Infeasible;

                result_.statusMessage =
                    "SuperADMM primal "
                    "infeasibility detected";

                break;
            }

            if (checkDualInfeasibility()) {
                result_.status =
                    QpStatus::Unbounded;

                result_.statusMessage =
                    "SuperADMM dual "
                    "infeasibility detected";

                break;
            }
        }

        /*
         * The paper terminates when b < 1 because the
         * stability interval
         *
         *     1/b <= rho_i <= b
         *
         * can no longer be satisfied.
         */
        if (b_ < 1.0) {
            result_.status =
                QpStatus::NumericalFailure;

            result_.statusMessage =
                "SuperADMM stability bound fell below 1";

            break;
        }
    }

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
     * This implementation does not calculate a separate
     * dual objective, so do not fabricate one.
     */
    result_.dualObjective =
        -std::numeric_limits<double>::infinity();

    result_.solveTimeSeconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start
        ).count();

    return result_;
}

bool SuperAdmmSolver::iteration() {
    const int n =
        model_.numVariables();

    const int m =
        model_.numConstraints();

    const std::vector<double> oldX =
        x_;

    const std::vector<double> oldZ =
        z_;

    const std::vector<double> oldY =
        y_;

    /*
     * Store x^{k-1}, y^{k-1} before computing the
     * current iteration. These are required for the
     * infeasibility certificates.
     */
    if (havePreviousIterate_) {
        previousX_ = oldX;
        previousY_ = oldY;
    }

    /*
     * ----------------------------------------------------
     * Equation (17a)
     *
     * [ P + sigma I      A^T ] [x^{k+1} ] =
     * [ sigma x^k - q       ]
     *
     * [ A             -R^-1 ] [nu^{k+1}] =
     * [ z^k - R^-1 y^k     ]
     * ----------------------------------------------------
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

        const double rho =
            rho_[index];

        if (!std::isfinite(rho) ||
            rho <= 0.0) {
            return false;
        }

        rhsNu[index] =
            oldZ[index] -
            oldY[index] / rho;
    }

    std::vector<double> newX;
    std::vector<double> newNu;

    if (!kktSolver_.solve(
            rho_,
            options_.sigma,
            rhsX,
            rhsNu,
            newX,
            newNu
        )) {
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
     * ----------------------------------------------------
     * Equation (17b)
     *
     * ztilde^{k+1}
     *
     *     = z^k + R_k^{-1}(nu^{k+1} - y^k)
     * ----------------------------------------------------
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        zTilde_[index] =
            oldZ[index] +
            (
                nu_[index] -
                oldY[index]
            ) /
            rho_[index];
    }

    /*
     * ----------------------------------------------------
     * Equation (17c)
     *
     * z^{k+1}
     *
     *   = Pi_[l,u](
     *       ztilde^{k+1}
     *       + R_k^{-1} y^k
     *     )
     * ----------------------------------------------------
     */
    project(
        zTilde_,
        oldY,
        z_
    );

    /*
     * ----------------------------------------------------
     * Equation (17d)
     *
     * y^{k+1}
     *
     *   = y^k + R_k(
     *       ztilde^{k+1} - z^{k+1}
     *     )
     * ----------------------------------------------------
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
     * Equation (19): numerical KKT error.
     */
    const double epsilon =
        kktResidual(
            oldX,
            oldZ,
            oldY
        );

    const double rPrim =
        primalResidual();

    if (!updateBound(
            epsilon,
            rPrim
        )) {
        return false;
    }

    /*
     * Algorithm 1 line 12 uses z^k, i.e. oldZ,
     * to determine which constraints are active.
     */
    if (!updateWeights(oldZ)) {
        return false;
    }

    havePreviousIterate_ = true;

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
            yOld[index] /
                rho_[index];

        const double lower =
            model_.l[index];

        const double upper =
            model_.u[index];

        double projected =
            value;

        if (std::isfinite(lower) &&
            projected < lower) {
            projected = lower;
        }

        if (std::isfinite(upper) &&
            projected > upper) {
            projected = upper;
        }

        zNew[index] =
            projected;
    }
}

bool SuperAdmmSolver::updateBound(
    double epsilon,
    double primalResidualValue
) {
    if (!std::isfinite(epsilon) ||
        !std::isfinite(primalResidualValue)) {
        return false;
    }

    /*
     * Equation (18):
     *
     * b^{k+1} =
     *
     *     tau b^k,  if epsilon >= r_prim
     *
     *     b^k,      otherwise.
     */
    if (epsilon >= primalResidualValue) {
        b_ *= options_.tau;
    }

    if (!std::isfinite(b_)) {
        return false;
    }

    return true;
}

bool SuperAdmmSolver::updateWeights(
    const std::vector<double>& oldZ
) {
    const int m =
        model_.numConstraints();

    /*
     * Once b < 1, the stability interval is empty.
     * Do not generate an invalid R.
     */
    if (b_ < 1.0) {
        return true;
    }

    const double lowerRho =
        1.0 / b_;

    const double upperRho =
        b_;

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double lower =
            model_.l[index];

        const double upper =
            model_.u[index];

        const double zValue =
            oldZ[index];

        const bool atLower =
            std::isfinite(lower) &&
            zValue == lower;

        const bool atUpper =
            std::isfinite(upper) &&
            zValue == upper;

        const bool active =
            atLower || atUpper;

        double newRho;

        if (active) {
            newRho =
                std::min(
                    upperRho,
                    options_.alpha *
                        rho_[index]
                );
        } else {
            newRho =
                std::max(
                    lowerRho,
                    rho_[index] /
                        options_.alpha
                );
        }

        /*
         * Explicitly enforce the numerical stability
         * interval from the paper.
         */
        newRho =
            std::clamp(
                newRho,
                lowerRho,
                upperRho
            );

        if (!std::isfinite(newRho) ||
            newRho <= 0.0) {
            return false;
        }

        rho_[index] =
            newRho;
    }

    return true;
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
    std::vector<double> ATnu;
    std::vector<double> Ax;

    model_.P.multiply(
        x_,
        Px
    );

    model_.A.transposeMultiply(
        nu_,
        ATnu
    );

    model_.A.multiply(
        x_,
        Ax
    );

    double epsilon =
        0.0;

    /*
     * Top KKT equation:
     *
     * sigma*x^k - q
     *
     * -
     *
     * ((P + sigma I)x^{k+1}
     *      + A^T nu^{k+1})
     */
    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        const double residual =
            options_.sigma *
                oldX[index] -
            model_.q[index] -
            (
                Px[index] +
                options_.sigma *
                    x_[index] +
                ATnu[index]
            );

        epsilon =
            std::max(
                epsilon,
                std::abs(residual)
            );
    }

    /*
     * Bottom KKT equation:
     *
     * z^k - R^{-1}y^k
     *
     * -
     *
     * (Ax^{k+1} - R^{-1}nu^{k+1})
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double rho =
            rho_[index];

        const double residual =
            oldZ[index] -
            oldY[index] / rho -
            Ax[index] +
            nu_[index] / rho;

        epsilon =
            std::max(
                epsilon,
                std::abs(residual)
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

    if (previousY_.size() != y_.size()) {
        return false;
    }

    std::vector<double> deltaY(
        static_cast<std::size_t>(m),
        0.0
    );

    bool nonzero =
        false;

    double deltaYScale =
        0.0;

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        deltaY[index] =
            y_[index] -
            previousY_[index];

        deltaYScale =
            std::max(
                deltaYScale,
                std::abs(deltaY[index])
            );

        if (std::abs(deltaY[index]) >
            options_.infeasibilityTolerance) {
            nonzero = true;
        }
    }

    if (!nonzero) {
        return false;
    }

    std::vector<double> ATdeltaY;

    model_.A.transposeMultiply(
        deltaY,
        ATdeltaY
    );

    double aNorm =
        0.0;

    const auto& aValues =
        model_.A.csrValues();

    for (double value : aValues) {
        aNorm =
            std::max(
                aNorm,
                std::abs(value)
            );
    }

    const double directionTolerance =
        options_.infeasibilityTolerance *
        std::max(
            1.0,
            aNorm * deltaYScale
        );

    for (double value : ATdeltaY) {
        if (std::abs(value) >
            directionTolerance) {
            return false;
        }
    }

    /*
     * Certificate:
     *
     *     u^T deltaY_+
     *   + l^T deltaY_-
     *   < 0
     *
     * where deltaY_+ contains positive entries and
     * deltaY_- contains negative entries.
     */
    double certificate =
        0.0;

    double certificateScale =
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

            certificateScale +=
                std::abs(upper * dy);

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

            certificateScale +=
                std::abs(lower * dy);
        }
    }

    const double certificateTolerance =
        options_.infeasibilityTolerance *
        std::max(
            1.0,
            certificateScale
        );

    return certificate <
           -certificateTolerance;
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

    if (previousX_.size() != x_.size()) {
        return false;
    }

    std::vector<double> deltaX(
        static_cast<std::size_t>(n),
        0.0
    );

    bool nonzero =
        false;

    double deltaXScale =
        0.0;

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        deltaX[index] =
            x_[index] -
            previousX_[index];

        deltaXScale =
            std::max(
                deltaXScale,
                std::abs(deltaX[index])
            );

        if (std::abs(deltaX[index]) >
            options_.infeasibilityTolerance) {
            nonzero = true;
        }
    }

    if (!nonzero) {
        return false;
    }

    /*
     * q^T deltaX < 0.
     */
    double qDeltaX =
        0.0;

    double qScale =
        0.0;

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        qDeltaX +=
            model_.q[index] *
            deltaX[index];

        qScale +=
            std::abs(
                model_.q[index] *
                deltaX[index]
            );
    }

    const double qTolerance =
        options_.infeasibilityTolerance *
        std::max(
            1.0,
            qScale
        );

    if (qDeltaX >=
        -qTolerance) {
        return false;
    }

    /*
     * P deltaX = 0.
     */
    std::vector<double> PdeltaX;

    model_.P.multiply(
        deltaX,
        PdeltaX
    );

    double pNorm =
        0.0;

    for (double value :
         model_.P.csrValues()) {
        pNorm =
            std::max(
                pNorm,
                std::abs(value)
            );
    }

    const double pTolerance =
        options_.infeasibilityTolerance *
        std::max(
            1.0,
            pNorm * deltaXScale
        );

    for (double value :
         PdeltaX) {
        if (std::abs(value) >
            pTolerance) {
            return false;
        }
    }

    /*
     * Constraint-direction conditions:
     *
     * finite lower + finite upper:
     *     A_i deltaX = 0
     *
     * finite lower + infinite upper:
     *     A_i deltaX >= 0
     *
     * infinite lower + finite upper:
     *     A_i deltaX <= 0
     *
     * both infinite:
     *     no restriction.
     */
    if (m > 0) {
        std::vector<double> AdeltaX;

        model_.A.multiply(
            deltaX,
            AdeltaX
        );

        double aNorm =
            0.0;

        for (double value :
             model_.A.csrValues()) {
            aNorm =
                std::max(
                    aNorm,
                    std::abs(value)
                );
        }

        const double aTolerance =
            options_.infeasibilityTolerance *
            std::max(
                1.0,
                aNorm * deltaXScale
            );

        for (int i = 0; i < m; ++i) {
            const std::size_t index =
                static_cast<std::size_t>(i);

            const double value =
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

                if (std::abs(value) >
                    aTolerance) {
                    return false;
                }

            } else if (
                lowerFinite &&
                !upperFinite
            ) {

                if (value <
                    -aTolerance) {
                    return false;
                }

            } else if (
                !lowerFinite &&
                upperFinite
            ) {

                if (value >
                    aTolerance) {
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