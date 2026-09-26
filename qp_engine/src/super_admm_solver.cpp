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

    if (options_.timeLimitSeconds < 0.0) {
        throw std::invalid_argument(
            "SuperADMM timeLimitSeconds must be non-negative"
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
    result_.statusMessage = "SuperADMM iteration limit reached";
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

    if (m > 0) {
        model_.A.multiply(x_, z_);
    }

    result_.primalResidual = primalResidual();
    result_.dualResidual = dualResidual();

    result_.primalObjective = objective(x_);
    result_.bestObjective = result_.primalObjective;

    result_.finalRho =
        rho_.empty()
            ? 0.0
            : rho_.front();

    return true;
}

AdmmResult SuperAdmmSolver::solve() {
    const auto start =
        std::chrono::steady_clock::now();

    try {
        model_.validate();
    } catch (const std::exception& e) {
        result_.status = QpStatus::InvalidProblem;
        result_.statusMessage = e.what();
        return result_;
    }

    if (!initialize()) {
        result_.status = QpStatus::NumericalFailure;
        result_.statusMessage =
            "SuperADMM initialization failed";
        return result_;
    }

    const int n = model_.numVariables();
    const int m = model_.numConstraints();

    if (n == 0) {
        result_.status = QpStatus::Optimal;
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
        result_.solveTimeSeconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start
            ).count();

        return result_;
    }

    if (options_.iterationLimit == 0) {
        result_.status = QpStatus::IterationLimit;
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

            if (elapsed >= options_.timeLimitSeconds) {
                result_.status = QpStatus::TimeLimit;
                result_.statusMessage =
                    "SuperADMM time limit reached";
                break;
            }
        }

        if (!iteration()) {
            result_.status = QpStatus::NumericalFailure;
            result_.statusMessage =
                "SuperADMM KKT solve failed";
            break;
        }

        result_.iterations = iterationCount;

        result_.primal = x_;
        result_.constraintDual = y_;

        result_.primalResidual = primalResidual();
        result_.dualResidual = dualResidual();

        result_.primalObjective = objective(x_);
        result_.bestObjective = result_.primalObjective;

        if (!std::isfinite(result_.primalObjective) ||
            !std::isfinite(result_.primalResidual) ||
            !std::isfinite(result_.dualResidual)) {

            result_.status = QpStatus::NumericalFailure;
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
         * Check optimality.
         *
         * SuperADMM uses the same primal and dual residual
         * definitions as the underlying QP:
         *
         *   r_prim = ||Ax - z||_inf
         *   r_dual = ||Px + q + A^T y||_inf
         */
        if (result_.primalResidual <=
                options_.primalTolerance &&
            result_.dualResidual <=
                options_.dualTolerance) {

            result_.status = QpStatus::Optimal;
            result_.statusMessage =
                "SuperADMM converged";

            break;
        }

        /*
         * Infeasibility checks are performed at the
         * requested interval, following the paper.
         */
        if (iterationCount %
                options_.infeasibilityCheckInterval ==
            0) {

            if (checkPrimalInfeasibility()) {
                result_.status = QpStatus::Infeasible;
                result_.statusMessage =
                    "SuperADMM primal infeasibility detected";
                break;
            }

            if (checkDualInfeasibility()) {
                result_.status = QpStatus::Unbounded;
                result_.statusMessage =
                    "SuperADMM dual infeasibility detected";
                break;
            }
        }

        /*
         * Once b drops below 1, the numerical-stability
         * safeguard has reached its stopping condition.
         */
        if (b_ < 1.0) {
            result_.status = QpStatus::NumericalFailure;
            result_.statusMessage =
                "SuperADMM stability bound fell below 1";
            break;
        }
    }

    result_.primal = x_;
    result_.constraintDual = y_;

    result_.primalResidual = primalResidual();
    result_.dualResidual = dualResidual();

    result_.primalObjective = objective(x_);
    result_.bestObjective = result_.primalObjective;

    result_.finalRho =
        rho_.empty()
            ? 0.0
            : *std::max_element(
                  rho_.begin(),
                  rho_.end()
              );

    /*
     * SuperADMM does not construct a separate dual
     * objective in this implementation. Therefore we
     * explicitly leave the standard AdmmResult dual
     * objective at -infinity rather than fabricating a
     * value.
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
    const int n = model_.numVariables();
    const int m = model_.numConstraints();

    const std::vector<double> oldX = x_;
    const std::vector<double> oldZ = z_;
    const std::vector<double> oldY = y_;

    /*
     * Keep the previous iterate for infeasibility
     * certificates.
     */
    if (havePreviousIterate_) {
        previousX_ = oldX;
        previousY_ = oldY;
    }

    /*
     * ----------------------------------------------------
     * KKT system
     * ----------------------------------------------------
     *
     * [ P + sigma I     A^T      ] [x     ] =
     * [ sigma*x^k - q ]
     * [ A              -R^{-1}  ] [nu    ]   [z^k - R^{-1}y^k]
     */
    std::vector<double> rhsX(
        static_cast<std::size_t>(n),
        0.0
    );

    for (int j = 0; j < n; ++j) {
        rhsX[static_cast<std::size_t>(j)] =
            options_.sigma *
                oldX[static_cast<std::size_t>(j)] -
            model_.q[static_cast<std::size_t>(j)];
    }

    std::vector<double> rhsNu(
        static_cast<std::size_t>(m),
        0.0
    );

    for (int i = 0; i < m; ++i) {
        const double ri =
            rho_[static_cast<std::size_t>(i)];

        if (ri <= 0.0 ||
            !std::isfinite(ri)) {
            return false;
        }

        rhsNu[static_cast<std::size_t>(i)] =
            oldZ[static_cast<std::size_t>(i)] -
            oldY[static_cast<std::size_t>(i)] / ri;
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

    x_ = std::move(newX);
    nu_ = std::move(newNu);

    /*
     * ----------------------------------------------------
     * z-tilde update
     * ----------------------------------------------------
     *
     * ztilde^{k+1}
     *   = z^k + R^{-1}(nu^{k+1} - y^k)
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        zTilde_[index] =
            oldZ[index] +
            (nu_[index] - oldY[index]) /
                rho_[index];
    }

    /*
     * ----------------------------------------------------
     * Projection
     * ----------------------------------------------------
     *
     * z^{k+1}
     *   = Pi_[l,u](ztilde^{k+1} + R^{-1}y^k)
     */
    project(
        zTilde_,
        oldY,
        z_
    );

    /*
     * ----------------------------------------------------
     * Dual update
     * ----------------------------------------------------
     *
     * y^{k+1}
     *   = y^k + R(ztilde^{k+1} - z^{k+1})
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        y_[index] =
            oldY[index] +
            rho_[index] *
                (zTilde_[index] - z_[index]);
    }

    if (!finiteVector(x_) ||
        !finiteVector(z_) ||
        !finiteVector(y_) ||
        !finiteVector(nu_)) {
        return false;
    }

    /*
     * Compute the KKT residual before changing R.
     *
     * This epsilon corresponds to the residual of the
     * direct KKT system used by SuperADMM.
     */
    const double epsilon =
        kktResidual(
            oldX,
            oldZ,
            oldY
        );

    const double rPrim =
        primalResidual();

    /*
     * Algorithm 1:
     *
     * b^{k+1} =
     *     tau*b^k   if epsilon >= r_prim
     *     b^k       otherwise
     */
    updateBound(
        epsilon,
        rPrim
    );

    /*
     * IMPORTANT:
     *
     * The paper updates R using the OLD z^k,
     * not the newly projected z^{k+1}.
     */
    updateWeights(oldZ);

    havePreviousIterate_ = true;

    return true;
}

void SuperAdmmSolver::project(
    const std::vector<double>& zTilde,
    const std::vector<double>& yOld,
    std::vector<double>& zNew
) const {
    const int m = model_.numConstraints();

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

        double projected = value;

        if (std::isfinite(lower) &&
            projected < lower) {
            projected = lower;
        }

        if (std::isfinite(upper) &&
            projected > upper) {
            projected = upper;
        }

        zNew[index] = projected;
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

    /*
     * The paper's stability condition requires
     *
     *     1 / b <= R_i <= b.
     *
     * Once b becomes smaller than 1, this condition
     * cannot be satisfied. We therefore allow the main
     * loop to terminate on the next check.
     */
}

void SuperAdmmSolver::updateWeights(
    const std::vector<double>& oldZ
) {
    const int m = model_.numConstraints();

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double lower =
            model_.l[index];

        const double upper =
            model_.u[index];

        const double zi =
            oldZ[index];

        const bool atLower =
            std::isfinite(lower) &&
            zi == lower;

        const bool atUpper =
            std::isfinite(upper) &&
            zi == upper;

        const bool active =
            atLower || atUpper;

        double newRho;

        if (active) {
            newRho =
                std::min(
                    b_,
                    options_.alpha * rho_[index]
                );
        } else {
            newRho =
                std::max(
                    1.0 / b_,
                    rho_[index] / options_.alpha
                );
        }

        if (!std::isfinite(newRho) ||
            newRho <= 0.0) {
            /*
             * Do not allow an invalid R value to enter
             * the next KKT system.
             */
            newRho =
                std::max(
                    1.0 / std::max(b_, 1.0),
                    std::numeric_limits<double>::min()
                );
        }

        rho_[index] = newRho;
    }
}

double SuperAdmmSolver::primalResidual() const {
    const int m = model_.numConstraints();

    if (m == 0) {
        return 0.0;
    }

    std::vector<double> Ax;

    model_.A.multiply(
        x_,
        Ax
    );

    double residual = 0.0;

    for (int i = 0; i < m; ++i) {
        residual =
            std::max(
                residual,
                std::abs(
                    Ax[static_cast<std::size_t>(i)] -
                    z_[static_cast<std::size_t>(i)]
                )
            );
    }

    return residual;
}

double SuperAdmmSolver::dualResidual() const {
    const int n = model_.numVariables();

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

    double residual = 0.0;

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
    const int n = model_.numVariables();
    const int m = model_.numConstraints();

    std::vector<double> Px;
    std::vector<double> Aty;
    std::vector<double> Ax;

    model_.P.multiply(
        x_,
        Px
    );

    model_.A.transposeMultiply(
        nu_,
        Aty
    );

    model_.A.multiply(
        x_,
        Ax
    );

    double epsilon = 0.0;

    /*
     * Top KKT block:
     *
     * (sigma*x^k - q)
     * -
     * ((P + sigma I)x^{k+1} + A^T nu^{k+1})
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
                options_.sigma * x_[index] +
                Aty[index]
            );

        epsilon =
            std::max(
                epsilon,
                std::abs(residual)
            );
    }

    /*
     * Bottom KKT block:
     *
     * z^k - R^{-1}y^k
     * -
     * Ax^{k+1}
     * +
     * R^{-1}nu^{k+1}
     */
    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double ri =
            rho_[index];

        const double residual =
            oldZ[index] -
            oldY[index] / ri -
            Ax[index] +
            nu_[index] / ri;

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
    const int n = model_.numVariables();

    std::vector<double> Px;

    model_.P.multiply(
        x,
        Px
    );

    double value = 0.0;

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
    const int m = model_.numConstraints();

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

    bool hasNonzeroDirection = false;

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        deltaY[index] =
            y_[index] -
            previousY_[index];

        if (std::abs(deltaY[index]) >
            options_.infeasibilityTolerance) {
            hasNonzeroDirection = true;
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

    double certificate = 0.0;

    for (int i = 0; i < m; ++i) {
        const std::size_t index =
            static_cast<std::size_t>(i);

        const double dy =
            deltaY[index];

        if (dy > options_.infeasibilityTolerance) {
            const double upper =
                model_.u[index];

            if (!std::isfinite(upper)) {
                return false;
            }

            certificate +=
                upper * dy;
        } else if (
            dy < -options_.infeasibilityTolerance
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
    const int n = model_.numVariables();
    const int m = model_.numConstraints();

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

    bool hasNonzeroDirection = false;

    for (int j = 0; j < n; ++j) {
        const std::size_t index =
            static_cast<std::size_t>(j);

        deltaX[index] =
            x_[index] -
            previousX_[index];

        if (std::abs(deltaX[index]) >
            options_.infeasibilityTolerance) {
            hasNonzeroDirection = true;
        }
    }

    if (!hasNonzeroDirection) {
        return false;
    }

    /*
     * q^T delta_x < 0
     */
    double qDeltaX = 0.0;

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
     *     A_i delta_x = 0
     *
     * lower finite, upper infinite:
     *     A_i delta_x >= 0
     *
     * lower infinite, upper finite:
     *     A_i delta_x <= 0
     *
     * both infinite:
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
                std::isfinite(model_.l[index]);

            const bool upperFinite =
                std::isfinite(model_.u[index]);

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