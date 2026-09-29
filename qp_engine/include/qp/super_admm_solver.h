#pragma once

#include "qp/admm_solver.h"
#include "qp/super_admm_kkt_solver.h"

#include <cstdint>
#include <vector>

namespace qp {

/**
 * @brief Options controlling the SuperADMM solver.
 *
 * SuperADMM solves convex quadratic programs of the form
 *
 *     minimize    1/2 x^T P x + q^T x
 *     subject to  l <= A x <= u
 *
 * using a dynamically weighted ADMM iteration.
 */
struct SuperAdmmOptions {
    // ---------------------------------------------------------------------
    // Termination limits
    // ---------------------------------------------------------------------

    std::int64_t iterationLimit = 5000;

    // 0.0 means no time limit.
    double timeLimitSeconds = 0.0;

    double primalTolerance = 1e-6;
    double dualTolerance = 1e-6;

    // ---------------------------------------------------------------------
    // SuperADMM parameters
    // ---------------------------------------------------------------------

    // alpha controls how quickly individual constraint weights can change.
    double alpha = 500.0;

    // Positive regularization parameter used in the KKT system.
    double sigma = 1e-6;

    // Initial numerical-stability bound.
    double b0 = 1e8;

    // Stability-bound reduction factor.
    double tau = 0.5;

    // Initial value of every diagonal entry of R.
    double rho0 = 1.0;

    // ---------------------------------------------------------------------
    // Infeasibility detection
    // ---------------------------------------------------------------------

    // Check infeasibility certificates every N iterations.
    int infeasibilityCheckInterval = 10;

    // Numerical tolerance used by the infeasibility tests.
    double infeasibilityTolerance = 1e-8;
};

/**
 * @brief SuperADMM solver for convex quadratic programs.
 *
 * The implementation follows the SuperADMM formulation in which the
 * scalar ADMM penalty is replaced by a diagonal matrix
 *
 *     R = diag(rho_1, ..., rho_m)
 *
 * whose entries are dynamically increased for active constraints and
 * decreased for inactive constraints.
 *
 * This solver is implemented separately from the repository's existing
 * AdmmSolver. The existing ADMM implementation is not modified.
 */
class SuperAdmmSolver {
public:
    explicit SuperAdmmSolver(
        const QpModel& model,
        const SuperAdmmOptions& options = {}
    );

    /**
     * @brief Solve the quadratic program.
     */
    AdmmResult solve();

private:
    // ---------------------------------------------------------------------
    // Initialization and iteration
    // ---------------------------------------------------------------------

    bool initialize();

    bool iteration();

    // ---------------------------------------------------------------------
    // SuperADMM operations
    // ---------------------------------------------------------------------

    void project(
        const std::vector<double>& zTilde,
        const std::vector<double>& yOld,
        std::vector<double>& zNew
    ) const;

    void updateBound(
        double epsilon,
        double primalResidual
    );

    void updateWeights(
        const std::vector<double>& oldZ
    );

    // ---------------------------------------------------------------------
    // Residuals and objective
    // ---------------------------------------------------------------------

    double primalResidual() const;

    double dualResidual() const;

    double kktResidual(
        const std::vector<double>& oldX,
        const std::vector<double>& oldZ,
        const std::vector<double>& oldY
    ) const;

    double objective(
        const std::vector<double>& x
    ) const;

    // ---------------------------------------------------------------------
    // Infeasibility detection
    // ---------------------------------------------------------------------

    bool checkPrimalInfeasibility();

    bool checkDualInfeasibility();

    // ---------------------------------------------------------------------
    // Utility
    // ---------------------------------------------------------------------

    bool finiteVector(
        const std::vector<double>& values
    ) const;

private:
    const QpModel& model_;

    SuperAdmmOptions options_;

    SuperAdmmKktSolver kktSolver_;

    // ---------------------------------------------------------------------
    // SuperADMM iterates
    // ---------------------------------------------------------------------

    // Primal variable.
    std::vector<double> x_;

    // Projected constraint variable.
    std::vector<double> z_;

    // Scaled dual variable.
    std::vector<double> y_;

    // KKT multiplier.
    std::vector<double> nu_;

    // Intermediate constraint value before projection.
    std::vector<double> zTilde_;

    // Diagonal entries of R.
    std::vector<double> rho_;

    // Numerical-stability bound b.
    double b_ = 0.0;

    // ---------------------------------------------------------------------
    // Previous iterates used by infeasibility certificates.
    // ---------------------------------------------------------------------

    std::vector<double> previousX_;

    std::vector<double> previousY_;

    bool havePreviousIterate_ = false;

    // ---------------------------------------------------------------------
    // Solver result
    // ---------------------------------------------------------------------

    AdmmResult result_;
};

}  // namespace qp