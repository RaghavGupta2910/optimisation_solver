#pragma once

#include "qp/admm_solver.h"
#include "qp/super_admm_kkt_solver.h"

#include <cstdint>
#include <vector>

namespace qp {

// -----------------------------------------------------------------------------
// SuperADMM options.
//
// Defaults follow the SuperADMM paper:
//   alpha = 500
//   sigma = 1e-6
//   b0    = 1e8
//   tau   = 0.5
// -----------------------------------------------------------------------------
struct SuperAdmmOptions {
    std::int64_t iterationLimit = 5000;
    double timeLimitSeconds = 0.0;

    double primalTolerance = 1e-6;
    double dualTolerance = 1e-6;

    double alpha = 500.0;
    double sigma = 1e-6;
    double b0 = 1e8;
    double tau = 0.5;
    double rho0 = 1.0;

    int infeasibilityCheckInterval = 10;
    double infeasibilityTolerance = 1e-8;
};

// -----------------------------------------------------------------------------
// SuperADMM solver.
//
// Solves
//
//   minimize    0.5 x^T P x + q^T x
//   subject to  l <= A x <= u
//
// using the direct SuperADMM formulation.
// -----------------------------------------------------------------------------
class SuperAdmmSolver {
public:
    explicit SuperAdmmSolver(
        const QpModel& model,
        const SuperAdmmOptions& options = {}
    );

    AdmmResult solve();

private:
    bool initialize();

    bool iteration();

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

    bool checkPrimalInfeasibility();

    bool checkDualInfeasibility();

    bool finiteVector(
        const std::vector<double>& values
    ) const;

private:
    const QpModel& model_;
    SuperAdmmOptions options_;

    SuperAdmmKktSolver kktSolver_;

    std::vector<double> x_;
    std::vector<double> z_;
    std::vector<double> y_;
    std::vector<double> nu_;
    std::vector<double> zTilde_;

    // Diagonal entries of R.
    std::vector<double> rho_;

    // Numerical stability bound b^k.
    double b_ = 0.0;

    // Previous iterate for infeasibility certificates.
    std::vector<double> previousX_;
    std::vector<double> previousY_;

    bool havePreviousIterate_ = false;

    AdmmResult result_;
};

}  // namespace qp