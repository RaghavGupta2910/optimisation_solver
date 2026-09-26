#pragma once

#include "qp/qp_model.h"

#include <vector>

namespace qp {

// -----------------------------------------------------------------------------
// Direct SuperADMM KKT solver
//
// Solves
//
//   [ P + sigma I     A^T ] [ x  ]   [ rhsX  ]
//   [ A              -R^-1] [ nu ] = [ rhsNu ]
//
// where R is diagonal and supplied as a vector containing its diagonal.
//
// This first implementation deliberately uses a dense LU factorisation.
// It is correctness-first and is suitable for validating the algorithm.
// A sparse LDLT implementation can be added later without changing the
// SuperAdmmSolver interface.
// -----------------------------------------------------------------------------
class SuperAdmmKktSolver {
public:
    explicit SuperAdmmKktSolver(const QpModel& model);

    bool solve(
        const std::vector<double>& rho,
        double sigma,
        const std::vector<double>& rhsX,
        const std::vector<double>& rhsNu,
        std::vector<double>& x,
        std::vector<double>& nu
    ) const;

private:
    const QpModel* model_ = nullptr;

    int n_ = 0;
    int m_ = 0;
};

}  // namespace qp