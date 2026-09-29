#pragma once

#include "qp/qp_model.h"

#include <vector>

namespace qp {

/**
 * @brief Direct KKT solver used by SuperADMM.
 *
 * Solves the SuperADMM linear system
 *
 *     [ P + sigma I      A^T ] [x ] = [rhsX ]
 *     [ A             -R^-1 ] [nu]   [rhsNu]
 *
 * where R is a positive diagonal matrix containing
 * the per-constraint SuperADMM penalty weights.
 *
 * This implementation uses dense LU factorization with
 * partial pivoting. It is intended as a correctness-first
 * reference implementation.
 */
class SuperAdmmKktSolver {
public:
    explicit SuperAdmmKktSolver(const QpModel& model);

    /**
     * @brief Solve one SuperADMM KKT system.
     *
     * @param rho    Diagonal entries of R.
     * @param sigma  Positive proximal regularization.
     * @param rhsX   Right-hand side for the x block.
     * @param rhsNu  Right-hand side for the nu block.
     * @param x      Solution vector for x.
     * @param nu     Solution vector for nu.
     *
     * @return true if the system was solved successfully.
     */
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