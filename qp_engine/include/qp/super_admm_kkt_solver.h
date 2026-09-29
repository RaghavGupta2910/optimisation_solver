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
 * where R is a positive diagonal matrix containing the
 * per-constraint SuperADMM penalty weights.
 *
 * The current implementation uses a dense LU factorization with
 * partial pivoting. This is intended as a correctness-first
 * reference implementation. It is not intended to be the
 * high-performance sparse backend for large-scale QPs.
 */
class SuperAdmmKktSolver {
public:
    explicit SuperAdmmKktSolver(const QpModel& model);

    /**
     * @brief Solve one SuperADMM KKT system.
     *
     * @param rho    Diagonal entries of R.
     * @param sigma  Positive regularization parameter.
     * @param rhsX   Right-hand side corresponding to x.
     * @param rhsNu  Right-hand side corresponding to nu.
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