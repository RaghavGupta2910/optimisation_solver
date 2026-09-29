#include "qp/super_admm_kkt_solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace qp {
namespace {

constexpr double kPivotTolerance = 1e-12;

/**
 * Dense LU factorization with partial pivoting.
 *
 * The matrix is overwritten by its LU factors.
 * pivots[k] stores the row selected as the pivot at step k.
 */
bool luFactor(
    std::vector<double>& matrix,
    int n,
    std::vector<int>& pivots
) {
    if (n == 0) {
        pivots.clear();
        return true;
    }

    pivots.resize(static_cast<std::size_t>(n));

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;
        double pivotAbs = std::abs(
            matrix[static_cast<std::size_t>(k) * n + k]
        );

        for (int i = k + 1; i < n; ++i) {
            const double value = std::abs(
                matrix[static_cast<std::size_t>(i) * n + k]
            );

            if (value > pivotAbs) {
                pivotAbs = value;
                pivotRow = i;
            }
        }

        if (!std::isfinite(pivotAbs) ||
            pivotAbs <= kPivotTolerance) {
            return false;
        }

        pivots[static_cast<std::size_t>(k)] = pivotRow;

        if (pivotRow != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(
                    matrix[static_cast<std::size_t>(k) * n + j],
                    matrix[static_cast<std::size_t>(pivotRow) * n + j]
                );
            }
        }

        const double pivot =
            matrix[static_cast<std::size_t>(k) * n + k];

        for (int i = k + 1; i < n; ++i) {
            double& multiplier =
                matrix[static_cast<std::size_t>(i) * n + k];

            multiplier /= pivot;

            if (!std::isfinite(multiplier)) {
                return false;
            }

            for (int j = k + 1; j < n; ++j) {
                matrix[static_cast<std::size_t>(i) * n + j] -=
                    multiplier *
                    matrix[static_cast<std::size_t>(k) * n + j];
            }
        }
    }

    return true;
}

/**
 * Solve LU x = rhs after luFactor().
 */
bool luSolve(
    const std::vector<double>& lu,
    int n,
    const std::vector<int>& pivots,
    std::vector<double>& rhs
) {
    if (n == 0) {
        return true;
    }

    if (static_cast<int>(rhs.size()) != n ||
        static_cast<int>(pivots.size()) != n ||
        static_cast<int>(lu.size()) != n * n) {
        return false;
    }

    // Apply row permutations and solve Ly = Pb.
    for (int k = 0; k < n; ++k) {
        const int pivotRow = pivots[static_cast<std::size_t>(k)];

        if (pivotRow != k) {
            std::swap(
                rhs[static_cast<std::size_t>(k)],
                rhs[static_cast<std::size_t>(pivotRow)]
            );
        }

        for (int i = k + 1; i < n; ++i) {
            rhs[static_cast<std::size_t>(i)] -=
                lu[static_cast<std::size_t>(i) * n + k] *
                rhs[static_cast<std::size_t>(k)];
        }
    }

    // Solve Ux = y.
    for (int i = n - 1; i >= 0; --i) {
        double value = rhs[static_cast<std::size_t>(i)];

        for (int j = i + 1; j < n; ++j) {
            value -=
                lu[static_cast<std::size_t>(i) * n + j] *
                rhs[static_cast<std::size_t>(j)];
        }

        const double diagonal =
            lu[static_cast<std::size_t>(i) * n + i];

        if (!std::isfinite(diagonal) ||
            std::abs(diagonal) <= kPivotTolerance) {
            return false;
        }

        rhs[static_cast<std::size_t>(i)] =
            value / diagonal;

        if (!std::isfinite(rhs[static_cast<std::size_t>(i)])) {
            return false;
        }
    }

    return true;
}

}  // namespace

SuperAdmmKktSolver::SuperAdmmKktSolver(const QpModel& model)
    : model_(&model),
      n_(model.numVariables()),
      m_(model.numConstraints()) {
}

bool SuperAdmmKktSolver::solve(
    const std::vector<double>& rho,
    double sigma,
    const std::vector<double>& rhsX,
    const std::vector<double>& rhsNu,
    std::vector<double>& x,
    std::vector<double>& nu
) const {
    if (model_ == nullptr) {
        return false;
    }

    if (!std::isfinite(sigma) || sigma <= 0.0) {
        return false;
    }

    if (static_cast<int>(rho.size()) != m_) {
        return false;
    }

    if (static_cast<int>(rhsX.size()) != n_) {
        return false;
    }

    if (static_cast<int>(rhsNu.size()) != m_) {
        return false;
    }

    for (double value : rho) {
        if (!std::isfinite(value) || value <= 0.0) {
            return false;
        }
    }

    for (double value : rhsX) {
        if (!std::isfinite(value)) {
            return false;
        }
    }

    for (double value : rhsNu) {
        if (!std::isfinite(value)) {
            return false;
        }
    }

    const int dimension = n_ + m_;

    /*
     * There are no variables and no constraints.
     *
     * The KKT system is empty, so there is nothing to solve.
     */
    if (dimension == 0) {
        x.clear();
        nu.clear();
        return true;
    }

    /*
     * Dense KKT matrix:
     *
     *       [ P + sigma I       A^T       ]
     * K  =  [                         ]
     *       [ A                -R^{-1}    ]
     *
     * where R = diag(rho).
     */
    std::vector<double> kkt(
        static_cast<std::size_t>(dimension) * dimension,
        0.0
    );

    // ---------------------------------------------------------------------
    // Top-left block: P + sigma I
    // ---------------------------------------------------------------------
    const auto& pRowStart = model_->P.csrRowStart();
    const auto& pColumnIndex = model_->P.csrColumnIndex();
    const auto& pValues = model_->P.csrValues();

    for (int row = 0; row < n_; ++row) {
        const int begin = pRowStart[static_cast<std::size_t>(row)];
        const int end = pRowStart[static_cast<std::size_t>(row + 1)];

        for (int k = begin; k < end; ++k) {
            const int col =
                pColumnIndex[static_cast<std::size_t>(k)];

            if (col < 0 || col >= n_) {
                return false;
            }

            const double value =
                pValues[static_cast<std::size_t>(k)];

            if (!std::isfinite(value)) {
                return false;
            }

            kkt[
                static_cast<std::size_t>(row) * dimension + col
            ] += value;
        }

        kkt[
            static_cast<std::size_t>(row) * dimension + row
        ] += sigma;
    }

    // ---------------------------------------------------------------------
    // Top-right and bottom-left blocks: A^T and A
    // ---------------------------------------------------------------------
    const auto& aRowStart = model_->A.csrRowStart();
    const auto& aColumnIndex = model_->A.csrColumnIndex();
    const auto& aValues = model_->A.csrValues();

    for (int row = 0; row < m_; ++row) {
        const int begin = aRowStart[static_cast<std::size_t>(row)];
        const int end = aRowStart[static_cast<std::size_t>(row + 1)];

        for (int k = begin; k < end; ++k) {
            const int col =
                aColumnIndex[static_cast<std::size_t>(k)];

            if (col < 0 || col >= n_) {
                return false;
            }

            const double value =
                aValues[static_cast<std::size_t>(k)];

            if (!std::isfinite(value)) {
                return false;
            }

            // Bottom-left: A.
            kkt[
                static_cast<std::size_t>(n_ + row) * dimension + col
            ] += value;

            // Top-right: A^T.
            kkt[
                static_cast<std::size_t>(col) * dimension +
                (n_ + row)
            ] += value;
        }
    }

    // ---------------------------------------------------------------------
    // Bottom-right block: -R^{-1}
    // ---------------------------------------------------------------------
    for (int i = 0; i < m_; ++i) {
        const double inverseRho =
            1.0 / rho[static_cast<std::size_t>(i)];

        if (!std::isfinite(inverseRho)) {
            return false;
        }

        kkt[
            static_cast<std::size_t>(n_ + i) * dimension +
            (n_ + i)
        ] = -inverseRho;
    }

    // ---------------------------------------------------------------------
    // Right-hand side.
    // ---------------------------------------------------------------------
    std::vector<double> rhs(
        static_cast<std::size_t>(dimension),
        0.0
    );

    for (int i = 0; i < n_; ++i) {
        rhs[static_cast<std::size_t>(i)] =
            rhsX[static_cast<std::size_t>(i)];
    }

    for (int i = 0; i < m_; ++i) {
        rhs[static_cast<std::size_t>(n_ + i)] =
            rhsNu[static_cast<std::size_t>(i)];
    }

    // ---------------------------------------------------------------------
    // LU factorization.
    // ---------------------------------------------------------------------
    std::vector<int> pivots;

    if (!luFactor(kkt, dimension, pivots)) {
        return false;
    }

    // ---------------------------------------------------------------------
    // Solve.
    // ---------------------------------------------------------------------
    if (!luSolve(kkt, dimension, pivots, rhs)) {
        return false;
    }

    // ---------------------------------------------------------------------
    // Extract x and nu.
    // ---------------------------------------------------------------------
    x.resize(static_cast<std::size_t>(n_));
    nu.resize(static_cast<std::size_t>(m_));

    for (int i = 0; i < n_; ++i) {
        x[static_cast<std::size_t>(i)] =
            rhs[static_cast<std::size_t>(i)];

        if (!std::isfinite(x[static_cast<std::size_t>(i)])) {
            return false;
        }
    }

    for (int i = 0; i < m_; ++i) {
        nu[static_cast<std::size_t>(i)] =
            rhs[static_cast<std::size_t>(n_ + i)];

        if (!std::isfinite(nu[static_cast<std::size_t>(i)])) {
            return false;
        }
    }

    return true;
}

}  // namespace qp