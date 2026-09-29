#include "qp/super_admm_kkt_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace qp {
namespace {

bool luFactor(
    std::vector<double>& matrix,
    int n,
    std::vector<int>& pivot
) {
    pivot.resize(static_cast<std::size_t>(n));

    double scale = 0.0;

    for (double value : matrix) {
        if (!std::isfinite(value)) {
            return false;
        }

        scale = std::max(scale, std::abs(value));
    }

    if (n == 0) {
        return true;
    }

    const double tolerance =
        1e-14 * std::max(1.0, scale);

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;

        double pivotAbs =
            std::abs(
                matrix[
                    static_cast<std::size_t>(k) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(k)
                ]
            );

        for (int i = k + 1; i < n; ++i) {
            const double candidate =
                std::abs(
                    matrix[
                        static_cast<std::size_t>(i) *
                            static_cast<std::size_t>(n) +
                        static_cast<std::size_t>(k)
                    ]
                );

            if (candidate > pivotAbs) {
                pivotAbs = candidate;
                pivotRow = i;
            }
        }

        if (!std::isfinite(pivotAbs) ||
            pivotAbs <= tolerance) {
            return false;
        }

        pivot[static_cast<std::size_t>(k)] = pivotRow;

        if (pivotRow != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(
                    matrix[
                        static_cast<std::size_t>(k) *
                            static_cast<std::size_t>(n) +
                        static_cast<std::size_t>(j)
                    ],
                    matrix[
                        static_cast<std::size_t>(pivotRow) *
                            static_cast<std::size_t>(n) +
                        static_cast<std::size_t>(j)
                    ]
                );
            }
        }

        const double diagonal =
            matrix[
                static_cast<std::size_t>(k) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(k)
            ];

        if (!std::isfinite(diagonal) ||
            std::abs(diagonal) <= tolerance) {
            return false;
        }

        for (int i = k + 1; i < n; ++i) {
            const std::size_t ik =
                static_cast<std::size_t>(i) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(k);

            matrix[ik] /= diagonal;

            if (!std::isfinite(matrix[ik])) {
                return false;
            }

            for (int j = k + 1; j < n; ++j) {
                const std::size_t ij =
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j);

                const std::size_t kj =
                    static_cast<std::size_t>(k) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j);

                matrix[ij] -=
                    matrix[ik] * matrix[kj];

                if (!std::isfinite(matrix[ij])) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool luSolve(
    const std::vector<double>& lu,
    int n,
    const std::vector<int>& pivot,
    std::vector<double>& rhs
) {
    if (static_cast<int>(rhs.size()) != n ||
        static_cast<int>(pivot.size()) != n) {
        return false;
    }

    /*
     * Apply row permutations generated during
     * partial pivoting.
     */
    for (int k = 0; k < n; ++k) {
        const int pivotRow =
            pivot[static_cast<std::size_t>(k)];

        if (pivotRow < 0 || pivotRow >= n) {
            return false;
        }

        if (pivotRow != k) {
            std::swap(
                rhs[static_cast<std::size_t>(k)],
                rhs[static_cast<std::size_t>(pivotRow)]
            );
        }
    }

    /*
     * Forward substitution for L.
     */
    for (int i = 0; i < n; ++i) {
        double sum =
            rhs[static_cast<std::size_t>(i)];

        for (int j = 0; j < i; ++j) {
            sum -=
                lu[
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j)
                ] *
                rhs[static_cast<std::size_t>(j)];
        }

        rhs[static_cast<std::size_t>(i)] = sum;

        if (!std::isfinite(sum)) {
            return false;
        }
    }

    /*
     * Back substitution for U.
     */
    for (int i = n - 1; i >= 0; --i) {
        double sum =
            rhs[static_cast<std::size_t>(i)];

        for (int j = i + 1; j < n; ++j) {
            sum -=
                lu[
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j)
                ] *
                rhs[static_cast<std::size_t>(j)];
        }

        const double diagonal =
            lu[
                static_cast<std::size_t>(i) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(i)
            ];

        if (!std::isfinite(diagonal) ||
            diagonal == 0.0) {
            return false;
        }

        rhs[static_cast<std::size_t>(i)] =
            sum / diagonal;

        if (!std::isfinite(
                rhs[static_cast<std::size_t>(i)])) {
            return false;
        }
    }

    return true;
}

}  // namespace

SuperAdmmKktSolver::SuperAdmmKktSolver(
    const QpModel& model
)
    : model_(&model),
      n_(model.numVariables()),
      m_(model.numConstraints()) {
    model.validate();
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

    if (!std::isfinite(sigma) ||
        sigma <= 0.0) {
        return false;
    }

    if (static_cast<int>(rho.size()) != m_ ||
        static_cast<int>(rhsX.size()) != n_ ||
        static_cast<int>(rhsNu.size()) != m_) {
        return false;
    }

    for (double value : rho) {
        if (!std::isfinite(value) ||
            value <= 0.0) {
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

    /*
     * KKT dimension:
     *
     *     N = n + m
     */
    if (n_ > std::numeric_limits<int>::max() - m_) {
        return false;
    }

    const int N = n_ + m_;

    if (N == 0) {
        x.clear();
        nu.clear();
        return true;
    }

    /*
     * Guard against impossible / overflowing dense allocations.
     */
    const std::size_t sizeN =
        static_cast<std::size_t>(N);

    if (sizeN >
        std::numeric_limits<std::size_t>::max() / sizeN) {
        return false;
    }

    const std::size_t matrixSize =
        sizeN * sizeN;

    if (matrixSize >
        std::vector<double>().max_size()) {
        return false;
    }

    std::vector<double> K(
        matrixSize,
        0.0
    );

    /*
     * ----------------------------------------------------
     * Top-left block: P + sigma I
     * ----------------------------------------------------
     */
    const auto& pRows =
        model_->P.csrRowStart();

    const auto& pCols =
        model_->P.csrColumnIndex();

    const auto& pVals =
        model_->P.csrValues();

    for (int i = 0; i < n_; ++i) {
        const Offset begin =
            pRows[
                static_cast<std::size_t>(i)
            ];

        const Offset end =
            pRows[
                static_cast<std::size_t>(i + 1)
            ];

        for (Offset k = begin; k < end; ++k) {
            const int j =
                pCols[
                    static_cast<std::size_t>(k)
                ];

            if (j < 0 || j >= n_) {
                return false;
            }

            K[
                static_cast<std::size_t>(i) *
                    sizeN +
                static_cast<std::size_t>(j)
            ] +=
                pVals[
                    static_cast<std::size_t>(k)
                ];
        }

        K[
            static_cast<std::size_t>(i) *
                sizeN +
            static_cast<std::size_t>(i)
        ] += sigma;
    }

    /*
     * ----------------------------------------------------
     * Off-diagonal blocks:
     *
     *     [ A^T ]
     *     [ A   ]
     *
     * Bottom-right block:
     *
     *     -R^{-1}
     * ----------------------------------------------------
     */
    const auto& aRows =
        model_->A.csrRowStart();

    const auto& aCols =
        model_->A.csrColumnIndex();

    const auto& aVals =
        model_->A.csrValues();

    for (int i = 0; i < m_; ++i) {
        const Offset begin =
            aRows[
                static_cast<std::size_t>(i)
            ];

        const Offset end =
            aRows[
                static_cast<std::size_t>(i + 1)
            ];

        for (Offset k = begin; k < end; ++k) {
            const Index j =
                aCols[
                    static_cast<std::size_t>(k)
                ];

            if (j < 0 || j >= n_) {
                return false;
            }

            const double value =
                aVals[
                    static_cast<std::size_t>(k)
                ];

            const std::size_t rowA =
                static_cast<std::size_t>(n_ + i);

            const std::size_t colX =
                static_cast<std::size_t>(j);

            K[
                colX * sizeN +
                rowA
            ] = value;

            K[
                rowA * sizeN +
                colX
            ] = value;
        }

        K[
            static_cast<std::size_t>(n_ + i) *
                sizeN +
            static_cast<std::size_t>(n_ + i)
        ] =
            -1.0 /
            rho[
                static_cast<std::size_t>(i)
            ];
    }

    /*
     * Combined RHS:
     *
     * [ rhsX ]
     * [ rhsNu ]
     */
    std::vector<double> rhs(
        static_cast<std::size_t>(N),
        0.0
    );

    for (int i = 0; i < n_; ++i) {
        rhs[
            static_cast<std::size_t>(i)
        ] =
            rhsX[
                static_cast<std::size_t>(i)
            ];
    }

    for (int i = 0; i < m_; ++i) {
        rhs[
            static_cast<std::size_t>(n_ + i)
        ] =
            rhsNu[
                static_cast<std::size_t>(i)
            ];
    }

    std::vector<int> pivot;

    if (!luFactor(K, N, pivot)) {
        return false;
    }

    if (!luSolve(K, N, pivot, rhs)) {
        return false;
    }

    x.resize(
        static_cast<std::size_t>(n_)
    );

    nu.resize(
        static_cast<std::size_t>(m_)
    );

    for (int i = 0; i < n_; ++i) {
        x[
            static_cast<std::size_t>(i)
        ] =
            rhs[
                static_cast<std::size_t>(i)
            ];
    }

    for (int i = 0; i < m_; ++i) {
        nu[
            static_cast<std::size_t>(i)
        ] =
            rhs[
                static_cast<std::size_t>(n_ + i)
            ];
    }

    return true;
}

}  // namespace qp