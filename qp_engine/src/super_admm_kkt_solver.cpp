#include "qp/super_admm_kkt_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace qp {
namespace {

bool luFactor(
    std::vector<double>& a,
    int n,
    std::vector<int>& pivot
) {
    pivot.resize(static_cast<std::size_t>(n));

    double scale = 0.0;

    for (double value : a) {
        scale = std::max(
            scale,
            std::abs(value)
        );
    }

    const double tolerance =
        1e-14 * std::max(1.0, scale);

    for (int k = 0; k < n; ++k) {
        int pivotRow = k;

        double pivotAbs =
            std::abs(
                a[
                    static_cast<std::size_t>(k) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(k)
                ]
            );

        for (int i = k + 1; i < n; ++i) {
            const double candidate =
                std::abs(
                    a[
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

        pivot[
            static_cast<std::size_t>(k)
        ] = pivotRow;

        if (pivotRow != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(
                    a[
                        static_cast<std::size_t>(k) *
                            static_cast<std::size_t>(n) +
                        static_cast<std::size_t>(j)
                    ],
                    a[
                        static_cast<std::size_t>(pivotRow) *
                            static_cast<std::size_t>(n) +
                        static_cast<std::size_t>(j)
                    ]
                );
            }
        }

        const double diagonal =
            a[
                static_cast<std::size_t>(k) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(k)
            ];

        for (int i = k + 1; i < n; ++i) {
            const std::size_t ik =
                static_cast<std::size_t>(i) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(k);

            a[ik] /= diagonal;

            for (int j = k + 1; j < n; ++j) {
                a[
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j)
                ] -=
                    a[ik] *
                    a[
                        static_cast<std::size_t>(k) *
                            static_cast<std::size_t>(n) +
                        static_cast<std::size_t>(j)
                    ];
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
    for (int k = 0; k < n; ++k) {
        const int p =
            pivot[
                static_cast<std::size_t>(k)
            ];

        if (p != k) {
            std::swap(
                rhs[
                    static_cast<std::size_t>(k)
                ],
                rhs[
                    static_cast<std::size_t>(p)
                ]
            );
        }
    }

    // Forward substitution: L y = P b.
    for (int i = 0; i < n; ++i) {
        double sum =
            rhs[
                static_cast<std::size_t>(i)
            ];

        for (int j = 0; j < i; ++j) {
            sum -=
                lu[
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j)
                ] *
                rhs[
                    static_cast<std::size_t>(j)
                ];
        }

        rhs[
            static_cast<std::size_t>(i)
        ] = sum;
    }

    // Back substitution: U x = y.
    for (int i = n - 1; i >= 0; --i) {
        double sum =
            rhs[
                static_cast<std::size_t>(i)
            ];

        for (int j = i + 1; j < n; ++j) {
            sum -=
                lu[
                    static_cast<std::size_t>(i) *
                        static_cast<std::size_t>(n) +
                    static_cast<std::size_t>(j)
                ] *
                rhs[
                    static_cast<std::size_t>(j)
                ];
        }

        const double diagonal =
            lu[
                static_cast<std::size_t>(i) *
                    static_cast<std::size_t>(n) +
                static_cast<std::size_t>(i)
            ];

        if (diagonal == 0.0 ||
            !std::isfinite(diagonal)) {
            return false;
        }

        rhs[
            static_cast<std::size_t>(i)
        ] =
            sum / diagonal;
    }

    for (double value : rhs) {
        if (!std::isfinite(value)) {
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
    if (model_ == nullptr ||
        sigma <= 0.0 ||
        !std::isfinite(sigma)) {
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

    const int N = n_ + m_;

    if (N == 0) {
        x.clear();
        nu.clear();
        return true;
    }

    /*
     * Dense KKT matrix:
     *
     *     [ P + sigma I     A^T ]
     * K = [                    ]
     *     [ A              -R^-1]
     */
    std::vector<double> K(
        static_cast<std::size_t>(N) *
            static_cast<std::size_t>(N),
        0.0
    );

    const auto& pRows =
        model_->P.csrRowStart();

    const auto& pCols =
        model_->P.csrColumnIndex();

    const auto& pVals =
        model_->P.csrValues();

    // P + sigma I.
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
            const Index j =
                pCols[
                    static_cast<std::size_t>(k)
                ];

            K[
                static_cast<std::size_t>(i) *
                    static_cast<std::size_t>(N) +
                static_cast<std::size_t>(j)
            ] +=
                pVals[
                    static_cast<std::size_t>(k)
                ];
        }

        K[
            static_cast<std::size_t>(i) *
                static_cast<std::size_t>(N) +
            static_cast<std::size_t>(i)
        ] += sigma;
    }

    const auto& aRows =
        model_->A.csrRowStart();

    const auto& aCols =
        model_->A.csrColumnIndex();

    const auto& aVals =
        model_->A.csrValues();

    /*
     * Constraint blocks:
     *
     *     [ A^T     ]
     *     [ A  -R^-1]
     */
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

            const double value =
                aVals[
                    static_cast<std::size_t>(k)
                ];

            K[
                static_cast<std::size_t>(j) *
                    static_cast<std::size_t>(N) +
                static_cast<std::size_t>(n_ + i)
            ] = value;

            K[
                static_cast<std::size_t>(n_ + i) *
                    static_cast<std::size_t>(N) +
                static_cast<std::size_t>(j)
            ] = value;
        }

        K[
            static_cast<std::size_t>(n_ + i) *
                static_cast<std::size_t>(N) +
            static_cast<std::size_t>(n_ + i)
        ] =
            -1.0 /
            rho[
                static_cast<std::size_t>(i)
            ];
    }

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