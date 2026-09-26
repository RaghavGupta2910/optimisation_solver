#include "qp/super_admm_solver.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using qp::QpModel;
using qp::QpStatus;
using qp::SparseMatrix;
using qp::SuperAdmmOptions;
using qp::SuperAdmmSolver;

SparseMatrix makeSparse(
    int rows,
    int columns,
    const std::vector<int>& row,
    const std::vector<int>& column,
    const std::vector<double>& value
) {
    std::vector<double> rowDouble(row.begin(), row.end());
    std::vector<double> columnDouble(column.begin(), column.end());

    return SparseMatrix::fromTriplets(
        rows,
        columns,
        rowDouble,
        columnDouble,
        value
    );
}

void require(
    bool condition,
    const std::string& message
) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

void requireNear(
    double actual,
    double expected,
    double tolerance,
    const std::string& message
) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr
            << "FAIL: " << message
            << " expected=" << expected
            << " actual=" << actual
            << " tolerance=" << tolerance
            << '\n';

        std::exit(EXIT_FAILURE);
    }
}

// -----------------------------------------------------------------------------
// Test 1
//
// min 0.5 x^2
// s.t. 1 <= x <= 2
//
// Solution: x = 1.
// -----------------------------------------------------------------------------
void testSimpleBound() {
    QpModel model;

    model.P = makeSparse(
        1,
        1,
        {0},
        {0},
        {1.0}
    );

    model.A = makeSparse(
        1,
        1,
        {0},
        {0},
        {1.0}
    );

    model.q = {0.0};
    model.l = {1.0};
    model.u = {2.0};

    SuperAdmmOptions options;
    options.iterationLimit = 100;
    options.primalTolerance = 1e-5;
    options.dualTolerance = 1e-5;

    SuperAdmmSolver solver(model, options);
    const auto result = solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "simple bound problem should be optimal"
    );

    requireNear(
        result.primal[0],
        1.0,
        2e-5,
        "simple bound solution"
    );
}

// -----------------------------------------------------------------------------
// Test 2
//
// min 0.5(x1^2 + x2^2)
// s.t. x1 + x2 = 1
//
// Solution: x1 = x2 = 0.5.
// -----------------------------------------------------------------------------
void testEqualityConstraint() {
    QpModel model;

    model.P = makeSparse(
        2,
        2,
        {0, 0, 1, 1},
        {0, 1, 0, 1},
        {1.0, 0.0, 0.0, 1.0}
    );

    model.A = makeSparse(
        1,
        2,
        {0, 0},
        {0, 1},
        {1.0, 1.0}
    );

    model.q = {0.0, 0.0};
    model.l = {1.0};
    model.u = {1.0};

    SuperAdmmOptions options;
    options.iterationLimit = 100;
    options.primalTolerance = 1e-6;
    options.dualTolerance = 1e-6;

    SuperAdmmSolver solver(model, options);
    const auto result = solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "equality constrained problem should be optimal"
    );

    requireNear(
        result.primal[0],
        0.5,
        2e-5,
        "equality solution x1"
    );

    requireNear(
        result.primal[1],
        0.5,
        2e-5,
        "equality solution x2"
    );
}

// -----------------------------------------------------------------------------
// Test 3
//
// min 0.5(x1^2 + x2^2)
// s.t.
//     x1 + x2 >= 1
//     x1 >= 0
//
// Solution: x1 = x2 = 0.5.
// -----------------------------------------------------------------------------
void testMixedInequalities() {
    QpModel model;

    const double inf =
        std::numeric_limits<double>::infinity();

    model.P = makeSparse(
        2,
        2,
        {0, 0, 1, 1},
        {0, 1, 0, 1},
        {1.0, 0.0, 0.0, 1.0}
    );

    model.A = makeSparse(
        2,
        2,
        {0, 0, 1},
        {0, 1, 0},
        {1.0, 1.0, 1.0}
    );

    model.q = {0.0, 0.0};

    model.l = {
        1.0,
        0.0
    };

    model.u = {
        inf,
        inf
    };

    SuperAdmmOptions options;
    options.iterationLimit = 100;
    options.primalTolerance = 1e-6;
    options.dualTolerance = 1e-6;

    SuperAdmmSolver solver(model, options);
    const auto result = solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "mixed inequality problem should be optimal"
    );

    requireNear(
        result.primal[0],
        0.5,
        2e-5,
        "mixed inequality x1"
    );

    requireNear(
        result.primal[1],
        0.5,
        2e-5,
        "mixed inequality x2"
    );
}

// -----------------------------------------------------------------------------
// Test 4
//
// PSD objective:
//
// min 0.5*x1^2
// s.t.
//     x1 + x2 = 1
//     x2 >= 0
//
// Solution: x1 = 0, x2 = 1.
// -----------------------------------------------------------------------------
void testPositiveSemidefiniteObjective() {
    QpModel model;

    const double inf =
        std::numeric_limits<double>::infinity();

    model.P = makeSparse(
        2,
        2,
        {0},
        {0},
        {1.0}
    );

    model.A = makeSparse(
        2,
        2,
        {0, 0, 1},
        {0, 1, 1},
        {1.0, 1.0, 1.0}
    );

    model.q = {0.0, 0.0};

    model.l = {
        1.0,
        0.0
    };

    model.u = {
        1.0,
        inf
    };

    SuperAdmmOptions options;
    options.iterationLimit = 100;
    options.primalTolerance = 1e-6;
    options.dualTolerance = 1e-6;

    SuperAdmmSolver solver(model, options);
    const auto result = solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "PSD objective problem should be optimal"
    );

    requireNear(
        result.primal[0],
        0.0,
        2e-5,
        "PSD problem x1"
    );

    requireNear(
        result.primal[1],
        1.0,
        2e-5,
        "PSD problem x2"
    );
}

// -----------------------------------------------------------------------------
// Test 5
//
// Infeasible:
//
//     x >= 2
//     x <= 1
//
// Expected status: Infeasible.
// -----------------------------------------------------------------------------
void testInfeasibleProblem() {
    const double inf =
        std::numeric_limits<double>::infinity();

    QpModel model;

    model.P = makeSparse(
        1,
        1,
        {0},
        {0},
        {1.0}
    );

    model.A = makeSparse(
        2,
        1,
        {0, 1},
        {0, 0},
        {1.0, 1.0}
    );

    model.q = {0.0};

    model.l = {
        2.0,
        -inf
    };

    model.u = {
        inf,
        1.0
    };

    SuperAdmmOptions options;
    options.iterationLimit = 100;
    options.primalTolerance = 1e-6;
    options.dualTolerance = 1e-6;

    SuperAdmmSolver solver(model, options);
    const auto result = solver.solve();

    require(
        result.status == QpStatus::Infeasible,
        "infeasible problem should receive Infeasible status"
    );
}

}  // namespace

int main() {
    testSimpleBound();
    std::cout << "PASS: simple bound\n";

    testEqualityConstraint();
    std::cout << "PASS: equality constraint\n";

    testMixedInequalities();
    std::cout << "PASS: mixed inequalities\n";

    testPositiveSemidefiniteObjective();
    std::cout << "PASS: PSD objective\n";

    testInfeasibleProblem();
    std::cout << "PASS: infeasible problem\n";

    std::cout << "\nAll SuperADMM tests passed.\n";

    return EXIT_SUCCESS;
}