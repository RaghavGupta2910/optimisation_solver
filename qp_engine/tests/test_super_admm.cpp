#include "qp/qp_model.h"
#include "qp/super_admm_solver.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using qp::AdmmResult;
using qp::QpModel;
using qp::QpStatus;
using qp::SparseMatrix;
using qp::SuperAdmmOptions;
using qp::SuperAdmmSolver;

constexpr double kTestTolerance = 1e-4;

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
            << " | expected = " << expected
            << ", actual = " << actual
            << '\n';

        std::exit(EXIT_FAILURE);
    }
}

SparseMatrix makeMatrix(
    int rows,
    int columns,
    const std::vector<int>& row,
    const std::vector<int>& column,
    const std::vector<double>& values
) {
    std::vector<double> rowDouble(
        row.begin(),
        row.end()
    );

    std::vector<double> columnDouble(
        column.begin(),
        column.end()
    );

    return SparseMatrix::fromTriplets(
        rows,
        columns,
        rowDouble,
        columnDouble,
        values
    );
}

QpModel makeModel(
    int n,
    int m,
    const std::vector<int>& pRow,
    const std::vector<int>& pCol,
    const std::vector<double>& pVal,
    const std::vector<int>& aRow,
    const std::vector<int>& aCol,
    const std::vector<double>& aVal,
    const std::vector<double>& q,
    const std::vector<double>& l,
    const std::vector<double>& u
) {
    QpModel model;

    model.P = makeMatrix(
        n,
        n,
        pRow,
        pCol,
        pVal
    );

    model.A = makeMatrix(
        m,
        n,
        aRow,
        aCol,
        aVal
    );

    model.q = q;
    model.l = l;
    model.u = u;

    return model;
}

SuperAdmmOptions defaultOptions() {
    SuperAdmmOptions options;

    options.iterationLimit = 5000;
    options.timeLimitSeconds = 0.0;

    options.primalTolerance = 1e-6;
    options.dualTolerance = 1e-6;

    options.alpha = 500.0;
    options.sigma = 1e-6;
    options.b0 = 1e8;
    options.tau = 0.5;
    options.rho0 = 1.0;

    options.infeasibilityCheckInterval = 10;
    options.infeasibilityTolerance = 1e-8;

    return options;
}

/*
 * --------------------------------------------------------------------------
 * 1. Simple bound-constrained QP
 *
 *     minimize 1/2 x^2
 *     subject to 1 <= x <= 2
 *
 * Expected solution:
 *
 *     x = 1
 * --------------------------------------------------------------------------
 */
void testSimpleBound() {
    QpModel model = makeModel(
        1,
        1,

        // P = [1]
        {0},
        {0},
        {1.0},

        // A = [1]
        {0},
        {0},
        {1.0},

        // q
        {0.0},

        // l <= Ax <= u
        {1.0},
        {2.0}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "simple bound should converge"
    );

    requireNear(
        result.primal[0],
        1.0,
        kTestTolerance,
        "simple bound solution"
    );
}

/*
 * --------------------------------------------------------------------------
 * 2. Equality constraint
 *
 *     minimize 1/2(x1^2 + x2^2)
 *     subject to x1 + x2 = 1
 *
 * Expected:
 *
 *     x1 = x2 = 0.5
 * --------------------------------------------------------------------------
 */
void testEqualityConstraint() {
    QpModel model = makeModel(
        2,
        1,

        // P = I
        {0, 1},
        {0, 1},
        {1.0, 1.0},

        // A = [1 1]
        {0, 0},
        {0, 1},
        {1.0, 1.0},

        {0.0, 0.0},

        {1.0},
        {1.0}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "equality constraint should converge"
    );

    requireNear(
        result.primal[0],
        0.5,
        kTestTolerance,
        "equality constraint x1"
    );

    requireNear(
        result.primal[1],
        0.5,
        kTestTolerance,
        "equality constraint x2"
    );
}

/*
 * --------------------------------------------------------------------------
 * 3. Mixed inequalities
 *
 *     minimize 1/2(x1^2 + x2^2)
 *
 *     subject to
 *
 *         x1 + x2 >= 1
 *         x1 >= 0
 *
 * Expected solution:
 *
 *     x1 = x2 = 0.5
 * --------------------------------------------------------------------------
 */
void testMixedInequalities() {
    QpModel model = makeModel(
        2,
        2,

        {0, 1},
        {0, 1},
        {1.0, 1.0},

        // A:
        // [1 1]
        // [1 0]
        {0, 0, 1},
        {0, 1, 0},
        {1.0, 1.0, 1.0},

        {0.0, 0.0},

        // x1 + x2 >= 1
        // x1 >= 0
        {1.0, 0.0},

        {
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()
        }
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "mixed inequalities should converge"
    );

    requireNear(
        result.primal[0],
        0.5,
        kTestTolerance,
        "mixed inequalities x1"
    );

    requireNear(
        result.primal[1],
        0.5,
        kTestTolerance,
        "mixed inequalities x2"
    );
}

/*
 * --------------------------------------------------------------------------
 * 4. Positive-semidefinite objective
 *
 *     P = diag(1, 0)
 *
 *     minimize 1/2 x1^2
 *
 *     subject to
 *
 *         x1 + x2 = 1
 *         x2 >= 0
 *
 * Expected:
 *
 *     x1 = 0
 *     x2 = 1
 *
 * This verifies that SuperADMM accepts PSD, not only positive-definite,
 * quadratic matrices.
 * --------------------------------------------------------------------------
 */
void testPsdObjective() {
    QpModel model = makeModel(
        2,
        2,

        // P = diag(1, 0)
        {0},
        {0},
        {1.0},

        // A:
        // [1 1]
        // [0 1]
        {0, 0, 1},
        {0, 1, 1},
        {1.0, 1.0, 1.0},

        {0.0, 0.0},

        {1.0, 0.0},

        {1.0,
         std::numeric_limits<double>::infinity()}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "PSD objective should converge"
    );

    requireNear(
        result.primal[0],
        0.0,
        kTestTolerance,
        "PSD objective x1"
    );

    requireNear(
        result.primal[1],
        1.0,
        kTestTolerance,
        "PSD objective x2"
    );
}

/*
 * --------------------------------------------------------------------------
 * 5. Unconstrained convex QP
 *
 *     minimize 1/2 x^2 - x
 *
 * The solution is x = 1.
 *
 * This verifies the m = 0 path.
 * --------------------------------------------------------------------------
 */
void testUnconstrainedQp() {
    QpModel model = makeModel(
        1,
        0,

        // P = [1]
        {0},
        {0},
        {1.0},

        // A has zero rows
        {},
        {},
        {},

        // q = [-1]
        {-1.0},

        {},
        {}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "unconstrained QP should converge"
    );

    requireNear(
        result.primal[0],
        1.0,
        kTestTolerance,
        "unconstrained QP solution"
    );
}

/*
 * --------------------------------------------------------------------------
 * 6. Infeasible QP
 *
 *     x >= 2
 *     x <= 1
 *
 * The model itself is valid because every individual bound satisfies
 * l <= u; the contradiction comes from two separate constraints.
 * --------------------------------------------------------------------------
 */
void testInfeasibleProblem() {
    QpModel model = makeModel(
        1,
        2,

        {0},
        {0},
        {1.0},

        // A:
        // [1]
        // [1]
        {0, 1},
        {0, 0},
        {1.0, 1.0},

        {0.0},

        {2.0, -std::numeric_limits<double>::infinity()},

        {std::numeric_limits<double>::infinity(), 1.0}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Infeasible,
        "infeasible QP should be detected"
    );
}

/*
 * --------------------------------------------------------------------------
 * 7. Unbounded QP
 *
 *     minimize -x
 *
 *     P = 0
 *     q = -1
 *
 * There are no constraints.
 *
 * The objective is unbounded below.
 * --------------------------------------------------------------------------
 */
void testUnboundedProblem() {
    QpModel model = makeModel(
        1,
        0,

        // P = [0]
        {},
        {},
        {},

        {},
        {},
        {},

        {-1.0},

        {},
        {}
    );

    SuperAdmmOptions options =
        defaultOptions();

    options.infeasibilityCheckInterval =
        10;

    SuperAdmmSolver solver(
        model,
        options
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Unbounded,
        "unbounded QP should be detected"
    );
}

/*
 * --------------------------------------------------------------------------
 * 8. Zero-variable feasible problem
 *
 * There is one constraint:
 *
 *     0 <= 0 <= 1
 *
 * Since x has zero variables, x is necessarily empty and Ax = 0.
 *
 * This should be feasible.
 * --------------------------------------------------------------------------
 */
void testZeroVariableFeasible() {
    QpModel model = makeModel(
        0,
        1,

        {},
        {},
        {},

        {},
        {},
        {},

        {},

        {0.0},
        {1.0}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "zero-variable feasible problem should be optimal"
    );

    require(
        result.primal.empty(),
        "zero-variable solution should be empty"
    );
}

/*
 * --------------------------------------------------------------------------
 * 9. Zero-variable infeasible problem
 *
 * x has zero variables, so Ax = 0.
 *
 * Constraint:
 *
 *     1 <= 0
 *
 * is impossible.
 *
 * This specifically checks the zero-variable bug that previously caused
 * SuperADMM to incorrectly return Optimal.
 * --------------------------------------------------------------------------
 */
void testZeroVariableInfeasible() {
    QpModel model = makeModel(
        0,
        1,

        {},
        {},
        {},

        {},
        {},
        {},

        {},

        {1.0},
        {std::numeric_limits<double>::infinity()}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Infeasible,
        "zero-variable infeasible problem should be detected"
    );
}

/*
 * --------------------------------------------------------------------------
 * 10. Indefinite quadratic matrix
 *
 *     P = [-1]
 *
 * This is not a convex QP and therefore must be rejected by SuperADMM.
 * --------------------------------------------------------------------------
 */
void testIndefiniteP() {
    QpModel model = makeModel(
        1,
        0,

        {0},
        {0},
        {-1.0},

        {},
        {},
        {},

        {0.0},

        {},
        {}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::InvalidProblem,
        "indefinite P should be rejected"
    );
}

/*
 * --------------------------------------------------------------------------
 * 11. Non-symmetric quadratic matrix
 *
 *     P = [1  2]
 *         [0  1]
 *
 * This must be rejected instead of silently symmetrizing the problem.
 * --------------------------------------------------------------------------
 */
void testNonSymmetricP() {
    QpModel model = makeModel(
        2,
        0,

        {0, 0, 1},
        {0, 1, 1},
        {1.0, 2.0, 1.0},

        {},
        {},
        {},

        {0.0, 0.0},

        {},
        {}
    );

    SuperAdmmSolver solver(
        model,
        defaultOptions()
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::InvalidProblem,
        "non-symmetric P should be rejected"
    );
}

/*
 * --------------------------------------------------------------------------
 * 12. Iteration limit
 *
 * Force the solver to stop before convergence.
 * --------------------------------------------------------------------------
 */
void testIterationLimit() {
    QpModel model = makeModel(
        1,
        1,

        {0},
        {0},
        {1.0},

        {0},
        {0},
        {1.0},

        {0.0},

        {1.0},
        {2.0}
    );

    SuperAdmmOptions options =
        defaultOptions();

    options.iterationLimit = 1;

    /*
     * Use very tight tolerances so that the test cannot accidentally
     * classify the first iteration as converged.
     */
    options.primalTolerance = 1e-14;
    options.dualTolerance = 1e-14;

    SuperAdmmSolver solver(
        model,
        options
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::IterationLimit,
        "iteration limit should be reported"
    );

    require(
        result.iterations == 1,
        "iteration count should equal the iteration limit"
    );
}

/*
 * --------------------------------------------------------------------------
 * 13. Time limit
 *
 * A practically immediate time limit should terminate the solve.
 * --------------------------------------------------------------------------
 */
void testTimeLimit() {
    QpModel model = makeModel(
        1,
        1,

        {0},
        {0},
        {1.0},

        {0},
        {0},
        {1.0},

        {0.0},

        {1.0},
        {2.0}
    );

    SuperAdmmOptions options =
        defaultOptions();

    options.timeLimitSeconds =
        1e-12;

    SuperAdmmSolver solver(
        model,
        options
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::TimeLimit,
        "time limit should be reported"
    );
}

/*
 * --------------------------------------------------------------------------
 * 14. Dynamic weighting
 *
 * Constraint:
 *
 *     0 <= x <= 1
 *
 * Objective:
 *
 *     minimize 1/2 x^2
 *
 * The initial z^0 = 0 is exactly on the lower bound, so the constraint
 * is classified as active when the SuperADMM weight update is performed.
 *
 * Therefore rho should increase from its initial value of 1.
 * --------------------------------------------------------------------------
 */
void testDynamicWeighting() {
    QpModel model = makeModel(
        1,
        1,

        {0},
        {0},
        {1.0},

        {0},
        {0},
        {1.0},

        {0.0},

        {0.0},
        {1.0}
    );

    SuperAdmmOptions options =
        defaultOptions();

    options.iterationLimit =
        50;

    SuperAdmmSolver solver(
        model,
        options
    );

    const AdmmResult result =
        solver.solve();

    require(
        result.status == QpStatus::Optimal,
        "dynamic-weighting problem should converge"
    );

    require(
        result.finalRho > options.rho0,
        "active constraint should increase its SuperADMM weight"
    );

    requireNear(
        result.primal[0],
        0.0,
        kTestTolerance,
        "dynamic-weighting solution"
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

    testPsdObjective();
    std::cout << "PASS: PSD objective\n";

    testUnconstrainedQp();
    std::cout << "PASS: unconstrained QP\n";

    testInfeasibleProblem();
    std::cout << "PASS: infeasible problem\n";

    testUnboundedProblem();
    std::cout << "PASS: unbounded problem\n";

    testZeroVariableFeasible();
    std::cout << "PASS: zero-variable feasible problem\n";

    testZeroVariableInfeasible();
    std::cout << "PASS: zero-variable infeasible problem\n";

    testIndefiniteP();
    std::cout << "PASS: indefinite P rejection\n";

    testNonSymmetricP();
    std::cout << "PASS: non-symmetric P rejection\n";

    testIterationLimit();
    std::cout << "PASS: iteration limit\n";

    testTimeLimit();
    std::cout << "PASS: time limit\n";

    testDynamicWeighting();
    std::cout << "PASS: dynamic weighting\n";

    std::cout << "\nAll SuperADMM tests passed.\n";

    return EXIT_SUCCESS;
}