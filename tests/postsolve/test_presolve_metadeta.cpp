#include "model/model.h"
#include "presolve/presolver.h"
#include "test_helpers.h"

#include <iostream>

using namespace test_helpers;

// ============================================================
// Task #3 (presolve side): Fixed variable round trip
//
// Verifies Presolver correctly records a fixed variable in
// PostsolveMetadata.fixedVariables — this is the data Postsolver
// will need to reconstruct x = 5. Does not call Postsolver.
// ============================================================

void test_fixed_variable_recorded() {
    const std::string testName = "test_fixed_variable_recorded";

    model::Model original;
    original.name = "FixedVariableModel";

    // x fixed at 5 via bounds [5, 5]
    original.variables.push_back(makeVariable("x", 5.0, 5.0));
    original.variables.push_back(makeVariable("y", 0.0, 10.0));

    // A constraint referencing both, so removing x must be tracked
    original.constraints.push_back(
        makeConstraint("C1", -100.0, 100.0, { term(0, 1.0), term(1, 1.0) })
    );

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "model is not infeasible", !result.infeasible);

    checkTrue(
        testName,
        "at least one fixed variable was recorded",
        !result.postsolve.fixedVariables.empty()
    );

    bool foundX = false;
    for (const auto& fv : result.postsolve.fixedVariables) {
        if (fv.originalIndex == 0) {
            foundX = true;
            checkNearlyEqual(testName, "fixedValue for x", 5.0, fv.fixedValue);
            checkTrue(testName, "recorded name is 'x'", fv.name == "x");
        }
    }

    checkTrue(testName, "fixed variable record for original index 0 (x) found",
              foundX);

    // Mapping check: original index 0 should map to -1 (eliminated)
    checkTrue(
        testName,
        "originalToPresolvedVar has an entry for x",
        result.postsolve.originalToPresolvedVar.size() > 0
    );

    if (!result.postsolve.originalToPresolvedVar.empty()) {
        checkTrue(
            testName,
            "x (original index 0) maps to -1 (eliminated) in originalToPresolvedVar",
            result.postsolve.originalToPresolvedVar[0] == -1
        );
    }

    std::cout << "[PASS] " << testName << "\n";
}

// ============================================================
// Task #4 (presolve side): Singleton equality
//
// 2x = 10 should presolve to a fixed variable x = 5.
// -2x = -10 should also presolve to x = 5 (sign check).
// ============================================================

void test_singleton_equality_positive_coefficient() {
    const std::string testName = "test_singleton_equality_positive_coefficient";

    model::Model original;
    original.name = "SingletonPositive";
    original.variables.push_back(makeVariable("x", -1000.0, 1000.0));
    original.constraints.push_back(
        makeConstraint("C1", 10.0, 10.0, { term(0, 2.0) })  // 2x = 10
    );

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "model is not infeasible", !result.infeasible);

    bool found = false;
    for (const auto& fv : result.postsolve.fixedVariables) {
        if (fv.originalIndex == 0) {
            found = true;
            checkNearlyEqual(testName, "x from 2x=10", 5.0, fv.fixedValue);
        }
    }
    checkTrue(testName, "x was fixed via singleton row", found);

    std::cout << "[PASS] " << testName << "\n";
}

void test_singleton_equality_negative_coefficient() {
    const std::string testName = "test_singleton_equality_negative_coefficient";

    model::Model original;
    original.name = "SingletonNegative";
    original.variables.push_back(makeVariable("x", -1000.0, 1000.0));
    original.constraints.push_back(
        makeConstraint("C1", -10.0, -10.0, { term(0, -2.0) })  // -2x = -10
    );

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "model is not infeasible", !result.infeasible);

    bool found = false;
    for (const auto& fv : result.postsolve.fixedVariables) {
        if (fv.originalIndex == 0) {
            found = true;
            checkNearlyEqual(
                testName,
                "x from -2x=-10 (sign check)",
                5.0,
                fv.fixedValue
            );
        }
    }
    checkTrue(testName, "x was fixed via singleton row with negative coefficient",
              found);

    std::cout << "[PASS] " << testName << "\n";
}

// ============================================================
// Task #7 (presolve side): Bound tightening
//
// Original: 0 <= x <= 20
// A constraint should tighten this to 5 <= x <= 20.
// Verifies Presolver actually records the tightened bound,
// and that the ORIGINAL bounds are still what the original
// model says (0 to 20) — Postsolve must validate against
// original bounds, not presolved ones.
// ============================================================

void test_bound_tightening_recorded() {
    const std::string testName = "test_bound_tightening_recorded";

    model::Model original;
    original.name = "BoundTighteningModel";
    original.variables.push_back(makeVariable("x", 0.0, 20.0));

    // Constraint that forces x >= 5: e.g. x - y = 5, y in [0,0] fixed... 
    // Simpler: a direct constraint x >= 5 (as a range constraint)
    original.constraints.push_back(
        makeConstraint("C1", 5.0, 20.0, { term(0, 1.0) })
    );

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "model is not infeasible", !result.infeasible);

    // Confirm the ORIGINAL model's bounds are unchanged —
    // this is what Postsolve must validate the final x against.
    checkNearlyEqual(
        testName,
        "original.variables[0].lowerBound unchanged",
        0.0,
        original.variables[0].lowerBound
    );
    checkNearlyEqual(
        testName,
        "original.variables[0].upperBound unchanged",
        20.0,
        original.variables[0].upperBound
    );

    // Look for a TightenLowerBound transformation referencing x
    bool foundTightening = false;
    for (const auto& t : result.transformations) {
        if (
            t.type == presolve::TransformationType::TightenLowerBound &&
            t.originalVariableIndex == 0
        ) {
            foundTightening = true;
            checkNearlyEqual(
                testName,
                "tightened lower bound new value",
                5.0,
                t.newValue
            );
        }
    }

    checkTrue(
        testName,
        "a TightenLowerBound transformation for x was recorded",
        foundTightening
    );

    std::cout << "[PASS] " << testName << "\n";
}

// ============================================================
// Task #26 (presolve side): No-transformation baseline
//
// A model with nothing for Presolve to simplify should produce
// zero transformations and presolved dimensions equal to
// original dimensions.
// ============================================================

void test_no_transformations_baseline() {
    const std::string testName = "test_no_transformations_baseline";

    model::Model original;
    original.name = "NoOpModel";
    original.variables.push_back(makeVariable("x", 1.0, 2.0));
    original.variables.push_back(makeVariable("y", 1.0, 2.0));

    // A constraint that is neither redundant nor a singleton nor tightenable
    original.constraints.push_back(
        makeConstraint("C1", 0.0, 100.0, { term(0, 1.0), term(1, 3.0) })
    );

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "model is not infeasible", !result.infeasible);

    checkTrue(
        testName,
        "presolvedVariables equals originalVariables (no elimination expected)",
        result.presolvedVariables == result.originalVariables
    );

    checkTrue(
        testName,
        "presolvedConstraints equals originalConstraints (no elimination expected)",
        result.presolvedConstraints == result.originalConstraints
    );

    std::cout << "[PASS] " << testName << "\n";
}

int main() {
    test_fixed_variable_recorded();
    test_singleton_equality_positive_coefficient();
    test_singleton_equality_negative_coefficient();
    test_bound_tightening_recorded();
    test_no_transformations_baseline();

    std::cout << "All presolve metadata tests passed successfully!" << std::endl;
    return 0;
}