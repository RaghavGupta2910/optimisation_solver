#include "model/model.h"
#include "presolve/presolver.h"
#include "test_helpers.h"

#include <iostream>

using namespace test_helpers;

// ============================================================
// Task #24: Empty model cases
//
// Uses only Model + Presolver — no Postsolver dependency.
// ============================================================

void test_zero_variables() {
    const std::string testName = "test_zero_variables";

    model::Model original;
    original.name = "EmptyVarsModel";
    // no variables, but keep a trivial constraint with no terms
    original.constraints.push_back(
        makeConstraint("C1", 0.0, 10.0, {})
    );

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "presolve did not crash and returned a result",
              true);
    checkTrue(testName, "model is not reported infeasible",
              !result.infeasible);
    checkTrue(testName, "originalVariables == 0",
              result.originalVariables == 0);

    std::cout << "[PASS] " << testName << "\n";
}

void test_zero_constraints() {
    const std::string testName = "test_zero_constraints";

    model::Model original;
    original.name = "EmptyConstraintsModel";
    original.variables.push_back(makeVariable("x", 0.0, 10.0));

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "presolve did not crash and returned a result",
              true);
    checkTrue(testName, "model is not reported infeasible",
              !result.infeasible);
    checkTrue(testName, "originalConstraints == 0",
              result.originalConstraints == 0);

    std::cout << "[PASS] " << testName << "\n";
}

void test_zero_variables_and_constraints() {
    const std::string testName = "test_zero_variables_and_constraints";

    model::Model original;
    original.name = "CompletelyEmptyModel";

    presolve::Presolver presolver;
    presolve::PresolveResult result = presolver.run(original);

    checkTrue(testName, "presolve did not crash and returned a result",
              true);
    checkTrue(testName, "model is not reported infeasible",
              !result.infeasible);
    checkTrue(testName, "originalVariables == 0",
              result.originalVariables == 0);
    checkTrue(testName, "originalConstraints == 0",
              result.originalConstraints == 0);

    std::cout << "[PASS] " << testName << "\n";
}

int main() {
    test_zero_variables();
    test_zero_constraints();
    test_zero_variables_and_constraints();

    std::cout << "All empty model tests passed successfully!" << std::endl;
    return 0;
}