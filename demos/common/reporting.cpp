#include "reporting.h"

#include <iomanip>
#include <iostream>

namespace demos {

void printSolveResult(
    const solver::SolveResult& result,
    const std::string& title) {

    std::cout << "\n=== " << title << " ===\n";

    std::cout << "Status: "
              << solver::toString(result.status) << '\n';

    std::cout << "Message: "
              << result.message << '\n';

    std::cout << "Selected engine: "
              << solver::toString(result.engine) << '\n';

    std::cout << "Executed engine: "
              << solver::toString(result.executedEngine) << '\n';

    std::cout << std::setprecision(12);

    std::cout << "Objective: "
              << result.objectiveValue << '\n';

    std::cout << "Iterations: "
              << result.iterations << '\n';

    std::cout << "Node count: "
              << result.nodeCount << '\n';

    std::cout << "Solve time: "
              << result.solveSeconds << " s\n";

    std::cout << "Has primal solution: "
              << (result.hasPrimal ? "yes" : "no") << '\n';

    std::cout << "Integrality respected: "
              << (result.integralityRespected ? "yes" : "no") << '\n';
}

void printPdlpResult(
    const pdlp::PdlpResult& result,
    const std::string& title) {

    std::cout << "\n=== " << title << " ===\n";

    std::cout << "Status: "
              << pdlp::toString(result.status) << '\n';

    std::cout << "Status message: "
              << result.statusMessage << '\n';

    std::cout << std::setprecision(12);

    std::cout << "Primal objective: "
              << result.primalObjective << '\n';

    std::cout << "Dual objective: "
              << result.dualObjective << '\n';

    std::cout << "Primal residual: "
              << result.primalResidual << '\n';

    std::cout << "Dual residual: "
              << result.dualResidual << '\n';

    std::cout << "Relative gap: "
              << result.relativeGap << '\n';

    std::cout << "Iterations: "
              << result.iterations << '\n';

    std::cout << "Step trials: "
              << result.stepTrials << '\n';
}

void printQpResult(
    const qp::AdmmResult& result,
    const std::string& title) {

    std::cout << "\n=== " << title << " ===\n";

    std::cout << "Status: "
              << qp::toString(result.status) << '\n';

    std::cout << "Status message: "
              << result.statusMessage << '\n';

    std::cout << std::setprecision(12);

    std::cout << "Primal objective: "
              << result.primalObjective << '\n';

    std::cout << "Dual objective: "
              << result.dualObjective << '\n';

    std::cout << "Best objective: "
              << result.bestObjective << '\n';

    std::cout << "Primal residual: "
              << result.primalResidual << '\n';

    std::cout << "Dual residual: "
              << result.dualResidual << '\n';

    std::cout << "Iterations: "
              << result.iterations << '\n';

    std::cout << "Solve time: "
              << result.solveTimeSeconds << " s\n";

    std::cout << "Final rho: "
              << result.finalRho << '\n';

    std::cout << "Factorizations: "
              << result.factorizations << '\n';
}

bool solveSucceeded(const solver::SolveResult& result) {
    return result.status == solver::SolveStatus::Optimal;
}

bool solveSucceeded(const pdlp::PdlpResult& result) {
    return result.status == pdlp::PdlpStatus::Optimal;
}

bool solveSucceeded(const qp::AdmmResult& result) {
    return result.status == qp::QpStatus::Optimal;
}

} // namespace demos
