#include "../../common/data_utils.h"
#include "../../common/demo_utils.h"
#include "../../common/reporting.h"

#include "solver/orchestrator.h"

#include <iostream>

int main(int argc, char* argv[]) {
    demos::printHeader("MILP - Supply Chain - MIPLIB shs1023");

    if (argc < 2) {
        std::cerr << "Usage: supply_chain_milp <path-to-shs1023.mps>\n";
        return 2;
    }

    try {
        auto loaded = demos::loadMps(argv[1]);

        demos::printModelSummary(loaded.model);

        if (!demos::validateModel(
                loaded.model,
                "MILP supply-chain model")) {
            return 1;
        }

        if (demos::countIntegerVariables(loaded.model) == 0) {
            std::cerr
                << "ERROR: expected an integer or mixed-integer model.\n";
            return 1;
        }

        solver::SolverOptions options;
        options.forceEngine = solver::Engine::BranchAndCut;

        const auto result =
            solver::solve(loaded.model, options);

        demos::printSolveResult(
            result,
            "MILP Supply Chain Result");

        return demos::solveSucceeded(result) ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << '\n';
        return 1;
    }
}
