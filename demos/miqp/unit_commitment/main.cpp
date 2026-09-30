#include "../../common/data_utils.h"
#include "../../common/demo_utils.h"
#include "../../common/reporting.h"

#include "solver/orchestrator.h"

#include <iostream>

int main(int argc, char* argv[]) {
    demos::printHeader("MIQP - Unit Commitment");

    if (argc < 2) {
        std::cerr
            << "Usage: unit_commitment_miqp "
               "<path-to-unit-commitment-instance>\n";
        return 2;
    }

    try {
        auto loaded = demos::loadMps(argv[1]);

        demos::printModelSummary(loaded.model);

        if (!demos::validateModel(
                loaded.model,
                "MIQP unit-commitment model")) {
            return 1;
        }

        if (demos::countIntegerVariables(loaded.model) == 0) {
            std::cerr
                << "ERROR: expected an integer or mixed-integer model.\n";
            return 1;
        }

        if (loaded.model.objective.quadraticTerms.empty()) {
            std::cerr
                << "ERROR: expected quadratic objective terms.\n";
            return 1;
        }

        const auto result =
            solver::solve(loaded.model);

        demos::printSolveResult(
            result,
            "MIQP Unit Commitment Result");

        if (!demos::solveSucceeded(result)) {
            std::cerr
                << "\nMIQP backend boundary reached. "
                   "The instance was not silently relaxed.\n";
            return 1;
        }

        return 0;
    }
    catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << '\n';
        return 1;
    }
}
