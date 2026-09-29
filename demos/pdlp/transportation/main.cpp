#include "../../common/data_utils.h"
#include "../../common/demo_utils.h"
#include "../../common/reporting.h"

#include "solver/orchestrator.h"

#include <iostream>

int main(int argc, char* argv[]) {
    demos::printHeader("PDLP - Transportation");

    if (argc < 2) {
        std::cerr
            << "Usage: transportation_pdlp "
               "<path-to-transportation.mps>\n";
        return 2;
    }

    try {
        auto loaded = demos::loadMps(argv[1]);

        demos::printModelSummary(loaded.model);

        if (!demos::validateModel(
                loaded.model,
                "PDLP transportation model")) {
            return 1;
        }

        if (demos::countIntegerVariables(loaded.model) != 0) {
            std::cerr
                << "ERROR: PDLP demo requires a continuous LP. "
                   "Integer variables were detected.\n";
            return 1;
        }

        if (!loaded.model.objective.quadraticTerms.empty()) {
            std::cerr
                << "ERROR: PDLP demo requires a linear objective. "
                   "Quadratic terms were detected.\n";
            return 1;
        }

        solver::SolverOptions options;
        options.forceEngine = solver::Engine::Pdlp;

        const auto result =
            solver::solve(loaded.model, options);

        demos::printSolveResult(
            result,
            "PDLP Transportation Result");

        return demos::solveSucceeded(result) ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << '\n';
        return 1;
    }
}
