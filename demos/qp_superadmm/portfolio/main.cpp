#include "../../common/demo_utils.h"
#include "../../common/qp3_loader.h"
#include "../../common/reporting.h"

#include "qp/qp_adapter.h"
#include "qp/super_admm_solver.h"

#include <iostream>

int main(int argc, char* argv[]) {
    demos::printHeader(
        "SuperADMM - Portfolio Optimization - MINLPLib qp3");

    if (argc < 2) {
        std::cerr
            << "Usage: portfolio_superadmm <path-to-qp3.lp>\n";
        return 2;
    }

    try {
        const model::Model model =
            demos::loadQp3Lp(argv[1]);

        demos::printModelSummary(model);

        if (!demos::validateModel(
                model,
                "SuperADMM portfolio model")) {
            return 1;
        }

        if (demos::countIntegerVariables(model) != 0) {
            std::cerr
                << "ERROR: SuperADMM demo requires a continuous "
                   "quadratic model. Integer variables were detected.\n";
            return 1;
        }

        if (model.objective.quadraticTerms.empty()) {
            std::cerr
                << "ERROR: expected a quadratic objective.\n";
            return 1;
        }

        qp::QpTranslation translation;

        const qp::QpModel problem =
            qp::fromModel(model, translation);

        qp::SuperAdmmOptions options;
        options.iterationLimit = 10000;
        options.primalTolerance = 1e-8;
        options.dualTolerance = 1e-8;
        options.alpha = 500.0;
        options.sigma = 1e-6;
        options.b0 = 1e8;
        options.tau = 0.5;
        options.rho0 = 1.0;

        const auto result =
            qp::SuperAdmmSolver(problem, options).solve();

        demos::printQpResult(
            result,
            "SuperADMM Portfolio Result");

        return demos::solveSucceeded(result) ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << '\n';
        return 1;
    }
}
