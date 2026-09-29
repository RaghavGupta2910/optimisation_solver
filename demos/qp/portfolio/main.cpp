#include "../../common/demo_utils.h"
#include "../../common/qp3_loader.h"
#include "../../common/reporting.h"

#include "qp/qp_adapter.h"
#include "qp/qp_solver.h"

#include <iostream>

int main(int argc, char* argv[]) {
    demos::printHeader("QP - Portfolio Optimization - MINLPLib qp3");

    if (argc < 2) {
        std::cerr
            << "Usage: portfolio_qp <path-to-qp3.lp>\n";
        return 2;
    }

    try {
        const model::Model model =
            demos::loadQp3Lp(argv[1]);

        demos::printModelSummary(model);

        if (!demos::validateModel(
                model,
                "QP portfolio model")) {
            return 1;
        }

        if (demos::countIntegerVariables(model) != 0) {
            std::cerr
                << "ERROR: QP demo requires a continuous model. "
                   "Integer variables were detected.\n";
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

        qp::AdmmOptions options;
        options.iterationLimit = 10000;
        options.primalTolerance = 1e-8;
        options.dualTolerance = 1e-8;

        const auto result =
            qp::QpSolver{}.solve(problem, options);

        demos::printQpResult(
            result,
            "QP Portfolio Result");

        return demos::solveSucceeded(result) ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << '\n';
        return 1;
    }
}
