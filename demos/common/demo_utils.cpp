#include "demo_utils.h"

#include <iomanip>
#include <iostream>

namespace demos {

void printHeader(const std::string& title) {
    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << title << "\n";
    std::cout << "============================================================\n";
}

std::size_t countNonzeros(const model::Model& model) {
    std::size_t count = 0;

    for (const auto& constraint : model.constraints) {
        count += constraint.linearTerms.size();
    }

    count += model.objective.linearTerms.size();
    count += model.objective.quadraticTerms.size();

    return count;
}

std::size_t countIntegerVariables(const model::Model& model) {
    std::size_t count = 0;

    for (const auto& variable : model.variables) {
        if (variable.type == model::VariableType::Integer) {
            ++count;
        }
    }

    return count;
}

std::size_t countBinaryVariables(const model::Model& model) {
    std::size_t count = 0;

    for (const auto& variable : model.variables) {
        if (variable.type == model::VariableType::Binary) {
            ++count;
        }
    }

    return count;
}

std::size_t countContinuousVariables(const model::Model& model) {
    std::size_t count = 0;

    for (const auto& variable : model.variables) {
        if (variable.type == model::VariableType::Continuous) {
            ++count;
        }
    }

    return count;
}

void printModelSummary(const model::Model& model) {
    std::cout << "\nModel summary\n";
    std::cout << "-------------\n";

    std::cout << "Name:                 " << model.name << '\n';
    std::cout << "Variables:            " << model.variables.size() << '\n';
    std::cout << "Constraints:          " << model.constraints.size() << '\n';
    std::cout << "Continuous variables: " << countContinuousVariables(model) << '\n';
    std::cout << "Integer variables:    " << countIntegerVariables(model) << '\n';
    std::cout << "Binary variables:     " << countBinaryVariables(model) << '\n';
    std::cout << "Objective linear terms: "
              << model.objective.linearTerms.size() << '\n';
    std::cout << "Objective quadratic terms: "
              << model.objective.quadraticTerms.size() << '\n';
    std::cout << "Approx. nonzeros:     " << countNonzeros(model) << '\n';

    std::cout << "Objective sense:      "
              << (model.objective.sense == model::ObjectiveSense::Minimize
                      ? "minimize"
                      : "maximize")
              << '\n';

    std::cout << std::setprecision(12);
    std::cout << "Objective offset:     "
              << model.objective.offset << '\n';
}

bool validateModel(
    const model::Model& model,
    const std::string& context
) {
    const bool valid = model.validate();

    if (valid) {
        return true;
    }

    std::cerr << "\nERROR: model validation failed";

    if (!context.empty()) {
        std::cerr << " [" << context << "]";
    }

    std::cerr << "\n";

    return false;
}

} // namespace demos
