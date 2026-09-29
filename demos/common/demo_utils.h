#pragma once

#include "model/model.h"

#include <cstddef>
#include <string>

namespace demos {

void printHeader(const std::string& title);

void printModelSummary(const model::Model& model);

bool validateModel(
    const model::Model& model,
    const std::string& context = {}
);

std::size_t countNonzeros(const model::Model& model);

std::size_t countIntegerVariables(const model::Model& model);

std::size_t countBinaryVariables(const model::Model& model);

std::size_t countContinuousVariables(const model::Model& model);

} // namespace demos
