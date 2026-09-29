#pragma once

#include "model/model.h"

#include <string>
#include <vector>

namespace demos {

struct LoadedModel {
    model::Model model;
    std::vector<std::string> warnings;
};

LoadedModel loadMps(const std::string& filepath);

bool fileExists(const std::string& filepath);

} // namespace demos
