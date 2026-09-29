#include "data_utils.h"

#include "mps/mps_reader.h"

#include <filesystem>
#include <stdexcept>

namespace demos {

bool fileExists(const std::string& filepath) {
    return std::filesystem::exists(filepath);
}

LoadedModel loadMps(const std::string& filepath) {
    if (!fileExists(filepath)) {
        throw std::runtime_error(
            "Dataset file does not exist: " + filepath
        );
    }

    mps::MpsReader reader;

    LoadedModel loaded;
    loaded.model = reader.read(filepath);
    loaded.warnings = reader.warnings();

    return loaded;
}

} // namespace demos
