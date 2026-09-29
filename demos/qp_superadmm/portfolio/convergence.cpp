#include <fstream>
#include <string>
#include <vector>

namespace demos {

void writeConvergenceCsv(
    const std::string& filepath,
    const std::vector<double>& primalResidual,
    const std::vector<double>& dualResidual) {

    std::ofstream out(filepath);

    if (!out) {
        return;
    }

    out << "iteration,primal_residual,dual_residual\n";

    const std::size_t count =
        primalResidual.size() < dualResidual.size()
            ? primalResidual.size()
            : dualResidual.size();

    for (std::size_t i = 0; i < count; ++i) {
        out << (i + 1) << ','
            << primalResidual[i] << ','
            << dualResidual[i] << '\n';
    }
}

} // namespace demos
