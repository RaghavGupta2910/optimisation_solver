#include "../../common/demo_utils.h"

#include <iostream>

int main() {
    demos::printHeader(
        "NLP - Reactor Optimization - MINLPLib ex8_3_13");

    std::cout
        << "Dataset: MINLPLib ex8_3_13\n"
        << "Problem class: NLP\n\n"
        << "Status: NLP backend is not currently implemented "
           "in the repository solver stack.\n\n"
        << "The nonlinear model is therefore not silently "
           "relaxed or converted into another problem class.\n";

    return 1;
}
