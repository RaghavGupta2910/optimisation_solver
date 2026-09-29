#pragma once

#include "pdlp/pdlp_result.h"
#include "qp/admm_solver.h"
#include "solver/solve_result.h"

#include <string>

namespace demos {

void printSolveResult(
    const solver::SolveResult& result,
    const std::string& title = "Solve Result");

void printPdlpResult(
    const pdlp::PdlpResult& result,
    const std::string& title = "PDLP Result");

void printQpResult(
    const qp::AdmmResult& result,
    const std::string& title = "QP Result");

bool solveSucceeded(const solver::SolveResult& result);

bool solveSucceeded(const pdlp::PdlpResult& result);

bool solveSucceeded(const qp::AdmmResult& result);

} // namespace demos
