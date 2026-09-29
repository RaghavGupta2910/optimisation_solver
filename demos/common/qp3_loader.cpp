#include "qp3_loader.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <regex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace demos {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

struct ParsedTerm {
    int variable = -1;
    double coefficient = 0.0;
};

std::string trim(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");

    if (first == std::string::npos) {
        return {};
    }

    const auto last = value.find_last_not_of(" \t\r\n");

    return value.substr(first, last - first + 1);
}

int variableIndex(
    const std::string& name,
    std::unordered_map<std::string, int>& indices,
    model::Model& model) {

    const auto it = indices.find(name);

    if (it != indices.end()) {
        return it->second;
    }

    const int index =
        static_cast<int>(model.variables.size());

    model::Variable variable;
    variable.name = name;
    variable.type = model::VariableType::Continuous;
    variable.lowerBound = 0.0;
    variable.upperBound = kInf;

    model.variables.push_back(variable);
    indices.emplace(name, index);

    return index;
}

/*
 * Parse linear LP terms of the form:
 *
 *   2.5 x2
 *   - 2.5 x2
 *   x2
 *   - x2
 *   + x2
 *
 * Scientific notation is supported.
 */
std::vector<ParsedTerm> parseLinearExpression(
    const std::string& expression,
    std::unordered_map<std::string, int>& indices,
    model::Model& model) {

    /*
     * A term consists of:
     *
     *   optional sign
     *   optional coefficient
     *   variable x<number>
     *
     * The spaces in expressions such as
     *
     *   - 2.815e-2 x2
     *
     * are deliberately allowed.
     */
    static const std::regex termRegex(
        R"(([+-])?\s*((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?\s*(x\d+))");

    std::vector<ParsedTerm> terms;

    auto begin = std::sregex_iterator(
        expression.begin(),
        expression.end(),
        termRegex);

    const auto end = std::sregex_iterator();

    for (auto it = begin; it != end; ++it) {
        const std::smatch& match = *it;

        double coefficient = 1.0;

        const std::string sign = match[1].str();
        const std::string coefficientText = match[2].str();
        const std::string variableName = match[3].str();

        if (!coefficientText.empty()) {
            coefficient = std::stod(coefficientText);
        }

        if (sign == "-") {
            coefficient = -coefficient;
        }

        const int index =
            variableIndex(variableName, indices, model);

        terms.push_back({index, coefficient});
    }

    return terms;
}

void addLinearTerms(
    const std::vector<ParsedTerm>& terms,
    model::Constraint& constraint) {

    for (const auto& term : terms) {
        constraint.linearTerms.push_back(
            {term.variable, term.coefficient});
    }
}

/*
 * Parse quadratic terms from the qp3 objective.
 *
 * Examples:
 *
 *   2 x2 * x52
 *   x10 * x60
 *   - 2.5 x20 * x70
 */
void parseQuadraticExpression(
    const std::string& expression,
    std::unordered_map<std::string, int>& indices,
    model::Model& model) {

    static const std::regex quadraticRegex(
        R"(([+-])?\s*((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?\s*(x\d+)\s*\*\s*(x\d+))");

    auto begin = std::sregex_iterator(
        expression.begin(),
        expression.end(),
        quadraticRegex);

    const auto end = std::sregex_iterator();

    for (auto it = begin; it != end; ++it) {
        const std::smatch& match = *it;

        double coefficient = 1.0;

        const std::string sign = match[1].str();
        const std::string coefficientText = match[2].str();
        const std::string leftVariable = match[3].str();
        const std::string rightVariable = match[4].str();

        if (!coefficientText.empty()) {
            coefficient = std::stod(coefficientText);
        }

        if (sign == "-") {
            coefficient = -coefficient;
        }

        const int i =
            variableIndex(leftVariable, indices, model);

        const int j =
            variableIndex(rightVariable, indices, model);

        model.objective.quadraticTerms.push_back(
            {i, j, coefficient});
    }
}

} // namespace

model::Model loadQp3Lp(const std::string& filepath) {

    std::ifstream input(filepath);

    if (!input) {
        throw std::runtime_error(
            "Unable to open qp3 LP file: " + filepath);
    }

    model::Model model;
    model.name = "qp3";
    model.objective.sense =
        model::ObjectiveSense::Minimize;
    model.objective.offset = 0.0;

    std::unordered_map<std::string, int> variableIndices;

    /*
     * qp3 contains x2 ... x101.
     *
     * Create all 100 variables up front so variables which only
     * appear in the objective/bounds are retained.
     */
    for (int number = 2; number <= 101; ++number) {
        const std::string name =
            "x" + std::to_string(number);

        variableIndex(name, variableIndices, model);
    }

    std::string line;

    bool inObjective = false;
    bool inConstraints = false;
    bool inBounds = false;

    std::string objectiveText;

    while (std::getline(input, line)) {

        line = trim(line);

        if (line.empty() || line[0] == '\\') {
            continue;
        }

        if (line == "Minimize") {
            inObjective = true;
            inConstraints = false;
            inBounds = false;
            continue;
        }

        if (line == "Subject To") {
            inObjective = false;
            inConstraints = true;
            inBounds = false;
            continue;
        }

        if (line == "Bounds") {
            inObjective = false;
            inConstraints = false;
            inBounds = true;
            continue;
        }

        if (line == "End") {
            break;
        }

        if (inObjective) {
            objectiveText += " ";
            objectiveText += line;
            continue;
        }

        if (inBounds) {

            /*
             * qp3 uses lines such as:
             *
             *   x52 Free
             */
            const std::size_t space =
                line.find_first_of(" \t");

            if (space == std::string::npos) {
                continue;
            }

            const std::string variableName =
                trim(line.substr(0, space));

            const std::string boundType =
                trim(line.substr(space));

            if (variableName.empty()) {
                continue;
            }

            const int index =
                variableIndex(
                    variableName,
                    variableIndices,
                    model);

            if (boundType == "Free") {
                model.variables[index].lowerBound = -kInf;
                model.variables[index].upperBound = kInf;
            }

            continue;
        }

        if (!inConstraints) {
            continue;
        }

        /*
         * Each qp3 constraint starts with a name such as:
         *
         *   e2:
         *
         * and may continue over multiple lines.
         */
        const auto colon = line.find(':');

        if (colon == std::string::npos) {
            continue;
        }

        const std::string constraintName =
            trim(line.substr(0, colon));

        std::string fullConstraint =
            trim(line.substr(colon + 1));

        /*
         * Keep reading continuation lines until a relation
         * followed by a numeric RHS is encountered.
         */
        static const std::regex relationRegex(
            R"((<=|>=|=)\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\s*$)");

        while (!std::regex_search(
            fullConstraint,
            relationRegex)) {

            if (!std::getline(input, line)) {
                break;
            }

            line = trim(line);

            if (!line.empty() && line[0] != '\\') {
                fullConstraint += " ";
                fullConstraint += line;
            }
        }

        std::smatch relationMatch;

        if (!std::regex_search(
            fullConstraint,
            relationMatch,
            relationRegex)) {

            throw std::runtime_error(
                "Could not determine relation for constraint " +
                constraintName);
        }

        const std::string relation =
            relationMatch[1].str();

        const std::size_t relationPosition =
            relationMatch.position(1);

        const std::string lhs =
            trim(fullConstraint.substr(
                0,
                relationPosition));

        const std::string rhs =
            trim(fullConstraint.substr(
                relationPosition + relation.size()));

        model::Constraint constraint;
        constraint.name = constraintName;

        addLinearTerms(
            parseLinearExpression(
                lhs,
                variableIndices,
                model),
            constraint);

        const double rhsValue =
            std::stod(rhs);

        if (relation == "=") {
            constraint.lowerBound = rhsValue;
            constraint.upperBound = rhsValue;
        }
        else if (relation == ">=") {
            constraint.lowerBound = rhsValue;
            constraint.upperBound = kInf;
        }
        else {
            constraint.lowerBound = -kInf;
            constraint.upperBound = rhsValue;
        }

        model.constraints.push_back(
            std::move(constraint));
    }

    /*
     * qp3 objective format:
     *
     *   linear_part + [ quadratic_part ] / 2
     *
     * The model's QP convention is:
     *
     *   0.5 x^T P x + q^T x
     *
     * Therefore the coefficients appearing inside the
     * official [ ... ] block are stored directly as
     * QuadraticTerm values.
     */
    const auto bracketBegin =
        objectiveText.find('[');

    const auto bracketEnd =
        objectiveText.find(']');

    if (bracketBegin == std::string::npos ||
        bracketEnd == std::string::npos ||
        bracketEnd <= bracketBegin) {

        throw std::runtime_error(
            "qp3 objective does not contain the expected quadratic block");
    }

    const std::string linearPart =
        objectiveText.substr(
            0,
            bracketBegin);

    const std::string quadraticPart =
        objectiveText.substr(
            bracketBegin + 1,
            bracketEnd - bracketBegin - 1);

    const auto linearTerms =
        parseLinearExpression(
            linearPart,
            variableIndices,
            model);

    for (const auto& term : linearTerms) {
        model.objective.linearTerms.push_back(
            {term.variable, term.coefficient});
    }

    parseQuadraticExpression(
        quadraticPart,
        variableIndices,
        model);

    if (!model.validate()) {
        throw std::runtime_error(
            "Parsed qp3 model failed model::Model::validate()");
    }

    return model;
}

} // namespace demos
