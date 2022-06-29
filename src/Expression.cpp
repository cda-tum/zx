#include "Expression.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>

namespace sym {
    std::unordered_map<std::string, std::size_t> Variable::registered = std::unordered_map<std::string, std::size_t>();
    std::unordered_map<std::size_t, std::string> Variable::names      = std::unordered_map<std::size_t, std::string>();
    std::size_t                                  Variable::nextId;

    Variable::Variable(const std::string& name) {
        auto it = registered.find(name);
        if (it != registered.end()) {
            id = it->second;
        } else {
            registered[name] = nextId;
            names[nextId]    = name;
            id               = nextId;
            nextId++;
        }
    }

    std::string Variable::getName() const {
        return names[id];
    }

    // bool Expression::isPauli() const {
    //     return isConstant() && constant.isInteger();
    // }
    // bool Expression::isClifford() const {
    //     return isConstant() && (constant.isInteger() || constant.getDenom() == 2);
    // }
    // bool Expression::isProperClifford() const {
    //     return isConstant() && constant.getDenom() == 2;
    // }
    // void Expression::roundToClifford(fp tolerance) {
    //     if (!isConstant())
    //         return;

    //     if (constant.isCloseDivPi(0, tolerance)) {
    //         constant = PiRational(0, 1);
    //     } else if (constant.isCloseDivPi(0.5, tolerance)) {
    //         constant = PiRational(1, 2);
    //     } else if (constant.isCloseDivPi(-0.5, tolerance)) {
    //         constant = PiRational(-1, 2);
    //     } else if (constant.isCloseDivPi(1, tolerance)) {
    //         constant = PiRational(1, 1);
    //     }
    // }

} // namespace sym
