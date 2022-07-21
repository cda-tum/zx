#pragma once

#include "Definitions.hpp"
#include "Rational.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace zx {
    struct Variable {
        Variable(std::int32_t id, std::string name):
            id(id), name(std::move(name)){};
        std::int32_t id;
        std::string  name;
    };

    inline bool operator==(const Variable& lhs, const Variable& rhs) {
        return lhs.id == rhs.id;
    }
    class Term {
    public:
        [[nodiscard]] Variable getVar() const { return var; }
        [[nodiscard]] double   getCoeff() const { return coeff; }
        [[nodiscard]] bool     hasZeroCoeff() const { return std::abs(coeff) < TOLERANCE; }

        void addCoeff(double r);
        Term(double coeff, Variable var):
            coeff(coeff), var(std::move(var)){};
        explicit Term(Variable var):
            coeff(1), var(std::move(var)){};

        Term  operator-() const { return Term(-coeff, var); }
        Term& operator*=(double rhs);
        Term& operator/=(double rhs);

    private:
        double   coeff;
        Variable var;
    };

    inline Term operator*(Term lhs, double rhs) {
        lhs *= rhs;
        return lhs;
    }
    inline Term operator/(Term lhs, double rhs) {
        lhs /= rhs;
        return lhs;
    }
    inline Term operator*(double lhs, const Term& rhs) {
        return rhs * lhs;
    }

    inline Term operator/(double lhs, const Term& rhs) {
        return rhs / lhs;
    }

    class Expression {
    public:
        using iterator       = std::vector<Term>::iterator;
        using const_iterator = std::vector<Term>::const_iterator;

        template<typename... Args>
        explicit Expression(Term t, Args... ms) {
            terms.emplace_back(t);
            (terms.emplace_back(std::forward<Args>(ms)), ...);
            sortTerms();
            aggregateEqualTerms();
        }

        template<typename... Args>
        explicit Expression(Variable v, Args... ms) {
            terms.emplace_back(Term(1, std::move(v)));
            (terms.emplace_back(std::forward<Args>(ms)), ...);
            sortTerms();
            aggregateEqualTerms();
        }

        Expression():
            constant(PiRational(0, 1)){};
        explicit Expression(PiRational r):
            constant(std::move(r)){};

        iterator                     begin() { return terms.begin(); }
        iterator                     end() { return terms.end(); }
        [[nodiscard]] const_iterator begin() const { return terms.cbegin(); }
        [[nodiscard]] const_iterator end() const { return terms.cend(); }
        [[nodiscard]] const_iterator cbegin() const { return terms.cbegin(); }
        [[nodiscard]] const_iterator cend() const { return terms.cend(); }

        [[nodiscard]] bool isZero() const;
        [[nodiscard]] bool isConstant() const;
        [[nodiscard]] bool isPauli() const;
        [[nodiscard]] bool isClifford() const;
        [[nodiscard]] bool isProperClifford() const;

        void        roundToClifford(fp tolerance);
        Expression& operator+=(const Expression& rhs);
        Expression& operator+=(const Term& rhs);
        Expression& operator+=(const PiRational& rhs);

        Expression&              operator-=(const Expression& rhs);
        Expression&              operator-=(const Term& rhs);
        Expression&              operator-=(const PiRational& rhs);
        [[nodiscard]] Expression operator-() const;

        [[nodiscard]] const Term& operator[](std::size_t i) const { return terms[i]; }
        [[nodiscard]] PiRational  getConst() const { return constant; }
        [[nodiscard]] auto        numTerms() const { return terms.size(); }

    private:
        std::vector<Term> terms;
        PiRational        constant;
        void              sortTerms();
        void              aggregateEqualTerms();
    };

    inline Expression operator+(Expression lhs, const Expression& rhs) {
        lhs += rhs;
        return lhs;
    }
    inline Expression operator+(Expression lhs, const Term& rhs) {
        lhs += rhs;
        return lhs;
    }
    inline Expression operator+(Expression lhs, const PiRational& rhs) {
        lhs += rhs;
        return lhs;
    }
    inline Expression operator-(Expression lhs, const Expression& rhs) {
        lhs -= rhs;
        return lhs;
    }
    inline Expression operator-(Expression lhs, const Term& rhs) {
        lhs -= rhs;
        return lhs;
    }
    inline Expression operator-(Expression lhs, const PiRational& rhs) {
        lhs -= rhs;
        return lhs;
    }

    bool operator==(const Expression& lhs, const Expression& rhs);
} // namespace zx

inline std::ostream& operator<<(std::ostream& os, const zx::Variable& rhs) {
    os << rhs.name;
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const zx::Term& rhs) {
    os << rhs.getCoeff() << "*" << rhs.getVar();
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const zx::Expression& rhs) {
    for (auto& t: rhs) {
        os << t << " + ";
    }
    os << rhs.getConst();
    return os;
}
