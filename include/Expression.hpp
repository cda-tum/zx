#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sym {
    static constexpr double TOLERANCE = 1e-9;

    struct Variable {
        static std::unordered_map<std::string, std::size_t> registered;
        static std::unordered_map<std::size_t, std::string> names;
        static std::size_t                                  nextId;

        explicit Variable(const std::string& name);
        std::size_t id{};

        [[nodiscard]] std::string getName() const;
    };

    inline bool operator==(const Variable& lhs, const Variable& rhs) {
        return lhs.id == rhs.id;
    }
    class Term {
    public:
        [[nodiscard]] Variable getVar() const { return var; }
        [[nodiscard]] double   getCoeff() const { return coeff; }
        [[nodiscard]] bool     hasZeroCoeff() const {
                return std::abs(coeff) < TOLERANCE;
        }

        void addCoeff(double r);
        Term(double coeff, Variable var):
            coeff(coeff), var(var){};
        explicit Term(Variable var):
            coeff(1), var(var){};

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

    template<class T>
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
            terms.emplace_back(Term(1, v));
            (terms.emplace_back(std::forward<Args>(ms)), ...);
            sortTerms();
            aggregateEqualTerms();
        }

        Expression():
            constant{0.0} {};
        explicit Expression(T r):
            constant{std::move(r)} {};

        iterator                     begin() { return terms.begin(); }
        iterator                     end() { return terms.end(); }
        [[nodiscard]] const_iterator begin() const { return terms.cbegin(); }
        [[nodiscard]] const_iterator end() const { return terms.cend(); }
        [[nodiscard]] const_iterator cbegin() const { return terms.cbegin(); }
        [[nodiscard]] const_iterator cend() const { return terms.cend(); }

        [[nodiscard]] bool isZero() const { return terms.empty() && constant == 0; }
        [[nodiscard]] bool isConstant() const { return terms.empty(); }
        // [[nodiscard]] bool isPauli() const;
        // [[nodiscard]] bool isClifford() const;
        // [[nodiscard]] bool isProperClifford() const;

        // void        roundToClifford(fp tolerance);
        Expression& operator+=(const Expression& rhs) {
            if (this->isZero()) {
                *this = rhs;
                return *this;
            }

            if (rhs.isZero())
                return *this;

            auto t = rhs.begin();

            while (t != rhs.end()) {
                auto insert_pos = std::lower_bound(
                        terms.begin(), terms.end(), *t, [&](const Term& lhs, const Term& rhs) {
                            return lhs.getVar().id < rhs.getVar().id;
                        });
                if (insert_pos != terms.end() && insert_pos->getVar() == t->getVar()) {
                    if (insert_pos->getCoeff() == -t->getCoeff()) {
                        terms.erase(insert_pos);
                    } else {
                        insert_pos->addCoeff(t->getCoeff());
                    }
                } else {
                    terms.insert(insert_pos, *t);
                }
                ++t;
            }
            constant += rhs.constant;
            return *this;
        }

        Expression& operator+=(const Term& rhs) {
            return *this += Expression(rhs);
        }

        Expression& operator+=(const T& rhs) {
            constant += rhs;
            return *this;
        }

        Expression& operator-=(const Expression& rhs) {
            return *this += -rhs;
        }

        Expression& operator-=(const Term& rhs) {
            return *this += -rhs;
        }
        Expression& operator-=(const T& rhs) {
            return *this += -rhs;
        }

        [[nodiscard]] Expression operator-() const {
            Expression e;
            e.terms.reserve(terms.size());
            for (auto& t: terms)
                e.terms.push_back(-t);
            e.constant = -constant;
            return e;
        }

        [[nodiscard]] const Term& operator[](std::size_t i) const { return terms[i]; }
        [[nodiscard]] T           getConst() const { return constant; }
        void                      setConst(const T& val) { constant = val; }
        [[nodiscard]] auto        numTerms() const { return terms.size(); }

    private:
        std::vector<Term> terms;
        T                 constant;

        void sortTerms() {
            std::sort(terms.begin(), terms.end(), [&](const Term& lhs, const Term& rhs) {
                return lhs.getVar().id < rhs.getVar().id;
            });
        }
        void aggregateEqualTerms() {
            for (auto t = terms.begin(); t != terms.end();) {
                auto next = std::next(t);
                while (next != terms.end() && t->getVar() == next->getVar()) {
                    t->addCoeff(next->getCoeff());
                    next = terms.erase(next);
                }
                if (t->hasZeroCoeff()) {
                    t = terms.erase(t);
                } else {
                    t = next;
                }
            }
        }
    };

    template<class T>
    inline Expression<T> operator+(Expression<T> lhs, const Expression<T>& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<class T>
    inline Expression<T> operator+(Expression<T> lhs, const Term& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<class T>
    inline Expression<T> operator+(Expression<T> lhs, const T& rhs) {
        lhs += rhs;
        return lhs;
    }
    template<class T>
    inline Expression<T> operator-(Expression<T> lhs, const Expression<T>& rhs) {
        lhs -= rhs;
        return lhs;
    }
    template<class T>
    inline Expression<T> operator-(Expression<T> lhs, const Term& rhs) {
        lhs -= rhs;
        return lhs;
    }
    template<class T>
    inline Expression<T> operator-(Expression<T> lhs, const T& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<class T>
    inline bool operator==(const Expression<T>& lhs, const Expression<T>& rhs) {
        if (lhs.numTerms() != rhs.numTerms() || lhs.getConst() != rhs.getConst())
            return false;

        for (size_t i = 0; i < lhs.numTerms(); ++i) {
            if (std::abs(lhs[i].getCoeff() - rhs[i].getCoeff()) >= TOLERANCE)
                return false;
        }
        return true;
    }

    inline std::ostream& operator<<(std::ostream& os, const Variable& rhs) {
        os << rhs.getName();
        return os;
    }

    inline std::ostream& operator<<(std::ostream& os, const Term& rhs) {
        // os << rhs.getCoeff() << "*" << rhs.getVar();
        os << rhs.getVar();
        return os;
    }

    template<class T>
    inline std::ostream& operator<<(std::ostream& os, const Expression<T>& rhs) {
        for (auto& t: rhs) {
            os << t << " + ";
        }
        os << rhs.getConst();
        return os;
    }

} // namespace sym
