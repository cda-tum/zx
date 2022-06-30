#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sym {
    static constexpr double TOLERANCE = 1e-9;

    class SymbolicException: public std::invalid_argument {
        std::string msg;

    public:
        explicit SymbolicException(std::string msg):
            std::invalid_argument("Symbolic Exception"), msg(std::move(msg)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return msg.c_str();
        }
    };

    struct Variable {
        static std::unordered_map<std::string, std::size_t> registered;
        static std::unordered_map<std::size_t, std::string> names;
        static std::size_t                                  nextId;

        explicit Variable(const std::string& name);

        [[nodiscard]] std::string getName() const;

        inline bool operator==(const Variable& rhs) const {
            return id == rhs.id;
        }

        inline bool operator!=(const Variable& rhs) const {
            return !((*this) == rhs);
        }

        inline bool operator<(const Variable& rhs) const {
            return id < rhs.id;
        }

        inline bool operator>(const Variable& rhs) const {
            return id > rhs.id;
        }

    private:
        std::size_t id{};
    };
} // namespace sym

namespace std {
    template<>
    struct hash<sym::Variable> {
        std::size_t operator()(const sym::Variable& var) const {
            return std::hash<std::string>()(var.getName());
        }
    };
} // namespace std

namespace sym {
    using VariableAssignment = std::unordered_map<Variable, double>;

    template<typename T, typename = std::enable_if<std::is_constructible<int, T>::value && std::is_constructible<T, double>::value>>
    class Term {
    public:
        [[nodiscard]] Variable getVar() const { return var; }
        [[nodiscard]] T        getCoeff() const { return coeff; }
        [[nodiscard]] bool     hasZeroCoeff() const {
                return std::abs(static_cast<double>(coeff)) < TOLERANCE;
        }

        Term(T coeff, Variable var):
            coeff(coeff), var(var){};
        explicit Term(Variable var):
            coeff(1), var(var){};

        Term operator-() const { return Term(-coeff, var); }

        void addCoeff(const T& r) {
            coeff += r;
        }
        Term& operator*=(const T& rhs) {
            coeff *= rhs;
            return *this;
        }

        Term& operator/=(const T& rhs) {
            coeff /= rhs;
            return *this;
        }
        [[nodiscard]] bool totalAssignment(const VariableAssignment& assignment) const {
            return assignment.find(getVar()) != assignment.end();
        }

        [[nodiscard]] double evaluate(const VariableAssignment& assignment) const {
            if (!totalAssignment(assignment))
                throw SymbolicException("Cannot instantiate variable " + getVar().getName() + ". No value given.");
            return assignment.at(getVar()) * getCoeff();
        }

    private:
        T        coeff;
        Variable var;
    };
    template<typename T, typename = std::enable_if<std::is_constructible<int, T>::value>>
    inline Term<T> operator*(Term<T> lhs, double rhs) {
        lhs *= rhs;
        return lhs;
    }
    template<typename T, typename = std::enable_if<std::is_constructible<int, T>::value>>
    inline Term<T> operator/(Term<T> lhs, double rhs) {
        lhs /= rhs;
        return lhs;
    }
    template<typename T, typename = std::enable_if<std::is_constructible<int, T>::value>>
    inline Term<T> operator*(double lhs, const Term<T>& rhs) {
        return rhs * lhs;
    }
    template<typename T, typename = std::enable_if<std::is_constructible<int, T>::value>>
    inline Term<T> operator/(double lhs, const Term<T>& rhs) {
        return rhs / lhs;
    }

    template<typename T, typename U, typename = std::enable_if<std::is_constructible<T, U>::value && std::is_constructible<U, T>::value && std::is_constructible<int, T>::value && std::is_constructible<T, double>::value && std::is_constructible<U, double>::value>>
    class Expression {
    public:
        using iterator       = typename std::vector<Term<T>>::iterator;
        using const_iterator = typename std::vector<Term<T>>::const_iterator;

        template<typename... Args>
        explicit Expression(Term<T> t, Args... ms) {
            terms.emplace_back(t);
            (terms.emplace_back(std::forward<Args>(ms)), ...);
            sortTerms();
            aggregateEqualTerms();
        }

        template<typename... Args>
        explicit Expression(Variable v, Args... ms) {
            terms.emplace_back(Term(T{1}, v));
            (terms.emplace_back(std::forward<Args>(ms)), ...);
            sortTerms();
            aggregateEqualTerms();
        }

        Expression(const std::vector<Term<T>>& terms, const U& constant):
            terms(terms), constant(constant){};

        Expression() = default;

        explicit Expression(const U& r):
            constant(r){};

        iterator                     begin() { return terms.begin(); }
        iterator                     end() { return terms.end(); }
        [[nodiscard]] const_iterator begin() const { return terms.cbegin(); }
        [[nodiscard]] const_iterator end() const { return terms.cend(); }
        [[nodiscard]] const_iterator cbegin() const { return terms.cbegin(); }
        [[nodiscard]] const_iterator cend() const { return terms.cend(); }

        [[nodiscard]] bool isZero() const { return terms.empty() && constant == U{T{0}}; }
        [[nodiscard]] bool isConstant() const { return terms.empty(); }

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
                        terms.begin(), terms.end(), *t, [&](const Term<T>& lhs, const Term<T>& rhs) {
                            return lhs.getVar() < rhs.getVar();
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

        Expression<T, U>& operator+=(const Term<T>& rhs) {
            return *this += Expression(rhs);
        }

        Expression<T, U>& operator+=(const U& rhs) {
            constant += rhs;
            return *this;
        }

        Expression<T, U>& operator-=(const Expression<T, U>& rhs) {
            return *this += -rhs;
        }

        Expression<T, U>& operator-=(const Term<T>& rhs) {
            return *this += -rhs;
        }
        Expression<T, U>& operator-=(const U& rhs) {
            return *this += -rhs;
        }

        Expression<T, U>& operator*=(const T& rhs) {
            if (std::abs(static_cast<double>(rhs)) < TOLERANCE) {
                terms.clear();
                constant = U{T{0}};
                return *this;
            }
            std::for_each(terms.begin(), terms.end(), [&](auto& term) { term *= rhs; });
            constant *= U{rhs};
            return *this;
        }

        template<typename = std::enable_if<!std::is_same<T, U>::value>>
        Expression<T, U>& operator*=(const U& rhs) {
            if (std::abs(static_cast<double>(T{rhs})) < TOLERANCE) {
                terms.clear();
                constant = U{T{0}};
                return *this;
            }
            std::for_each(terms.begin(), terms.end(), [&](auto& term) { term *= T{rhs}; });
            constant *= rhs;
            return *this;
        }

        Expression<T, U>& operator/=(const T& rhs) {
            if (std::abs(static_cast<double>(T{rhs})) < TOLERANCE) {
                throw std::runtime_error("Trying to divide expression by 0!");
            }
            std::for_each(terms.begin(), terms.end(), [&](auto& term) { term /= rhs; });
            constant /= U{rhs};
            return *this;
        }

        template<typename = std::enable_if<!std::is_same<T, U>::value>>
        Expression<T, U>& operator/=(const U& rhs) {
            if (std::abs(static_cast<double>(T{rhs})) < TOLERANCE) {
                throw std::runtime_error("Trying to divide expression by 0!");
            }
            std::for_each(terms.begin(), terms.end(), [&](auto& term) { term /= T{rhs}; });
            constant /= rhs;
            return *this;
        }

        [[nodiscard]] Expression<T, U> operator-() const {
            Expression<T, U> e;
            e.terms.reserve(terms.size());
            for (auto& t: terms)
                e.terms.push_back(-t);
            e.constant = -constant;
            return e;
        }

        [[nodiscard]] const Term<T>& operator[](std::size_t i) const { return terms[i]; }
        [[nodiscard]] U              getConst() const { return constant; }
        void                         setConst(const U& val) { constant = val; }
        [[nodiscard]] auto           numTerms() const { return terms.size(); }

        [[nodiscard]] const std::vector<Term<T>>& getTerms() const { return terms; }

        template<typename V, typename std::enable_if<std::is_constructible<U, V>::value>::type* = nullptr>
        Expression<T, V> convert() const {
            return Expression<T, V>(terms, V{constant});
        }

        double evaluate(const VariableAssignment& assignment) {
            auto initial = static_cast<double>(constant);
            return std::accumulate(terms.begin(), terms.end(), initial, [&](double sum, const auto& term) { return term.evaluate(assignment) + sum; });
        }

    private:
        std::vector<Term<T>> terms;
        U                    constant{T{0.0}};

        void sortTerms() {
            std::sort(terms.begin(), terms.end(), [&](const Term<T>& lhs, const Term<T>& rhs) {
                return lhs.getVar() < rhs.getVar();
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

    template<typename T, typename U>
    inline Expression<T, U> operator+(Expression<T, U>        lhs,
                                      const Expression<T, U>& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator+(Expression<T, U> lhs, const Term<T>& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator+(const Term<T>& lhs, Expression<T, U> rhs) {
        rhs += lhs;
        return rhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator+(const U& lhs, Expression<T, U> rhs) {
        rhs += lhs;
        return rhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator+(Expression<T, U> lhs, const U& rhs) {
        lhs += rhs;
        return lhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator+(const T& lhs, Expression<T, U> rhs) {
        rhs += rhs;
        return rhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator-(Expression<T, U> lhs, const Expression<T, U>& rhs) {
        lhs -= rhs;
        return lhs;
    }
    template<typename T, typename U>
    inline Expression<T, U> operator-(Expression<T, U> lhs, const Term<T>& rhs) {
        lhs -= rhs;
        return lhs;
    }
    template<typename T, typename U>
    inline Expression<T, U> operator-(const Term<T>& lhs, Expression<T, U> rhs) {
        rhs -= lhs;
        return rhs;
    }
    template<typename T, typename U>
    inline Expression<T, U> operator-(const U& lhs, Expression<T, U> rhs) {
        rhs -= lhs;
        return rhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator-(Expression<T, U> lhs, const U& rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator*(Expression<T, U> lhs, const T& rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename T, typename U, typename std::enable_if<!std::is_same<T, U>::value>::type* = nullptr>
    inline Expression<T, U> operator*(Expression<T, U> lhs, const U& rhs) {
        lhs *= rhs;
        return lhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator/(Expression<T, U> lhs, const T& rhs) {
        lhs /= rhs;
        return lhs;
    }

    template<typename T, typename U, typename std::enable_if<!std::is_same<T, U>::value>::type* = nullptr>
    inline Expression<T, U> operator/(Expression<T, U> lhs, const U& rhs) {
        lhs /= rhs;
        return lhs;
    }

    template<typename T, typename U>
    inline Expression<T, U> operator*(const T& lhs, Expression<T, U> rhs) {
        return rhs * lhs;
    }

    template<typename T, typename U, typename std::enable_if<!std::is_same<T, U>::value>::type* = nullptr>
    inline Expression<T, U> operator*(const U& lhs, Expression<T, U> rhs) {
        return rhs * lhs;
    }

    template<typename T, typename U>
    inline bool operator==(const Expression<T, U>& lhs, const Expression<T, U>& rhs) {
        if (lhs.numTerms() != rhs.numTerms() || lhs.getConst() != rhs.getConst())
            return false;

        for (size_t i = 0; i < lhs.numTerms(); ++i) {
            if (std::abs(lhs[i].getCoeff() - rhs[i].getCoeff()) >= TOLERANCE)
                return false;
        }
        return true;
    }

    std::ostream& operator<<(std::ostream& os, const Variable& var);

    template<typename T>
    std::ostream& operator<<(std::ostream& os, const Term<T>& term) {
        os << term.getCoeff() << "*" << term.getVar().getName();
        return os;
    }

    template<typename T, typename U>
    std::ostream& operator<<(std::ostream& os, const Expression<T, U>& expr) {
        std::for_each(expr.begin(), expr.end(), [&](const auto& term) { os << term << " + "; });
        os << expr.getConst();
        return os;
    }
} // namespace sym
