#include "Expression.hpp"
#include "Rational.hpp"
#include "ZXDiagram.hpp"

#include <gtest/gtest.h>
#include <stdexcept>

class ExpressionTest: public ::testing::Test {
public:
    // const sym::Variable x_var{0, "x"};
    // const sym::Variable y_var{1, "y"};
    // const sym::Variable z_var{2, "z"};

    sym::Term<double> x{1.0, sym::Variable("x")};
    sym::Term<double> y{sym::Variable("y")};
    sym::Term<double> z{sym::Variable("z")};

protected:
    virtual void SetUp() {}
};

TEST_F(ExpressionTest, basic_ops_1) {
    sym::Expression<double, zx::PiRational> e(x);

    EXPECT_EQ(1, e.numTerms());
    EXPECT_EQ(zx::PiRational(0, 1), e.getConst());

    e += x; // sym::Term(x);

    EXPECT_EQ(1, e.numTerms());
    EXPECT_EQ(zx::PiRational(0, 1), e.getConst());
    EXPECT_PRED_FORMAT2(testing::FloatLE, e[0].getCoeff(), 2.0);

    e += y;
    EXPECT_EQ(2, e.numTerms());
    EXPECT_PRED_FORMAT2(testing::FloatLE, e[0].getCoeff(), 2.0);
    EXPECT_PRED_FORMAT2(testing::FloatLE, e[1].getCoeff(), 1.0);
    EXPECT_EQ(e[0].getVar().getName(), "x");
    EXPECT_EQ(e[1].getVar().getName(), "y");
}

TEST_F(ExpressionTest, basic_ops_2) {
    sym::Expression<double, zx::PiRational> e1;
    e1 += x;
    e1 += 10.0 * y;
    e1 += 5.0 * z;
    e1 += zx::PiRational(1, 2);

    sym::Expression<double, zx::PiRational> e2;
    e2 += -5.0 * x;
    e2 += -10.0 * y;
    e2 += -4.9 * z;
    e2 += zx::PiRational(3, 2);

    auto sum = e1 + e2;

    EXPECT_EQ(2, sum.numTerms());
    EXPECT_PRED_FORMAT2(testing::FloatLE, sum[0].getCoeff(), -4.0);
    EXPECT_PRED_FORMAT2(testing::FloatLE, sum[1].getCoeff(), 0.1);
    EXPECT_EQ(sum[0].getVar().getName(), "x");
    EXPECT_EQ(sum[1].getVar().getName(), "z");
    EXPECT_EQ(sum.getConst(), zx::PiRational(0, 1));
}

TEST_F(ExpressionTest, mult) {
    sym::Expression<double, zx::PiRational> e(x);

    e = e * 2.0;

    EXPECT_PRED_FORMAT2(testing::FloatLE, e[0].getCoeff(), 2);

    e = e * zx::PiRational(0.5);

    EXPECT_PRED_FORMAT2(testing::FloatLE, e[0].getCoeff(), 1);

    e += sym::Expression<double, zx::PiRational>{};

    EXPECT_PRED_FORMAT2(testing::FloatLE, e[0].getCoeff(), 1);

    e = e * 0.0;

    EXPECT_TRUE(e.isZero());
}

TEST_F(ExpressionTest, div) {
    sym::Expression<double, zx::PiRational> e(x);

    e = e / 2.0;

    EXPECT_PRED_FORMAT2(testing::FloatLE, e[0].getCoeff(), 0.5);

    e = e / zx::PiRational(0.5);

    EXPECT_PRED_FORMAT2(testing::FloatLE, e[0].getCoeff(), 1);

    EXPECT_THROW(e = e / 0.0, std::runtime_error);
}

TEST_F(ExpressionTest, Commutativity) {
    sym::Expression<double, double> e1(x, y);
    sym::Expression<double, double> e2(z);
    e2.setConst(1.0);

    EXPECT_EQ(e1 + e2, e2 + e1);
    EXPECT_EQ(e1 * 2.0, 2.0 * e1);
}

TEST_F(ExpressionTest, Associativity) {
    sym::Expression<double, double> e1(x, y);
    sym::Expression<double, double> e2(z);
    sym::Expression<double, double> e3(1.0);

    EXPECT_EQ(e1 + (e2 + e3), (e1 + e2) + e3);
    EXPECT_EQ(e1 * (2.0 * 4.0), (e1 * 2.0) * 4.0);
}

TEST_F(ExpressionTest, Distributive) {
    sym::Expression<double, double> e1(x, y);
    sym::Expression<double, double> e2(z);

    EXPECT_EQ((e1 + e2) * 2.0, (e1 * 2.0) + (e2 * 2.0));
    EXPECT_EQ((e1 - e2) * 2.0, (e1 * 2.0) - (e2 * 2.0));
    std::cout << (e1 + e2) / 2.0 << std::endl;
    std::cout << (e1 / 2.0) + (e2 / 2.0) << std::endl;
    EXPECT_EQ((e1 + e2) / 2.0, (e1 / 2.0) + (e2 / 2.0));
    EXPECT_EQ((e1 - e2) / 2.0, (e1 / 2.0) - (e2 / 2.0));
}

TEST_F(ExpressionTest, Variable) {
    EXPECT_NE(sym::Variable{"x"}, sym::Variable{"y"});
    EXPECT_EQ(sym::Variable{"x"}, sym::Variable{"x"});
    EXPECT_TRUE(sym::Variable{"x"} < sym::Variable{"y"});
    EXPECT_TRUE(sym::Variable{"z"} > sym::Variable{"y"});
}

TEST_F(ExpressionTest, SumNegation) {
    sym::Expression<double, double> e1(x, y);
    sym::Expression<double, double> e2(z, y);

    EXPECT_EQ(e1 - e2, e1 + (-e2));
}

TEST_F(ExpressionTest, SumMult) {
    sym::Expression<double, double> e1(x, y);
    EXPECT_EQ(e1 + e1, e1 * 2.0);
}
