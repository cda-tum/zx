#include "Expression.hpp"
#include "Rational.hpp"
#include "ZXDiagram.hpp"

#include <gtest/gtest.h>

class ExpressionTest: public ::testing::Test {
public:
    // const sym::Variable x_var{0, "x"};
    // const sym::Variable y_var{1, "y"};
    // const sym::Variable z_var{2, "z"};

    sym::Term x{sym::Variable("x")};
    sym::Term y{sym::Variable("y")};
    sym::Term z{sym::Variable("z")};

protected:
    virtual void SetUp() {}
};

TEST_F(ExpressionTest, basic_ops_1) {
    sym::Expression<double> e(x);

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
    sym::Expression<zx::PiRational> e1;
    e1 += x;
    e1 += 10.0 * y;
    e1 += 5.0 * z;
    e1 += zx::PiRational(1, 2);

    sym::Expression<zx::PiRational> e2;
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
