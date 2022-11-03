#include "Definitions.hpp"
#include "MBQC.hpp"
#include "Rational.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace zx;
using ::testing::_;
using ::testing::ElementsAre;

class MBQCTest: public ::testing::Test {
};

static ZXDiagram makeIdentityDiagram(const std::size_t nqubits,
                                     const std::size_t spidersPerQubit) {
    ZXDiagram           diag(nqubits);
    std::vector<Vertex> rightmostVertices = diag.getInputs();

    for (std::size_t i = 0; i < nqubits; ++i) {
        diag.removeEdge(i, i + nqubits);
    }

    // add identity spiders
    for (Qubit qubit = 0; static_cast<std::size_t>(qubit) < nqubits; ++qubit) {
        for (std::size_t j = 0; j < spidersPerQubit; ++j) {
            const Vertex v = diag.addVertex(qubit);
            diag.addEdge(rightmostVertices[qubit], v);
            rightmostVertices[qubit] = v;
        }
    }

    for (std::size_t qubit = 0; qubit < nqubits; ++qubit) {
        diag.addEdge(rightmostVertices[qubit], qubit + nqubits);
    }

    return diag;
}

static ZXDiagram makeEmptyDiagram(const std::size_t nqubits) {
    auto diag = ::makeIdentityDiagram(nqubits, 0);
    for (std::size_t i = 0; i < nqubits; ++i) {
        diag.removeEdge(i, i + nqubits);
    }
    return diag;
}

TEST_F(MBQCTest, Init) {
    auto diag = makeEmptyDiagram(1);

    const auto& PI_2 = PiExpression(PiRational(1, 2));
    const auto& PI_4 = PiExpression(PiRational(1, 4));

    //add XY Measurement
    diag.addVertex(0, 0); // 2
    diag.addEdge(0, 2);

    //another XY Measurement
    diag.addVertex(0, 0);        // 3
    diag.addVertex(0, -1);       // 4
    diag.addVertex(0, -2, PI_4); // 5
    diag.addHadamardEdge(2, 3);
    diag.addHadamardEdge(3, 4);
    diag.addHadamardEdge(4, 5);

    // YZ Measurement
    diag.addVertex(0, 0);        // 6
    diag.addVertex(0, -1, PI_4); // 7
    diag.addHadamardEdge(3, 6);
    diag.addHadamardEdge(6, 7);

    // XZ Measurement
    diag.addVertex(0, 0, PI_2);  // 8
    diag.addVertex(0, -1, PI_4); // 9
    diag.addHadamardEdge(6, 8);
    diag.addHadamardEdge(8, 9);

    diag.addEdge(8, 1);

    const MBQCPattern pattern(diag);

    const auto& ANGLE = PiRational(1, 4);

    EXPECT_EQ(pattern.getNEdges(), 3);
    EXPECT_EQ(pattern.getNVertices(), 4);
    EXPECT_EQ(pattern[0].plane, MeasurementType::XY);
    EXPECT_EQ(pattern[0].angle, PiRational{});
    EXPECT_EQ(pattern[1].plane, MeasurementType::XY);
    EXPECT_EQ(pattern[1].angle, ANGLE);
    EXPECT_EQ(pattern[2].plane, MeasurementType::YZ);
    EXPECT_EQ(pattern[2].angle, ANGLE);
    EXPECT_EQ(pattern[3].plane, MeasurementType::XZ);
    EXPECT_EQ(pattern[3].angle, ANGLE);

    ASSERT_THAT(pattern.getInputs(), ElementsAre(0));
    ASSERT_THAT(pattern.getOutputs(), ElementsAre(3));
}

TEST_F(MBQCTest, gFlowLine) {
    MBQCPattern pattern({Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}}, {0}, {2}, {{0, 1}, {1, 2}});
    const auto& gFlowOpt = pattern.computeGFlow();
    EXPECT_TRUE(gFlowOpt.has_value());
    const auto& [ord, g] = gFlowOpt.value();

    ASSERT_THAT(g[0], ElementsAre(1));
    ASSERT_THAT(g[1], ElementsAre(2));

    EXPECT_EQ(ord.size(), 3);
    ASSERT_THAT(ord[0], ElementsAre(0));
    ASSERT_THAT(ord[1], ElementsAre(1));
    ASSERT_THAT(ord[2], ElementsAre(2));
}

TEST_F(MBQCTest, gFlow3Qubits) {
    MBQCPattern pattern({Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}}, {0, 1, 2}, {3, 4, 5}, {{0, 3}, {0, 5}, {1, 3}, {1, 4}, {1, 5}, {2, 4}, {2, 5}});
    EXPECT_EQ(pattern.getNEdges(), 7);
    EXPECT_EQ(pattern.getNVertices(), 6);

    const auto& gFlowOpt = pattern.computeGFlow();

    EXPECT_TRUE(gFlowOpt.has_value());
    const auto& [ord, g] = gFlowOpt.value();
    ASSERT_THAT(g[0], ElementsAre(4, 5));
    ASSERT_THAT(g[1], ElementsAre(3, 4, 5));
    ASSERT_THAT(g[2], ElementsAre(3, 5));

    EXPECT_EQ(ord.size(), 2);
    ASSERT_THAT(ord[0], ElementsAre(0, 1, 2));
    ASSERT_THAT(ord[1], ElementsAre(3, 4, 5));
}

TEST_F(MBQCTest, gFlowXY_YZ) {
    MBQCPattern pattern({Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::YZ, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}, Measurement{MeasurementType::XY, PiRational()}}, {0, 4}, {2, 6}, {{0, 1}, {1, 2}, {0, 3}, {3, 4}, {4, 5}, {5, 6}});

    const auto& gFlowOpt = pattern.computeGFlow();

    EXPECT_TRUE(gFlowOpt.has_value());
    // const auto& [opt, g] = gFlowOpt.value();
}
