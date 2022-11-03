#include "Definitions.hpp"
#include "Rational.hpp"
#include "Simplify.hpp"
#include "Utils.hpp"
#include "ZXDiagram.hpp"

#include <cstddef>
#include <cstdint>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using ::testing::_;
using ::testing::ElementsAre;
using namespace zx;

class ZXDiagramTest: public ::testing::Test {
public:
    ZXDiagram diag;

    /*
   * Diagram should look like this:
   * {0: (Boundary, Phase=0)} ---H--- {4: (Z, Phase=0)} ---- {5: (Z, Phase = 0)}
   * ----- {1: (Boundary, Phase = 0)}
   *                                                                 |
   *                                                                 |
   *                                                                 |
   * {2: (Boundary, Phase=0)} ------------------------------ {6: (X, Phase = 0)}
   * ----- {3: (Boundary, Phase = 0)}
   */
protected:
    virtual void SetUp() {
        diag = ZXDiagram();
        diag.addQubits(2);
        diag.addVertex(0, 0, PiExpression(), VertexType::Z);
        diag.addEdge(0, 4, EdgeType::Hadamard);
        diag.addVertex(0, 0, PiExpression(), VertexType::Z);
        diag.addEdge(4, 5);
        diag.addVertex(0, 0, PiExpression(), VertexType::X);
        diag.addEdge(2, 6);
        diag.addEdge(5, 6);
        diag.addEdge(5, 1);
        diag.addEdge(6, 3);
    }
};

TEST_F(ZXDiagramTest, create_diagram) {
    EXPECT_EQ(diag.getNVertices(), 7);
    EXPECT_EQ(diag.getNEdges(), 6);

    auto inputs = diag.getInputs();
    EXPECT_EQ(inputs[0], 0);
    EXPECT_EQ(inputs[1], 2);

    auto outputs = diag.getOutputs();
    EXPECT_EQ(outputs[0], 1);
    EXPECT_EQ(outputs[1], 3);

    EXPECT_EQ(diag.getEdge(0, 4).value().type, EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(4, 5).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(2, 6).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 6).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 1).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(6, 3).value().type, EdgeType::Simple);

    EXPECT_EQ(diag.getVData(0).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(1).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(2).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(3).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(4).value().type, VertexType::Z);
    EXPECT_EQ(diag.getVData(5).value().type, VertexType::Z);
    EXPECT_EQ(diag.getVData(6).value().type, VertexType::X);

    const auto nVerts = diag.getNVertices();
    for (std::size_t i = 0; i < nVerts; ++i) {
        EXPECT_TRUE(diag.getVData(6).value().phase.isZero());
    }
}

TEST_F(ZXDiagramTest, deletions) {
    diag.removeVertex(5);
    EXPECT_EQ(diag.getNVertices(), 6);
    EXPECT_EQ(diag.getNEdges(), 3);
    EXPECT_FALSE(diag.getVData(5).has_value());

    diag.removeEdge(0, 4);
    EXPECT_EQ(diag.getNVertices(), 6);
    EXPECT_EQ(diag.getNEdges(), 2);
}

TEST_F(ZXDiagramTest, graph_like) {
    diag.toGraphlike();

    EXPECT_EQ(diag.getEdge(0, 4).value().type, EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(5, 6).value().type, EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(2, 6).value().type, EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(3, 6).value().type, EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(4, 5).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 1).value().type, EdgeType::Simple);

    EXPECT_EQ(diag.getVData(0).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(1).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(2).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(3).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(4).value().type, VertexType::Z);
    EXPECT_EQ(diag.getVData(5).value().type, VertexType::Z);
    EXPECT_EQ(diag.getVData(6).value().type, VertexType::Z);

    const auto nVerts = diag.getNVertices();
    for (std::size_t i = 0; i < nVerts; ++i) {
        EXPECT_TRUE(diag.getVData(i).value().phase.isZero());
    }
}

// TEST_F(ZXDiagramTest, concat) {
//     auto copy = diag;
//     diag.concat(copy);

//     ASSERT_EQ(diag.getNEdges(), 10);
//     ASSERT_EQ(diag.getNVertices(), 10);

//     EXPECT_EQ(diag.getEdge(0, 4).value().type, EdgeType::Hadamard);
//     EXPECT_EQ(diag.getEdge(5, 6).value().type, EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(2, 6).value().type, EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(3, 6).value().type, EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(5, 9).value().type, EdgeType::Hadamard);
//     EXPECT_EQ(diag.getEdge(4, 9).value().type, EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(7, 8).value().type, EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(8, 10).value().type, EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(9, 11).value().type, EdgeType::Simple);

//     EXPECT_EQ(diag.getVData(0).value().type, VertexType::Boundary);
//     EXPECT_EQ(diag.getVData(1).value().type, VertexType::Boundary);
//     EXPECT_EQ(diag.getVData(2).value().type, VertexType::Z);
//     EXPECT_EQ(diag.getVData(3).value().type, VertexType::Z);
//     EXPECT_EQ(diag.getVData(4).value().type, VertexType::X);
//     EXPECT_EQ(diag.getVData(7).value().type, VertexType::Z);
//     EXPECT_EQ(diag.getVData(8).value().type, VertexType::Z);
//     EXPECT_EQ(diag.getVData(9).value().type, VertexType::X);
//     EXPECT_EQ(diag.getVData(10).value().type, VertexType::Boundary);
//     EXPECT_EQ(diag.getVData(11).value().type, VertexType::Boundary);

//     EXPECT_TRUE(diag.isDeleted(1));
//     EXPECT_TRUE(diag.isDeleted(3));
// }

TEST_F(ZXDiagramTest, adjoint) {
    diag = diag.adjoint();

    EXPECT_EQ(diag.getEdge(0, 4).value().type, EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(5, 6).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(2, 6).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(3, 6).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(4, 5).value().type, EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 1).value().type, EdgeType::Simple);

    EXPECT_EQ(diag.getVData(0).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(1).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(2).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(3).value().type, VertexType::Boundary);
    EXPECT_EQ(diag.getVData(4).value().type, VertexType::Z);
    EXPECT_EQ(diag.getVData(5).value().type, VertexType::Z);
    EXPECT_EQ(diag.getVData(6).value().type, VertexType::X);
}

TEST_F(ZXDiagramTest, approximate) {
    ZXDiagram almostId(3);

    almostId.removeEdge(0, 3);
    const auto v = almostId.addVertex(0, 1, PiExpression(PiRational(1e-8)));
    almostId.addEdge(0, v);
    almostId.addEdge(v, 3);

    EXPECT_FALSE(almostId.phase(v).isZero());

    almostId.approximateCliffords(1e-7);

    EXPECT_TRUE(almostId.phase(v).isZero());
}

TEST_F(ZXDiagramTest, ancilla) {
    ZXDiagram cx(2);
    cx.removeEdge(0, 2);
    cx.removeEdge(1, 3);
    const auto tar = cx.addVertex(0, 0, PiExpression{}, VertexType::X);

    const auto ctrl = cx.addVertex(1);

    cx.addEdge(tar, ctrl);
    cx.addEdge(0, tar);
    cx.addEdge(tar, 2);
    cx.addEdge(1, ctrl);
    cx.addEdge(ctrl, 3);
    EXPECT_EQ(cx.getInputs().size(), 2);
    EXPECT_EQ(cx.getOutputs().size(), 2);
    EXPECT_EQ(cx.getNVertices(), 6);

    cx.makeAncilla(1);
    EXPECT_EQ(cx.getInputs().size(), 1);
    EXPECT_EQ(cx.getOutputs().size(), 1);
    EXPECT_EQ(cx.getNVertices(), 6);

    fullReduce(cx);

    EXPECT_EQ(cx.getNEdges(), 1);
    EXPECT_EQ(cx.getNVertices(), 2);
    EXPECT_TRUE(cx.isIdentity());
}

TEST_F(ZXDiagramTest, RemoveScalarSubDiagram) {
    ZXDiagram idWithScal(1);

    const auto v = idWithScal.addVertex(1);
    const auto w = idWithScal.addVertex(2);
    idWithScal.addEdge(v, w);

    fullReduce(idWithScal);

    EXPECT_EQ(idWithScal.getNVertices(), 2);
    EXPECT_EQ(idWithScal.getNEdges(), 1);
    EXPECT_TRUE(idWithScal.isDeleted(v));
    EXPECT_TRUE(idWithScal.isDeleted(w));
}

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

TEST_F(ZXDiagramTest, testOutputs) {
    ZXDiagram   linear    = makeIdentityDiagram(1, 3);
    const auto& inSpiders = linear.getInputSpiders();
    EXPECT_EQ(inSpiders.size(), 1);
    EXPECT_EQ(inSpiders[0], 2);

    const auto& outSpiders = linear.getOutputSpiders();
    EXPECT_EQ(outSpiders.size(), 1);
    EXPECT_EQ(outSpiders[0], 4);
}

// TEST_F(ZXDiagramTest, testGFlow) {
//     ZXDiagram   linear   = makeIdentityDiagram(1, 3);
//     const auto& gFlowOpt = linear.computeGFlow();
//     EXPECT_TRUE(gFlowOpt.has_value());
//     const auto& [ord, g] = gFlowOpt.value();

//     ASSERT_THAT(g[2], ElementsAre(3));
//     ASSERT_THAT(g[3], ElementsAre(4));

//     ASSERT_THAT(ord, ElementsAre(_, _, 2, 1, 0));
// }

static ZXDiagram makeEmptyDiagram(const std::size_t nqubits) {
    auto diag = ::makeIdentityDiagram(nqubits, 0);
    for (std::size_t i = 0; i < nqubits; ++i) {
        diag.removeEdge(i, i + nqubits);
    }
    return diag;
}

// TEST_F(ZXDiagramTest, testGFlow2) {
//     ZXDiagram diag = makeEmptyDiagram(3);
//     diag.addVertex(0, 0); // 6
//     diag.addVertex(1, 0); // 7
//     diag.addVertex(2, 0); // 8
//     diag.addVertex(0, 0); // 9
//     diag.addVertex(1, 0); // 10
//     diag.addVertex(2, 0); // 11

//     diag.addEdge(0, 6, EdgeType::Hadamard);
//     diag.addEdge(1, 7, EdgeType::Hadamard);
//     diag.addEdge(2, 8, EdgeType::Hadamard);

//     diag.addEdge(9, 3, EdgeType::Hadamard);
//     diag.addEdge(10, 4, EdgeType::Hadamard);
//     diag.addEdge(11, 5, EdgeType::Hadamard);

//     diag.addEdge(6, 9, EdgeType::Hadamard);
//     diag.addEdge(6, 11, EdgeType::Hadamard);

//     diag.addEdge(7, 9, EdgeType::Hadamard);
//     diag.addEdge(7, 10, EdgeType::Hadamard);
//     diag.addEdge(7, 11, EdgeType::Hadamard);

//     diag.addEdge(8, 10, EdgeType::Hadamard);
//     diag.addEdge(8, 11, EdgeType::Hadamard);

//     const auto& gFlowOpt = diag.computeGFlow();

//     EXPECT_TRUE(gFlowOpt.has_value());
//     const auto& [ord, g] = gFlowOpt.value();
//     ASSERT_THAT(g[6], ElementsAre(10, 11));
//     ASSERT_THAT(g[7], ElementsAre(9, 10, 11));
//     ASSERT_THAT(g[8], ElementsAre(9, 11));

//     ASSERT_THAT(ord, ElementsAre(_, _, _, _, _, _, 1, 1, 1, 0, 0, 0));
// }

// TEST_F(ZXDiagramTest, noGFlow) {
//     ZXDiagram diag = makeIdentityDiagram(2, 2);
//     diag.addEdge(4, 6, EdgeType::Hadamard);
//     diag.addEdge(4, 7, EdgeType::Hadamard);
//     diag.addEdge(5, 6, EdgeType::Hadamard);
//     diag.addEdge(5, 7, EdgeType::Hadamard);

//     const auto& gFlowOpt = diag.computeGFlow();

//     EXPECT_FALSE(gFlowOpt.has_value());
// }

// TEST_F(ZXDiagramTest, gFlowBig) {
//     const auto& gFlowOpt = diag.computeGFlow();
//     diag.toGraphlike();

//     EXPECT_TRUE(gFlowOpt.has_value());
// }

TEST_F(ZXDiagramTest, MeasurementTest) {
    auto pattern = makeEmptyDiagram(1);

    const auto& PI_2 = PiExpression(PiRational(1, 2));
    const auto& PI_4 = PiExpression(PiRational(1, 4));

    //add XY Measurement
    pattern.addVertex(0, 0); // 2
    pattern.addEdge(0, 2);

    //another XY Measurement
    pattern.addVertex(0, 0);        // 3
    pattern.addVertex(0, -1);       // 4
    pattern.addVertex(0, -2, PI_4); // 5
    pattern.addHadamardEdge(2, 3);
    pattern.addHadamardEdge(3, 4);
    pattern.addHadamardEdge(4, 5);

    // YZ Measurement
    pattern.addVertex(0, 0);        // 6
    pattern.addVertex(0, -1, PI_4); // 7
    pattern.addHadamardEdge(3, 6);
    pattern.addHadamardEdge(6, 7);

    // XZ Measurement
    pattern.addVertex(0, 0, PI_2);  // 8
    pattern.addVertex(0, -1, PI_4); // 9
    pattern.addHadamardEdge(6, 8);
    pattern.addHadamardEdge(8, 9);

    pattern.addEdge(8, 1);

    auto measurementOpt = pattern.getMeasurementPlane(2);
    EXPECT_TRUE(measurementOpt.has_value());
    auto meas = measurementOpt.value();
    EXPECT_EQ(meas.plane, MeasurementType::XY);
    EXPECT_EQ(meas.angle, PiRational{});

    measurementOpt = pattern.getMeasurementPlane(3);
    EXPECT_TRUE(measurementOpt.has_value());
    meas = measurementOpt.value();
    EXPECT_EQ(meas.plane, MeasurementType::XY);
    EXPECT_EQ(meas.angle, PiRational(1, 4));

    measurementOpt = pattern.getMeasurementPlane(6);
    EXPECT_TRUE(measurementOpt.has_value());
    meas = measurementOpt.value();
    EXPECT_EQ(meas.plane, MeasurementType::YZ);
    EXPECT_EQ(meas.angle, PiRational(1, 4));

    measurementOpt = pattern.getMeasurementPlane(8);
    EXPECT_TRUE(measurementOpt.has_value());
    meas = measurementOpt.value();
    EXPECT_EQ(meas.plane, MeasurementType::XZ);
    EXPECT_EQ(meas.angle, PiRational(1, 4));

    // check that gadgets are not recognised as qubits
    EXPECT_FALSE(pattern.getMeasurementPlane(4).has_value());
    EXPECT_FALSE(pattern.getMeasurementPlane(5).has_value());
    EXPECT_FALSE(pattern.getMeasurementPlane(7).has_value());
    EXPECT_FALSE(pattern.getMeasurementPlane(9).has_value());
}
