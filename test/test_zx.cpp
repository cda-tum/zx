#include "Definitions.hpp"
#include "Rational.hpp"
#include "Simplify.hpp"
#include "ZXDiagram.hpp"

#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>

class ZXDiagramTest: public ::testing::Test {
public:
    zx::ZXDiagram diag;

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
        diag = zx::ZXDiagram();
        diag.addQubits(2);
        diag.addVertex(0, 0, zx::PiExpression(), zx::VertexType::Z);
        diag.addEdge(0, 4, zx::EdgeType::Hadamard);
        diag.addVertex(0, 0, zx::PiExpression(), zx::VertexType::Z);
        diag.addEdge(4, 5);
        diag.addVertex(0, 0, zx::PiExpression(), zx::VertexType::X);
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

    EXPECT_EQ(diag.getEdge(0, 4).value().type, zx::EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(4, 5).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(2, 6).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 6).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 1).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(6, 3).value().type, zx::EdgeType::Simple);

    EXPECT_EQ(diag.getVData(0).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(1).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(2).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(3).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(4).value().type, zx::VertexType::Z);
    EXPECT_EQ(diag.getVData(5).value().type, zx::VertexType::Z);
    EXPECT_EQ(diag.getVData(6).value().type, zx::VertexType::X);

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

    EXPECT_EQ(diag.getEdge(0, 4).value().type, zx::EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(5, 6).value().type, zx::EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(2, 6).value().type, zx::EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(3, 6).value().type, zx::EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(4, 5).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 1).value().type, zx::EdgeType::Simple);

    EXPECT_EQ(diag.getVData(0).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(1).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(2).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(3).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(4).value().type, zx::VertexType::Z);
    EXPECT_EQ(diag.getVData(5).value().type, zx::VertexType::Z);
    EXPECT_EQ(diag.getVData(6).value().type, zx::VertexType::Z);

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

//     EXPECT_EQ(diag.getEdge(0, 4).value().type, zx::EdgeType::Hadamard);
//     EXPECT_EQ(diag.getEdge(5, 6).value().type, zx::EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(2, 6).value().type, zx::EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(3, 6).value().type, zx::EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(5, 9).value().type, zx::EdgeType::Hadamard);
//     EXPECT_EQ(diag.getEdge(4, 9).value().type, zx::EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(7, 8).value().type, zx::EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(8, 10).value().type, zx::EdgeType::Simple);
//     EXPECT_EQ(diag.getEdge(9, 11).value().type, zx::EdgeType::Simple);

//     EXPECT_EQ(diag.getVData(0).value().type, zx::VertexType::Boundary);
//     EXPECT_EQ(diag.getVData(1).value().type, zx::VertexType::Boundary);
//     EXPECT_EQ(diag.getVData(2).value().type, zx::VertexType::Z);
//     EXPECT_EQ(diag.getVData(3).value().type, zx::VertexType::Z);
//     EXPECT_EQ(diag.getVData(4).value().type, zx::VertexType::X);
//     EXPECT_EQ(diag.getVData(7).value().type, zx::VertexType::Z);
//     EXPECT_EQ(diag.getVData(8).value().type, zx::VertexType::Z);
//     EXPECT_EQ(diag.getVData(9).value().type, zx::VertexType::X);
//     EXPECT_EQ(diag.getVData(10).value().type, zx::VertexType::Boundary);
//     EXPECT_EQ(diag.getVData(11).value().type, zx::VertexType::Boundary);

//     EXPECT_TRUE(diag.isDeleted(1));
//     EXPECT_TRUE(diag.isDeleted(3));
// }

TEST_F(ZXDiagramTest, adjoint) {
    diag = diag.adjoint();

    EXPECT_EQ(diag.getEdge(0, 4).value().type, zx::EdgeType::Hadamard);
    EXPECT_EQ(diag.getEdge(5, 6).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(2, 6).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(3, 6).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(4, 5).value().type, zx::EdgeType::Simple);
    EXPECT_EQ(diag.getEdge(5, 1).value().type, zx::EdgeType::Simple);

    EXPECT_EQ(diag.getVData(0).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(1).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(2).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(3).value().type, zx::VertexType::Boundary);
    EXPECT_EQ(diag.getVData(4).value().type, zx::VertexType::Z);
    EXPECT_EQ(diag.getVData(5).value().type, zx::VertexType::Z);
    EXPECT_EQ(diag.getVData(6).value().type, zx::VertexType::X);
}

TEST_F(ZXDiagramTest, approximate) {
    zx::ZXDiagram almostId(3);

    almostId.removeEdge(0, 3);
    const auto v = almostId.addVertex(0, 1, zx::PiExpression(zx::PiRational(1e-8)));
    almostId.addEdge(0, v);
    almostId.addEdge(v, 3);

    EXPECT_FALSE(almostId.phase(v).isZero());

    almostId.approximateCliffords(1e-7);

    EXPECT_TRUE(almostId.phase(v).isZero());
}

TEST_F(ZXDiagramTest, ancilla) {
    zx::ZXDiagram cx(2);
    cx.removeEdge(0, 2);
    cx.removeEdge(1, 3);
    const auto tar = cx.addVertex(0, 0, zx::PiExpression{}, zx::VertexType::X);

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

    zx::fullReduce(cx);

    EXPECT_EQ(cx.getNEdges(), 1);
    EXPECT_EQ(cx.getNVertices(), 2);
    EXPECT_TRUE(cx.isIdentity());
}

TEST_F(ZXDiagramTest, RemoveScalarSubDiagram) {
    zx::ZXDiagram idWithScal(1);

    const auto v = idWithScal.addVertex(1);
    const auto w = idWithScal.addVertex(2);
    idWithScal.addEdge(v, w);

    zx::fullReduce(idWithScal);

    EXPECT_EQ(idWithScal.getNVertices(), 2);
    EXPECT_EQ(idWithScal.getNEdges(), 1);
    EXPECT_TRUE(idWithScal.isDeleted(v));
    EXPECT_TRUE(idWithScal.isDeleted(w));
}

static zx::ZXDiagram makeIdentityDiagram(const std::size_t nqubits,
                                         const std::size_t spidersPerQubit) {
    zx::ZXDiagram           diag(nqubits);
    std::vector<zx::Vertex> rightmostVertices = diag.getInputs();

    for (std::size_t i = 0; i < nqubits; ++i) {
        diag.removeEdge(i, i + nqubits);
    }

    // add identity spiders
    for (zx::Qubit qubit = 0; static_cast<std::size_t>(qubit) < nqubits; ++qubit) {
        for (std::size_t j = 0; j < spidersPerQubit; ++j) {
            const zx::Vertex v = diag.addVertex(qubit);
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
    zx::ZXDiagram linear    = makeIdentityDiagram(1, 3);
    const auto&   inSpiders = linear.getInputSpiders();
    EXPECT_EQ(inSpiders.size(), 1);
    EXPECT_EQ(inSpiders[0], 2);

    const auto& outSpiders = linear.getOutputSpiders();
    EXPECT_EQ(outSpiders.size(), 1);
    EXPECT_EQ(outSpiders[0], 4);
}

TEST_F(ZXDiagramTest, testGFlow) {
    zx::ZXDiagram linear   = makeIdentityDiagram(1, 3);
    const auto&   gFlowOpt = linear.computeGFlow();
    if (gFlowOpt.has_value()) {
        const auto& [ord, g] = gFlowOpt.value();

        for (const auto& [v, _]: linear.getVertices()) {
            std::cout << "correcting " << v << " on vertices: ";
            for (const auto& cv: g[v]) std::cout << cv << ", ";
            std::cout << std::endl;
        }
    }
}
