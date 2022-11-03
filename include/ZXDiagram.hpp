#pragma once

#include "Definitions.hpp"
#include "Expression.hpp"
#include "Utils.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace zx {
    class ZXDiagram {
    public:
        ZXDiagram() = default;
        // create n_qubit identity_diagram
        explicit ZXDiagram(std::size_t nqubits);

        void addEdge(Vertex from, Vertex to, EdgeType type = EdgeType::Simple);
        void addHadamardEdge(const Vertex from, const Vertex to) {
            addEdge(from, to, EdgeType::Hadamard);
        };
        void addEdgeParallelAware(Vertex from, Vertex to,
                                  EdgeType eType = EdgeType::Simple);
        void removeEdge(Vertex from, Vertex to);

        Vertex addVertex(const VertexData& data);
        Vertex addVertex(Qubit qubit, Col col = 0,
                         const PiExpression& phase = PiExpression(),
                         VertexType          type  = VertexType::Z);
        void   addQubit();
        void   addQubits(zx::Qubit n);
        void   removeVertex(Vertex toRemove);

        [[nodiscard]] std::size_t getNdeleted() const { return deleted.size(); }
        [[nodiscard]] std::size_t getNVertices() const { return nvertices; }
        [[nodiscard]] std::size_t getNEdges() const { return nedges; }
        [[nodiscard]] std::size_t getNQubits() const { return inputs.size(); }

        [[nodiscard]] bool                     connected(Vertex from, Vertex to) const;
        [[nodiscard]] std::optional<Edge>      getEdge(Vertex from, Vertex to) const;
        [[nodiscard]] const std::vector<Edge>& incidentEdges(const Vertex v) const { return edges[v]; }
        [[nodiscard]] const Edge&              incidentEdge(const Vertex v, const std::size_t n) { return edges[v][n]; }

        [[nodiscard]] std::size_t degree(const Vertex v) const { return edges[v].size(); }

        [[nodiscard]] const PiExpression& phase(const Vertex v) const {
            return vertices[v].value().phase;
        }

        [[nodiscard]] Qubit qubit(const Vertex v) const {
            return vertices[v].value().qubit;
        }

        [[nodiscard]] VertexType type(const Vertex v) const {
            return vertices[v].value().type;
        }

        [[nodiscard]] std::optional<VertexData> getVData(const Vertex v) const {
            return vertices[v];
        }

        [[nodiscard]] std::vector<std::pair<Vertex, VertexData&>> getVertices();
        [[nodiscard]] std::vector<std::pair<Vertex, Vertex>>      getEdges();

        [[nodiscard]] const std::vector<Vertex>& getInputs() const {
            return inputs;
        }
        [[nodiscard]] Vertex getInput(const std::size_t i) const { return inputs[i]; }

        [[nodiscard]] const std::vector<Vertex>& getOutputs() const {
            return outputs;
        }

        [[nodiscard]] Vertex getOutput(const std::size_t i) const { return outputs[i]; }

        [[nodiscard]] bool isDeleted(const Vertex v) const {
            return !vertices[v].has_value();
        }

        [[nodiscard]] bool isBoundaryVertex(const Vertex v) const {
            return vertices[v].value().type == VertexType::Boundary;
        }

        [[nodiscard]] bool isInput(Vertex v) const;
        [[nodiscard]] bool isOutput(Vertex v) const;

        void addPhase(const Vertex v, const PiExpression& phase) {
            vertices[v].value().phase += phase;
        }

        void setPhase(const Vertex v, const PiExpression& phase) {
            vertices[v].value().phase = phase;
        }

        void setType(const Vertex v, const VertexType type) {
            vertices[v].value().type = type;
        }

        void toGraphlike();

        [[nodiscard]] bool isIdentity() const;

        [[nodiscard]] ZXDiagram adjoint() const;

        ZXDiagram& invert();

        ZXDiagram& concat(const ZXDiagram& rhs);
        ZXDiagram& operator+=(const ZXDiagram& rhs) { return this->concat(rhs); }

        void makeAncilla(Qubit qubit);
        void makeAncilla(Qubit in, Qubit out);

        void approximateCliffords(fp tolerance);

        void removeDisconnectedSpiders();

        void                       addGlobalPhase(const PiExpression& phase);
        [[nodiscard]] PiExpression getGlobalPhase() const {
            return globalPhase;
        }
        [[nodiscard]] bool globalPhaseIsZero() const {
            return globalPhase.isZero();
        }

        //implements the polynomial algorithm for finding a maximally delayed g-flow from https://arxiv.org/pdf/0709.2670.pdf
        // [[nodiscard]] std::optional<flow> computeGFlow();

        [[nodiscard]] std::vector<Vertex> getOutputSpiders() const;
        [[nodiscard]] std::vector<Vertex> getInputSpiders() const;

        // Results not guaranteed to make sense if diagram not graph-like and identity
        [[nodiscard]] std::optional<Measurement> getMeasurementPlane(Vertex v) const;

    private:
        std::vector<std::vector<Edge>>         edges;
        std::vector<std::optional<VertexData>> vertices;
        std::vector<Vertex>                    deleted;
        std::vector<Vertex>                    inputs;
        std::vector<Vertex>                    outputs;
        std::size_t                            nvertices   = 0;
        std::size_t                            nedges      = 0;
        PiExpression                           globalPhase = {};

        std::vector<Vertex> initGraph(std::size_t nqubits);
        void                closeGraph(const std::vector<Vertex>& qubitVertices);

        void removeHalfEdge(Vertex from, Vertex to);

        std::vector<Edge>::iterator getEdgePtr(Vertex from, Vertex to);

        // void gFlowAux(const gf2Mat& adjMat, std::vector<Vertex> in) const;
        gf2Mat                                              getAdjMat();
        std::pair<std::vector<Vertex>, std::vector<Vertex>> getNonProcessed(const std::vector<Vertex>& out) const;

        static bool isIn(const Vertex& v, const std::vector<Vertex>& vertices);

        gf2Mat constructLinearSystem(const gf2Mat& adjMat, std::vector<Vertex> out, std::vector<Vertex> out_prime, std::vector<Vertex> us) const;

        std::vector<Vertex> solutionFromTriangular(const gf2Mat& triu, const std::vector<Vertex>& us, std::size_t offset, const std::vector<Vertex>& out_prime, std::vector<std::vector<Vertex>>& g) const;
        std::vector<Vertex> getConnectedSet(const std::vector<Vertex> s) const;

        static std::vector<Vertex> fromIdxVec(const std::vector<bool>& indicator, const std::vector<Vertex>& set);
    };
} // namespace zx
