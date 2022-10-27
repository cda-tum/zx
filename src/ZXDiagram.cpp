#include "ZXDiagram.hpp"

#include "Definitions.hpp"
#include "Expression.hpp"
#include "Rational.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace zx {

    ZXDiagram::ZXDiagram(const std::size_t nqubits) {
        auto qubitVertices = initGraph(nqubits);
        closeGraph(qubitVertices);
    }

    void ZXDiagram::addEdge(const Vertex from, const Vertex to, const EdgeType type) {
        edges[from].emplace_back(to, type);
        edges[to].emplace_back(from, type);
        ++nedges;
    }

    void ZXDiagram::addEdgeParallelAware(const Vertex from, const Vertex to,
                                         const EdgeType eType) { // TODO: Scalars
        if (from == to) {
            if (type(from) != VertexType::Boundary && eType == EdgeType::Hadamard) {
                addPhase(from, PiExpression(PiRational(1, 1)));
            }
            return;
        }

        const auto edgeIt = getEdgePtr(from, to);

        if (edgeIt == edges[from].end()) {
            addEdge(from, to, eType);
            return;
        }

        if (type(from) == VertexType::Boundary || type(to) == VertexType::Boundary) {
            return;
        }

        if (type(from) == type(to)) {
            if (edgeIt->type == EdgeType::Hadamard && eType == EdgeType::Hadamard) {
                edges[from].erase(edgeIt);
                removeHalfEdge(to, from);
                --nedges;
            } else if (edgeIt->type == EdgeType::Hadamard &&
                       eType == EdgeType::Simple) {
                edgeIt->type = EdgeType::Simple;
                getEdgePtr(to, from)->toggle();
                addPhase(from, PiExpression(PiRational(1, 1)));
            } else if (edgeIt->type == EdgeType::Simple &&
                       eType == EdgeType::Hadamard) {
                addPhase(from, PiExpression(PiRational(1, 1)));
            }
        } else {
            if (edgeIt->type == EdgeType::Simple && eType == EdgeType::Simple) {
                edges[from].erase(edgeIt);
                removeHalfEdge(to, from);
                --nedges;
            } else if (edgeIt->type == EdgeType::Hadamard &&
                       eType == EdgeType::Simple) {
                addPhase(from, PiExpression(PiRational(1, 1)));
            } else if (edgeIt->type == EdgeType::Simple &&
                       eType == EdgeType::Hadamard) {
                edgeIt->type = EdgeType::Hadamard;
                getEdgePtr(to, from)->toggle();
                addPhase(from, PiExpression(PiRational(1, 1)));
            }
        }
    }

    void ZXDiagram::removeEdge(const Vertex from, const Vertex to) {
        removeHalfEdge(from, to);
        removeHalfEdge(to, from);
        --nedges;
    }

    void ZXDiagram::removeHalfEdge(const Vertex from, const Vertex to) {
        auto& incident = edges[from];
        incident.erase(std::remove_if(incident.begin(), incident.end(),
                                      [&](auto& edge) { return edge.to == to; }),
                       incident.end());
    }

    Vertex ZXDiagram::addVertex(const VertexData& data) {
        ++nvertices;
        Vertex v = 0;
        if (!deleted.empty()) {
            v = deleted.back();
            deleted.pop_back();
            vertices[v] = data;
            edges[v].clear();
            return v;
        }

        v = nvertices;
        vertices.emplace_back(data);
        edges.emplace_back();

        return nvertices - 1;
    }

    Vertex ZXDiagram::addVertex(const Qubit qubit, const Col col, const PiExpression& phase,
                                const VertexType type) {
        return addVertex({col, qubit, phase, type});
    }

    void ZXDiagram::addQubit() {
        auto in  = addVertex(static_cast<zx::Qubit>(getNQubits()) + 1, 0, PiExpression(), VertexType::Boundary);
        auto out = addVertex(static_cast<zx::Qubit>(getNQubits()) + 1, 0, PiExpression(), VertexType::Boundary);
        inputs.emplace_back(in);
        outputs.emplace_back(out);
    }
    void ZXDiagram::addQubits(const Qubit n) {
        for (zx::Qubit i = 0; i < n; ++i) {
            addQubit();
        }
    }

    void ZXDiagram::removeVertex(const Vertex toRemove) {
        deleted.push_back(toRemove);
        vertices[toRemove].reset();
        --nvertices;

        for (const auto& [to, _]: incidentEdges(toRemove)) {
            removeHalfEdge(to, toRemove);
            --nedges;
        }
    }

    [[nodiscard]] bool ZXDiagram::connected(const Vertex from, const Vertex to) const {
        if (isDeleted(from) || isDeleted(to)) {
            return false;
        }

        const auto& incident = edges[from];
        const auto  edge     = std::find_if(incident.begin(), incident.end(),
                                            [&](const auto& e) { return e.to == to; });
        return edge != incident.end();
    }

    [[nodiscard]] std::optional<Edge> ZXDiagram::getEdge(const Vertex from,
                                                         const Vertex to) const {
        std::optional<Edge> ret;
        const auto&         incident = edges[from];
        const auto          edge     = std::find_if(incident.begin(), incident.end(),
                                                    [&](const auto& e) { return e.to == to; });
        if (edge != incident.end()) {
            ret = *edge;
        }
        return ret;
    }

    std::vector<Edge>::iterator ZXDiagram::getEdgePtr(const Vertex from, const Vertex to) {
        auto& incident = edges[from];
        auto  edge     = std::find_if(incident.begin(), incident.end(),
                                      [&](const auto& e) { return e.to == to; });
        return edge;
    }

    [[nodiscard]] std::vector<std::pair<Vertex, VertexData&>>
    ZXDiagram::getVertices() {
        Vertices verts(vertices);
        return {verts.begin(), verts.end()};
    }

    [[nodiscard]] std::vector<std::pair<Vertex, Vertex>> ZXDiagram::getEdges() {
        Edges es(edges, vertices);
        return {es.begin(), es.end()};
    }

    bool ZXDiagram::isInput(const Vertex v) const {
        return std::find(inputs.begin(), inputs.end(), v) != inputs.end();
    }
    bool ZXDiagram::isOutput(const Vertex v) const {
        return std::find(outputs.begin(), outputs.end(), v) != outputs.end();
    }

    void ZXDiagram::toGraphlike() {
        const auto nverts = vertices.size();
        for (Vertex v = 0U; v < nverts; ++v) {
            if (!vertices[v].has_value()) {
                continue;
            }
            if (vertices[v].value().type == VertexType::X) {
                for (auto& edge: edges[v]) {
                    edge.toggle();
                    // toggle corresponding edge in other direction
                    getEdgePtr(edge.to, v)->toggle();
                }

                vertices[v].value().type = VertexType::Z;
            }
        }
    }

    [[nodiscard]] ZXDiagram ZXDiagram::adjoint() const {
        ZXDiagram copy = *this;
        copy.invert();
        return copy;
    }

    ZXDiagram& ZXDiagram::invert() {
        const auto h = inputs;
        inputs       = outputs;
        outputs      = h;

        for (auto& data: vertices) {
            if (data.has_value()) {
                data.value().phase = -data.value().phase;
            }
        }
        return *this;
    }

    ZXDiagram& ZXDiagram::concat(const ZXDiagram& rhs) {
        if (rhs.getNQubits() != this->getNQubits()) {
            throw ZXException(
                    "Cannot concatenate Diagrams with differing number of qubits!");
        }

        std::unordered_map<Vertex, Vertex> newVs;
        const auto                         nverts = rhs.vertices.size();
        for (std::size_t i = 0; i < nverts; ++i) {
            if (!rhs.vertices[i].has_value() || rhs.isInput(i)) {
                continue;
            }

            const auto newV = addVertex(rhs.vertices[i].value());
            newVs[i]        = newV;
        }

        for (std::size_t i = 0; i < nverts; ++i) { // add new edges
            if (!rhs.vertices[i].has_value() || rhs.isInput(i)) {
                continue;
            }

            for (const auto& [to, type]: rhs.edges[i]) {
                if (!rhs.isInput(to)) {
                    if (i < to) { // make sure not to add edge twice
                        addEdge(newVs[i], newVs[to], type);
                    }
                } else {
                    const auto outV = outputs[rhs.qubit(to)];
                    for (const auto& [interior_v, interior_type]:
                         edges[outV]) { // redirect edges going to outputs
                        if (interior_type == type) {
                            addEdge(interior_v, newVs[i], EdgeType::Simple);
                        } else {
                            addEdge(interior_v, newVs[i], EdgeType::Hadamard);
                        }
                    }
                }
            }
        } // add new edges

        const auto nOuptuts = outputs.size();
        for (size_t i = 0; i < nOuptuts; ++i) {
            removeVertex(outputs[i]);
            outputs[i] = newVs[rhs.outputs[i]];
        }

        this->addGlobalPhase(-rhs.globalPhase);
        return *this;
    }

    bool ZXDiagram::isIdentity() const {
        if (nedges != inputs.size() || !globalPhase.isZero()) {
            return false;
        }

        const auto nInputs = inputs.size();
        for (size_t i = 0; i < nInputs; ++i) {
            if (!connected(inputs[i], outputs[i])) {
                return false;
            }
        }
        return true;
    }

    std::vector<Vertex> ZXDiagram::initGraph(const std::size_t nqubits) {
        std::vector<Vertex> qubitVertices(nqubits, 0);

        const auto nVerts = qubitVertices.size();
        for (size_t i = 0; i < nVerts; ++i) {
            const auto v = addVertex(
                    {1, static_cast<Qubit>(i), PiExpression(), VertexType::Boundary});
            qubitVertices[i] = v;
            inputs.push_back(v);
        }

        return qubitVertices;
    }

    void ZXDiagram::closeGraph(const std::vector<Vertex>& qubitVertices) {
        for (const Vertex v: qubitVertices) {
            const VertexData vData = vertices[v].value();
            const Vertex     newV  = addVertex({vData.col + 1,
                                                vData.qubit,
                                                PiExpression(),
                                                VertexType::Boundary});
            addEdge(v, newV);
            outputs.push_back(newV);
        }
    }

    void ZXDiagram::makeAncilla(const Qubit qubit) {
        makeAncilla(qubit, qubit);
    }

    void ZXDiagram::makeAncilla(const Qubit in, const Qubit out) {
        const auto inV  = inputs[in];
        const auto outV = outputs[out];
        inputs.erase(inputs.begin() + in);
        outputs.erase(outputs.begin() + out);

        setType(inV, VertexType::X);
        setType(outV, VertexType::X);
    }

    void ZXDiagram::approximateCliffords(const fp tolerance) {
        for (auto& v: vertices) {
            if (v.has_value()) {
                roundToClifford(v.value().phase, tolerance);
            }
        }
    }

    void ZXDiagram::removeDisconnectedSpiders() {
        auto connectedToBoundary = [&](const Vertex v) {
            std::unordered_set<Vertex> visited{};
            std::vector<Vertex>        stack{};
            stack.push_back(v);

            while (!stack.empty()) {
                auto w = stack.back();
                stack.pop_back();

                if (visited.find(w) != visited.end()) {
                    continue;
                }

                visited.emplace(w);

                if (isInput(w) || isOutput(w)) {
                    return true;
                }

                for (const auto [to, _]: incidentEdges(w)) {
                    stack.push_back(to);
                }
            }
            return false;
        };

        const auto nVerts = vertices.size();
        for (Vertex v = 0; v < nVerts; ++v) {
            if (!isDeleted(v) && !connectedToBoundary(v)) {
                removeVertex(v);
            }
        }
    }

    void ZXDiagram::addGlobalPhase(const PiExpression& phase) {
        globalPhase += phase;
    }

    [[nodiscard]] std::optional<std::vector<std::uint32_t>> ZXDiagram::gFlow() {
        // TODO : check if diagram is graph-like
        std::vector<uint32_t>            partOrd(nvertices, 0U); //TODO: technically only the outputs are init to 0
        std::vector<std::vector<Vertex>> g(nvertices);
        const auto&                      adjMat = getAdjMat();
        std::vector<Vertex>              out_prime;
        std::vector<Vertex>              c;
        std::vector<Vertex>              out = outputs;

        do {
            for (const auto& v: outputs) {
                if (std::find(inputs.begin(), inputs.end(), v) == inputs.end())
                    out_prime.emplace_back(v);
            }
            std::vector<Vertex> us = getNonOutputs(out);

            auto system = constructLinearSystem(adjMat, out, out_prime, us);

            auto systemFlint = getFlintMatrix(system);
            systemFlint.set_rref(); // bring to upper triangular form
            c = solutionFromTriangular(getMatrixFromFlint(systemFlint), us, out_prime.size(), g);
        } while (!c.empty());

        return partOrd;
    }

    gf2Mat ZXDiagram::getAdjMat() {
        gf2Mat adjMat{nvertices, gf2Vec(nvertices, false)};
        for (const auto& [from, to]: getEdges()) {
            adjMat[from][to] = true;
            adjMat[to][from] = true;
        }
        return adjMat;
    }

    std::vector<Vertex> ZXDiagram::getNonOutputs(const std::vector<Vertex>& out) const {
        std::vector<Vertex> nonOuts;
        for (Vertex i = 0; i < vertices.size(); ++i) {
            if (vertices[i].has_value() && isIn(i, out))
                nonOuts.emplace_back(i);
        }
        return nonOuts;
    }

    bool ZXDiagram::isIn(const Vertex& v, const std::vector<Vertex>& vertices) {
        return std::find(vertices.begin(), vertices.end(), v) != vertices.end();
    }

    gf2Mat ZXDiagram::constructLinearSystem(const gf2Mat& adjMat, std::vector<Vertex> out, std::vector<Vertex> out_prime, std::vector<Vertex> us) const {
        gf2Mat system(nvertices - out.size(), gf2Vec(nvertices, false));

        for (std::size_t row = 0; row < us.size(); ++row) {
            for (std::size_t col = 0; col < out_prime.size(); ++col) {
                system[row][row] = adjMat[us[row]][out_prime[col]];
            }
        }

        for (std::size_t curr_col = out_prime.size(); curr_col < nvertices; ++curr_col) {
            auto idx              = us[curr_col - out_prime.size()];
            system[idx][curr_col] = true;
        }

        return system;
    }

    std::vector<Vertex> ZXDiagram::solutionFromTriangular(const gf2Mat& triu, const std::vector<Vertex>& us, std::size_t offset, std::vector<std::vector<Vertex>>& g) const {
        for (std::size_t col = offset; col < nvertices; ++col) {
            gf2Vec      sol(triu.size(), false);
            std::size_t maxNonZeroRow = 0;
            for (std::size_t i = 0; i < offset; ++i) {
                if (triu[offset - i][offset - i]) {
                    maxNonZeroRow = offset - i;
                    break;
                }
            }

            //backpropagation
            for (std::size_t col = offset; col < triu[0].size(); ++col) {
                bool hasSol = true;
                for (std::size_t row = triu.size(); row > maxNonZeroRow; --row) {
                    if (triu[row][col]) {
                        hasSol = false;
                        break;
                    }
                }
                if (!hasSol)
                    continue;

                // TODO: replace by bitwise ops
                for (std::size_t row_comp = 0; row_comp < offset; ++row_comp) {
                    auto row     = offset - row_comp - 1; //TODO
                    int  par_lhs = std::accumulate(trio[row][row].begin() +)
                }
            }
        }
    }
} // namespace zx
