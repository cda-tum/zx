#include "MBQC.hpp"

#include "Definitions.hpp"

#include <cstddef>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace zx {
    MBQCPattern::MBQCPattern(const ZXDiagram& diag) {
        std::vector<Vertex>                mbqcToZxVertexMap{};
        std::unordered_map<Vertex, Vertex> zxToMbqcVertexMap{};
        std::unordered_set<Vertex>         nonQubitZXVertices{};

        for (std::size_t i = 0; i < diag.getNVertices(); ++i) {
            const auto& measOpt = diag.getMeasurementPlane(i);
            if (measOpt.has_value()) {
                mbqcToZxVertexMap.emplace_back(i);
                zxToMbqcVertexMap[i] = mbqcToZxVertexMap.size() - 1;
                measurements.emplace_back(measOpt.value());
            } else {
                nonQubitZXVertices.insert(i);
            }
        }

        edges.reserve(measurements.size());
        for (std::size_t i = 0; i < measurements.size(); ++i) {
            for (const auto& [to, _]: diag.incidentEdges(mbqcToZxVertexMap[i])) {
                if (nonQubitZXVertices.find(to) == nonQubitZXVertices.end()) {
                    auto patternTo = zxToMbqcVertexMap[to];
                    if (patternTo >= i) {
                        edges[i].emplace_back(patternTo);
                        edges[patternTo].emplace_back(i);
                        nedges++;
                    }
                }
            }
        }

        for (const auto& inBoundary: diag.getInputs()) {
            for (const auto& [to, _]: diag.incidentEdges(inBoundary)) {
                inputs.emplace_back(zxToMbqcVertexMap[to]);
            }
        }
        for (const auto& outBoundary: diag.getOutputs()) {
            for (const auto& [to, _]: diag.incidentEdges(outBoundary)) {
                outputs.emplace_back(zxToMbqcVertexMap[to]);
            }
        }
    }

    MBQCPattern::MBQCPattern(const std::initializer_list<Measurement>& measurements, const std::initializer_list<Vertex>& inputs, const std::initializer_list<Vertex>& outputs, const std::initializer_list<std::pair<Vertex, Vertex>>& edges):
        measurements(measurements), inputs(inputs), outputs(outputs), edges(measurements.size()) {
        for (const auto& [from, to]: edges) addEdge(from, to);
    }

    std::optional<flow> MBQCPattern::computeGFlow() const {
        // TODO : check if diagram is graph-like
        std::vector<uint32_t>            partOrd(getNVertices(), 0U); //TODO: technically only the outputs are init to 0
        std::vector<std::vector<Vertex>> g(getNVertices());
        const auto&                      adjMat = getAdjMat();
        std::vector<Vertex>              c;
        std::vector<Vertex>              out = getOutputs();
        std::vector<Vertex>              in  = getInputs();

        uint32_t k = 1;

        do {
            c = {};

            const auto& [us, relevantProcessed] = getNonProcessed(out);
            if (relevantProcessed.empty()) break;

            auto system = constructLinearSystem(adjMat, out, relevantProcessed, us);

            auto systemFlint = getFlintMatrix(system);
            systemFlint.set_rref(); // bring to upper triangular form

            c = solutionFromTriangular(getMatrixFromFlint(systemFlint), us, relevantProcessed.size(), relevantProcessed, g); //TODO may be more efficient when filtering out vertices not connected to any output

            for (const auto& v: c) {
                partOrd[v] = k;
            }
            k++;
            out.insert(out.end(), c.begin(), c.end());
        } while (!c.empty());

        if (out.size() != getNVertices())
            return {};

        std::vector<std::vector<Vertex>> depthToVertex(k);
        for (std::size_t v = 0; v < partOrd.size(); ++v) {
            depthToVertex[k - partOrd[v] - 1].emplace_back(v);
        }
        return flow{depthToVertex, g};
    }

    gf2Mat MBQCPattern::getAdjMat() const {
        gf2Mat adjMat{getNVertices(), gf2Vec(getNVertices(), false)};
        for (std::size_t from = 0; from < edges.size(); ++from) {
            for (const auto& to: edges[from]) {
                adjMat[from][to] = true;
            }
        }

        for (std::size_t i = 0; i < adjMat.size(); ++i) {
            adjMat[i][i] = true;
        }
        return adjMat;
    }

    std::pair<std::vector<Vertex>, std::vector<Vertex>> MBQCPattern::getNonProcessed(const std::vector<Vertex>& processed) const {
        std::vector<Vertex> nonProcessed;
        for (const Vertex& i: getConnectedSet(processed)) {
            if (!isIn(i, processed))
                nonProcessed.emplace_back(i);
        }

        std::vector<Vertex> relevantProcessed;
        for (const Vertex& i: getConnectedSet(nonProcessed)) {
            if (isIn(i, processed))
                relevantProcessed.emplace_back(i);
        }
        return std::make_pair(nonProcessed, relevantProcessed);
    }

    bool MBQCPattern::isIn(const Vertex& v, const std::vector<Vertex>& vertices) {
        return std::find(vertices.begin(), vertices.end(), v) != vertices.end();
    }

    gf2Mat MBQCPattern::constructLinearSystem(const gf2Mat& adjMat, std::vector<Vertex> out, std::vector<Vertex> relevantProcessed, std::vector<Vertex> toProcess) const {
        const auto& nRelevantVertices = relevantProcessed.size() + toProcess.size();
        gf2Mat      system(nRelevantVertices - relevantProcessed.size(), gf2Vec(nRelevantVertices, false));
        for (std::size_t row = 0; row < toProcess.size(); ++row) {
            for (std::size_t col = 0; col < relevantProcessed.size(); ++col) {
                system[row][col] = adjMat[toProcess[row]][relevantProcessed[col]];
            }
        }

        for (std::size_t i = 0; i < toProcess.size(); ++i) {
            const auto  u        = toProcess[i];
            const auto& curr_col = relevantProcessed.size() + i;
            if (measurements[u].plane == MeasurementType::XY) {
                system[curr_col - toProcess.size()][curr_col] = true;
            } else {
                for (std::size_t j = 0; j < toProcess.size(); ++j) {
                    system[j][curr_col] = isIn(toProcess[j], edges[u]);
                }
                if (measurements[u].plane == MeasurementType::XZ) {
                    system[curr_col - toProcess.size()][curr_col] = !system[curr_col - relevantProcessed.size()][curr_col];
                }
            }
        }
        return system;
    }

    std::vector<Vertex> MBQCPattern::solutionFromTriangular(const gf2Mat& triu, const std::vector<Vertex>& us, std::size_t offset, const std::vector<Vertex>& relevantProcessed, std::vector<std::vector<Vertex>>& g) const {
        std::vector<Vertex> c;

        const auto& nRelevantVertices = getNVertices();

        std::size_t maxNonZeroRow = 0;
        for (std::size_t i = 0; i < offset; ++i) {
            if (triu[triu.size() - i - 1][triu.size() - i - 1]) {
                maxNonZeroRow = triu.size() - i - 1;
                break;
            }
        }

        //backsubstitution
        for (std::size_t col = offset; col < triu[0].size(); ++col) {
            gf2Vec sol(offset, false);
            bool   hasSol = true;
            for (std::size_t row_comp = 0; row_comp < triu.size() - maxNonZeroRow - 1; ++row_comp) {
                auto row = triu.size() - row_comp - 1;
                if (triu[row][col]) {
                    hasSol = false;
                    break;
                }
            }

            if (!hasSol)
                continue;

            // TODO: replace by bitwise ops
            for (std::size_t row_comp = 0; row_comp <= maxNonZeroRow; ++row_comp) {
                auto row = maxNonZeroRow - row_comp;
                int  sum = static_cast<int>(triu[row][row]);
                for (std::size_t i = row + 1; i < sol.size(); ++i) {
                    sum += static_cast<int>(triu[row][i]) * sol[i];
                }

                if (sum == 0 && triu[row][col]) {
                    hasSol = false;
                    break;
                }
                sol[row] = sum & (static_cast<int>(triu[row][col]));
            }

            if (hasSol) {
                const auto u = us[col - offset];
                g[u]         = fromIdxVec(sol, relevantProcessed);
                if (measurements[u].plane != MeasurementType::XY) {
                    g[u].emplace_back(u);
                }
                c.emplace_back(us[col - offset]);
            }
        }
        return c;
    }

    std::vector<Vertex> MBQCPattern::fromIdxVec(const std::vector<bool>& indicator, const std::vector<Vertex>& set) {
        std::vector<Vertex> ret;
        for (std::size_t i = 0; i < indicator.size(); ++i) {
            if (indicator[i])
                ret.emplace_back(set[i]);
        }
        return ret;
    }

    std::vector<Vertex> MBQCPattern::getConnectedSet(const std::vector<Vertex> s) const {
        std::vector<Vertex> connected;
        for (const auto v: s) {
            for (const auto to: edges[v]) {
                const auto& p = std::lower_bound(connected.begin(), connected.end(), to);
                if (p == connected.end()) {
                    connected.emplace_back(to);
                    continue;
                }
                if (*p != to) {
                    connected.insert(p, to);
                }
            }
        }
        return connected;
    }
} // namespace zx
