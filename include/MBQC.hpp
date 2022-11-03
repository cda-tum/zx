#pragma once

#include "Definitions.hpp"
#include "Utils.hpp"
#include "ZXDiagram.hpp"

#include <cmath>
#include <initializer_list>
#include <optional>
#include <vector>

namespace zx {

    class MBQCPattern {
    public:
        explicit MBQCPattern(const ZXDiagram& diag);
        MBQCPattern() = default;
        MBQCPattern(const std::initializer_list<Measurement>& measurements, const std::initializer_list<Vertex>& inputs, const std::initializer_list<Vertex>& outputs, const std::initializer_list<std::pair<Vertex, Vertex>>& edges);

        //implements the polynomial algorithm for finding a maximally delayed g-flow from https://arxiv.org/pdf/0709.2670.pdf
        [[nodiscard]] std::optional<flow> computeGFlow() const;
        [[nodiscard]] std::optional<flow> computeCausalFlow() const; // TODO

        [[nodiscard]] std::size_t getNVertices() const { return measurements.size(); }
        [[nodiscard]] std::size_t getNEdges() const { return nedges; }
        [[nodiscard]] std::size_t getNQubits() const { return inputs.size(); }

        [[nodiscard]] Measurement getMeasurement(const Vertex& v) const { return measurements[v]; }
        [[nodiscard]] Measurement operator[](std::size_t i) const { return measurements[i]; }

        [[nodiscard]] std::vector<Vertex> getInputs() const { return inputs; }
        [[nodiscard]] std::vector<Vertex> getOutputs() const { return outputs; }

        Vertex addQubit(const Measurement& meas) {
            measurements.emplace_back(meas);
            return measurements.size() - 1;
        }
        Vertex addInputQubit(Vertex v) {
            inputs.emplace_back(v);
            return v;
        }
        Vertex addInputQubit(const Measurement& meas) {
            const auto v = addQubit(meas);
            return addInputQubit(v);
        }
        Vertex addOutputQubit(Vertex v) {
            outputs.emplace_back(v);
            return v;
        }
        Vertex addOutputQubit(const Measurement& meas) {
            const auto v = addQubit(meas);
            return addOutputQubit(v);
        }
        void addEdge(const Vertex from, const Vertex to) {
            if (!isIn(from, edges[to])) {
                edges[from].emplace_back(to);
                edges[to].emplace_back(from);
                nedges++;
            }
        }

        [[nodiscard]] std::vector<Vertex> getNeighbourhood(const Vertex v) const {
            return edges[v];
        }

    private:
        std::vector<std::vector<Vertex>> edges{};
        std::vector<Measurement>         measurements{};
        std::vector<Vertex>              inputs{};
        std::vector<Vertex>              outputs{};
        std::size_t                      nedges = 0;

        gf2Mat                                              getAdjMat() const;
        std::pair<std::vector<Vertex>, std::vector<Vertex>> getNonProcessed(const std::vector<Vertex>& out) const;

        static bool isIn(const Vertex& v, const std::vector<Vertex>& vertices);

        gf2Mat constructLinearSystem(const gf2Mat& adjMat, std::vector<Vertex> out, std::vector<Vertex> out_prime, std::vector<Vertex> us) const;

        std::vector<Vertex> solutionFromTriangular(const gf2Mat& triu, const std::vector<Vertex>& us, std::size_t offset, const std::vector<Vertex>& out_prime, std::vector<std::vector<Vertex>>& g) const;
        std::vector<Vertex> getConnectedSet(const std::vector<Vertex> s) const;

        static std::vector<Vertex> fromIdxVec(const std::vector<bool>& indicator, const std::vector<Vertex>& set);
    };
} // namespace zx
