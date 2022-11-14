#pragma once

#include "Definitions.hpp"
#include "Expression.hpp"

#include <cstdint>
#include <iterator>
#include <optional>
#include <utility>
#include <vector>

namespace zx {

    struct Edge {
        Vertex   to;
        EdgeType type;

        Edge() = default;
        Edge(const Vertex to, const EdgeType type):
            to(to), type(type){};
        void toggle() {
            if (type == EdgeType::Simple) {
                type = EdgeType::Hadamard;
            } else {
                type = EdgeType::Simple;
            }
        }
    };

    struct VertexData {
        Col          col;
        Qubit        qubit;
        PiExpression phase;
        VertexType   type;
    };

    class Vertices {
    public:
        explicit Vertices(
                const std::vector<std::optional<VertexData>>& vertices):
            vertices(vertices){};

        class VertexIterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::int32_t;
            using value_type        = std::pair<Vertex, const VertexData&>;
            using pointer           = value_type*;
            using reference         = value_type&;

            explicit VertexIterator(const std::vector<std::optional<VertexData>>& vertices):
                currentPos(vertices.begin()), vertices(vertices) {
                nextValidVertex();
            }
            VertexIterator(const std::vector<std::optional<VertexData>>& vertices,
                           Vertex                                        v);

            value_type operator*() const { return {v, currentPos->value()}; }
            // pointer operator->() { return ptr; }

            // Prefix increment
            VertexIterator operator++();

            // Postfix increment
            const VertexIterator operator++(int);

            friend bool operator==(const VertexIterator& a,
                                   const VertexIterator& b);
            friend bool operator!=(const VertexIterator& a,
                                   const VertexIterator& b);

        private:
            Vertex                                                 v = 0;
            std::vector<std::optional<VertexData>>::const_iterator currentPos;
            const std::vector<std::optional<VertexData>>&          vertices;

            void nextValidVertex();
        };

        using iterator = VertexIterator;

        const iterator begin() { return VertexIterator(vertices); }
        const iterator end() { return {vertices, vertices.size()}; }

    private:
        const std::vector<std::optional<VertexData>>& vertices;
    };

    class Edges {
    public:
        Edges(const std::vector<std::vector<Edge>>&         edges,
              const std::vector<std::optional<VertexData>>& vertices):
            edges(edges),
            vertices(vertices){};

        class EdgeIterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::int32_t;
            using value_type        = std::pair<Vertex, Vertex>;
            using pointer           = value_type*;
            using reference         = value_type&;

            EdgeIterator(const std::vector<std::vector<Edge>>&         edges,
                         const std::vector<std::optional<VertexData>>& vertices);

            EdgeIterator(const std::vector<std::vector<Edge>>&         edges,
                         const std::vector<std::optional<VertexData>>& vertices, Vertex v);

            value_type operator*() const { return {v, currentPos->to}; }
            // pointer operator->() { return ptr; }

            // Prefix increment
            EdgeIterator operator++();

            // Postfix increment
            const EdgeIterator operator++(int);

            friend bool operator==(const EdgeIterator& a, const EdgeIterator& b);
            friend bool operator!=(const EdgeIterator& a, const EdgeIterator& b);

        private:
            Vertex                                         v;
            std::vector<Edge>::const_iterator              currentPos;
            std::vector<std::vector<Edge>>::const_iterator edgesPos;
            const std::vector<std::vector<Edge>>&          edges;
            const std::vector<std::optional<VertexData>>&  vertices;

            void checkNextVertex();
        };

        using iterator = EdgeIterator;

        const iterator begin() { return {edges, vertices}; }
        const iterator end() { return {edges, vertices, edges.size()}; }

    private:
        const std::vector<std::vector<Edge>>&         edges;
        const std::vector<std::optional<VertexData>>& vertices;
    };

    bool isPauli(const PiExpression& expr);
    bool isClifford(const PiExpression& expr);
    bool isProperClifford(const PiExpression& expr);

    void roundToClifford(PiExpression& expr, fp tolerance);
} // namespace zx
