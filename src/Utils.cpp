#include "Utils.hpp"

namespace zx {
    Vertices::VertexIterator::VertexIterator(
            std::vector<std::optional<VertexData>>& vertices, Vertex v):
        v(v),
        currentPos(vertices.begin()), vertices(vertices) {
        if ((size_t)v >= vertices.size()) {
            currentPos = vertices.end();
            v          = vertices.size();
        } else {
            currentPos = vertices.begin() + static_cast<int>(v);
            next_valid_vertex();
        }
    }
    // Prefix increment
    Vertices::VertexIterator Vertices::VertexIterator::operator++() {
        Vertices::VertexIterator it = *this;
        currentPos++;
        v++;
        next_valid_vertex();
        return it;
    }

    // Postfix increment
    const Vertices::VertexIterator Vertices::VertexIterator::operator++(int) {
        currentPos++;
        v++;
        next_valid_vertex();
        return *this;
    }

    bool operator==(const Vertices::VertexIterator& a,
                    const Vertices::VertexIterator& b) {
        return a.currentPos == b.currentPos;
    }
    bool operator!=(const Vertices::VertexIterator& a,
                    const Vertices::VertexIterator& b) {
        return !(a == b);
    }

    void Vertices::VertexIterator::next_valid_vertex() {
        while (currentPos != vertices.end() && !currentPos->has_value()) {
            v++;
            currentPos++;
        }
    }

    Edges::EdgeIterator::EdgeIterator(
            std::vector<std::vector<Edge>>&         edges,
            std::vector<std::optional<VertexData>>& vertices):
        v(0),
        currentPos(edges[0].begin()), edges(edges), vertices(vertices) {
        if (!vertices.empty()) {
            while ((size_t)v < edges.size() && !vertices[v].has_value())
                v++;
            currentPos = edges[v].begin();
            checkNextVertex();
        } else {
            currentPos = edges.back().end();
            v          = edges.size();
        }
    }

    Edges::EdgeIterator::EdgeIterator(
            std::vector<std::vector<Edge>>&         edges,
            std::vector<std::optional<VertexData>>& vertices, Vertex v):
        v(v),
        edges(edges), vertices(vertices) {
        if ((size_t)v >= edges.size()) {
            currentPos = edges.back().end();
            v          = edges.size();
        } else {
            currentPos = edges[v].begin();
        }
    }

    // Prefix increment
    Edges::EdgeIterator Edges::EdgeIterator::operator++() {
        Edges::EdgeIterator it = *this;
        currentPos++;
        checkNextVertex();
        return it;
    }

    void Edges::EdgeIterator::checkNextVertex() {
        while (currentPos != edges[v].end() &&
               currentPos->to < v) // make sure to not iterate over an edge twice
            currentPos++;

        while (currentPos == edges[v].end() && (size_t)v < edges.size()) {
            v++;
            while ((size_t)v < edges.size() && !vertices[v].has_value())
                v++;

            if ((size_t)v == edges.size()) {
                currentPos = edges.back().end();
                v--;
                return;
            }
            currentPos = edges[v].begin();
            while (currentPos != edges[v].end() &&
                   currentPos->to < v) // make sure to not iterate over an edge twice
                currentPos++;
        }
    }
    // Postfix increment
    const Edges::EdgeIterator Edges::EdgeIterator::operator++(int) {
        currentPos++;
        checkNextVertex();
        return *this;
    }

    bool operator==(const Edges::EdgeIterator& a, const Edges::EdgeIterator& b) {
        return a.currentPos == b.currentPos;
    }
    bool operator!=(const Edges::EdgeIterator& a, const Edges::EdgeIterator& b) {
        return !(a == b);
    }

    bool isPauli(const PiExpression& expr) {
        return expr.isConstant() && expr.getConst().isInteger();
    }
    bool isClifford(const PiExpression& expr) {
        return expr.isConstant() && (expr.getConst().isInteger() || expr.getConst().getDenom() == 2);
    }
    bool isProperClifford(const PiExpression& expr) {
        return expr.isConstant() && expr.getConst().getDenom() == 2;
    }

    void roundToClifford(PiExpression& expr, fp tolerance) {
        if (!expr.isConstant())
            return;

        if (expr.getConst().isCloseDivPi(0, tolerance)) {
            expr.setConst(PiRational(0, 1));
        } else if (expr.getConst().isCloseDivPi(0.5, tolerance)) {
            expr.setConst(PiRational(1, 2));
        } else if (expr.getConst().isCloseDivPi(-0.5, tolerance)) {
            expr.setConst(PiRational(-1, 2));
        } else if (expr.getConst().isCloseDivPi(1, tolerance)) {
            expr.setConst(PiRational(1, 1));
        }
    }
} // namespace zx
