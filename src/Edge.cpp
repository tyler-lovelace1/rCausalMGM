#include "Edge.hpp"

struct std::hash<Edge> {
    std::size_t operator()(const Edge& k) const {
        return std::hash<Variable*>()(k.node1) + std::hash<Variable*>()(k.node2);
    }
};

bool operator< (const Edge& e1, const Edge& e2) {
    if (e1.node1 != e2.node1) {
        return e1.node1 < e2.node1;
    }

    return e1.node2 < e2.node2;
}

bool operator>= (const Edge& e1, const Edge& e2) {
    return !(e1 < e2);
}

bool operator<= (const Edge& e1, const Edge& e2) {
    return (e1 < e2) || (e1 == e2);
}

bool operator> (const Edge& e1, const Edge& e2) {
    return !(e1 <= e2);
}

void Edge::sortEdges(std::vector<Edge> edges) {
    std::sort(edges.begin(), edges.end());
}