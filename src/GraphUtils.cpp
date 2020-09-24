#include "GraphUtils.hpp"

/**
 * Constructs a list of nodes from the given <code>nodes</code> list at the
 * given indices in that list.
 *
 * @param indices The indices of the desired nodes in <code>nodes</code>.
 * @param nodes The list of nodes from which we select a sublist.
 * @return the The sublist selected.
 */
std::vector<Variable*> GraphUtils::asList(std::vector<int>& indices, std::vector<Variable*>& nodes) {
    std::vector<Variable*> list;

    for (int i : indices) {
        list.push_back(nodes[i]);
    }

    return list;
}

std::unordered_set<Variable*> GraphUtils::asSet(std::vector<int>& indices, std::vector<Variable*>& nodes) {
    std::unordered_set<Variable*> set;

    for (int i : indices) {
        set.insert(nodes[i]);
    }

    return set;
}

EdgeListGraph GraphUtils::completeGraph(EdgeListGraph& graph) {
    EdgeListGraph graph2(graph.getNodes());

    graph2.removeEdges();

    std::vector<Variable*> nodes = graph2.getNodes();

    for (int i = 0; i < nodes.size(); i++) {
        for (int j = i+1; j < nodes.size(); j++) {
            Variable* node1 = nodes[i];
            Variable* node2 = nodes[j];
            graph2.addUndirectedEdge(node1, node2);
        }
    }

    return graph2;
}

EdgeListGraph GraphUtils::undirectedGraph(EdgeListGraph& graph) {
    EdgeListGraph graph2(graph.getNodes());

    for (Edge edge : graph.getEdges()) {
        graph2.addUndirectedEdge(edge.getNode1(), edge.getNode2());
    }

    return graph2;
}
