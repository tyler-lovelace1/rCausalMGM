#include "MeekRules.hpp"

#include "GraphUtils.hpp"

void MeekRules::orientImplied(EdgeListGraph& graph) {
    nodes = graph.getNodes();
    orientImplied(graph, nodes);
}

void MeekRules::orientImplied(EdgeListGraph& graph, std::vector<Node>& nodes) {
    this->nodes = nodes;
    visited.insert(nodes.begin(), nodes.end());

    orientUsingMeekRulesLocally(graph);
}

void MeekRules::orientUsingMeekRulesLocally(EdgeListGraph& graph) {
    
    oriented = {};

    if (shouldUndirectUnforcedEdges) {
        for (Node node : nodes) {
            undirectUnforcedEdges(node, graph);
            std::vector<Node> adjacentNodes = graph.getAdjacentNodes(node);
            for (Node adj : adjacentNodes) {
                directStack.push(adj);
            }
        }
    }

    for (Node node : nodes) {
        runMeekRules(node, graph);
    }

    while(!directStack.empty()) {
        Node node = directStack.top();
        directStack.pop();

        if (shouldUndirectUnforcedEdges) {
            undirectUnforcedEdges(node, graph);
        }

        runMeekRules(node, graph);
    }
}

void MeekRules::runMeekRules(Node node, EdgeListGraph& graph) {
    meekR1(node, graph);
    meekR2(node, graph);
    meekR3(node, graph);
    meekR4(node, graph);
}

// void MeekRules::meekR0(EdgeListGraph& graph) {
//     graph.reorientAllWith(ENDPOINT_TAIL);
//     SearchGraphUtils::pcOrientbk(knowledge, graph);
//     std::vector<Triple> orderedColliders = sepsets.getOrderedColliders();

//     for (Triple triple : orderedColliders) {
//         Node a = triple.getX();
//         Node b = triple.getY();
//         Node c = triple.getZ();

//         if (!(graph.getEndpoint(b, a) == ENDPOINT_ARROW
// 	      || graph.getEndpoint(b, c) == ENDPOINT_ARROW)) {
//             SearchGraphUtils::orientCollider(sepsets, graph, a, b, c);
//         }
//     }
// }

/**
 * Meek's rule R1: if a-->b, b---c, and a not adj to c, then a-->c
 */
void MeekRules::meekR1(Node b, EdgeListGraph& graph) {
    std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);

    if (adjacentNodes.size() < 2) {
        return;
    }

    ChoiceGenerator cg(adjacentNodes.size(), 2);
    std::vector<int> *choice;

    for (choice = cg.next(); choice != NULL; choice = cg.next()) {
        std::vector<Node> nodes = GraphUtils::asList(*choice, adjacentNodes);
        Node a = nodes[0];
        Node c = nodes[1];

        r1Helper(a, b, c, graph);
        r1Helper(c, b, a, graph);

    }

}

void MeekRules::r1Helper(Node a, Node b, Node c, EdgeListGraph& graph) {
    if (!graph.isAdjacentTo(a, c) && graph.isDirectedFromTo(a, b) && graph.isUndirectedFromTo(b, c)) {
        if (!isUnshieldedNoncollider(a, b, c, graph)) {
            return;
        }

        if (isArrowpointAllowed(b, c) && doesNotCreateCycle(b, c, graph)) {
            direct(b, c, graph);
            Edge edge = graph.getEdge(b, c);
            // Rcpp::Rcout << "Meek R1 edge: " << edge << std::endl;
            // Rcpp::Rcout << "Meek R1 triangle (" << b->getName() << "-->" << a->getName() << "---" << c->getName() << ")" << std::endl;
        }
    }
}

/**
 * If a-->b-->c, a--c, then b-->c.
 */
void MeekRules::meekR2(Node c, EdgeListGraph& graph) {
    std::vector<Node> adjacentNodes = graph.getAdjacentNodes(c);

    if (adjacentNodes.size() < 2) {
        return;
    }

    ChoiceGenerator cg(adjacentNodes.size(), 2);
    std::vector<int> *choice;

    for (choice = cg.next(); choice != NULL; choice = cg.next()) {
        std::vector<Node> nodes = GraphUtils::asList(*choice, adjacentNodes);
        Node a = nodes[0];
        Node b = nodes[1];

        r2Helper(a, b, c, graph);
        r2Helper(b, a, c, graph);
        r2Helper(a, c, b, graph);
        r2Helper(c, a, b, graph);
    }
}

void MeekRules::r2Helper(Node a, Node b, Node c, EdgeListGraph& graph) {
    if (graph.isDirectedFromTo(a, b) && graph.isDirectedFromTo(b, c) && graph.isUndirectedFromTo(a, c)) {
        if (isArrowpointAllowed(a, c) && doesNotCreateCycle(a, c, graph)) {
            direct(a, c, graph);
            Edge edge = graph.getEdge(a, c);
            // Rcpp::Rcout << "Meek R2 Edge: " << edge << std::endl;
        }
    }
}

/**
 * Meek's rule R3. If a--b, a--c, a--d, c-->b, d-->b, then orient a-->b.
 */
void MeekRules::meekR3(Node a, EdgeListGraph& graph) {
    std::vector<Node> adjacentNodes = graph.getAdjacentNodes(a);

    if (adjacentNodes.size() < 3) {
        return;
    }

    for (Node d : adjacentNodes) {
        if (Edge::isUndirectedEdge(graph.getEdge(a, d))) {
            std::vector<Node> otherAdjacents(adjacentNodes);
            otherAdjacents.erase(std::remove(otherAdjacents.begin(), otherAdjacents.end(), d), otherAdjacents.end());

            ChoiceGenerator cg(otherAdjacents.size(), 2);
            std::vector<int> *choice;

            for (choice = cg.next(); choice != NULL; choice = cg.next()) {
                std::vector<Node> nodes = GraphUtils::asList(*choice, otherAdjacents);
                Node b = nodes[0];
                Node c = nodes[1];

                if (isKite(a, d, b, c, graph)) {
                    if (isArrowpointAllowed(d, a) && doesNotCreateCycle(d, a, graph)) {
                        if (!isUnshieldedNoncollider(c, d, b, graph)) {
                            continue;
                        }

                        direct(d, a, graph);
                        Edge edge = graph.getEdge(d, a);
                        // Rcpp::Rcout << "Meek R3 Edge: " << edge << std::endl;
                    }
                }
            }

        }
    }
}

bool MeekRules::isKite(Node a, Node d, Node b, Node c, EdgeListGraph& graph) {
    return graph.isUndirectedFromTo(d, c) &&
           graph.isUndirectedFromTo(d, b) &&
           graph.isDirectedFromTo(b, a) &&
           graph.isDirectedFromTo(c, a) &&
           graph.isUndirectedFromTo(d, a);
}

void MeekRules::meekR4(Node a, EdgeListGraph& graph) {
    if (!useRule4) {
        return;
    }

    // TODO?
    Rcpp::Rcout << "Since rule4 requires knowlegde, this function should never be used" << std::endl;
}

void MeekRules::direct(Node a, Node c, EdgeListGraph& graph) {
    Edge before = graph.getEdge(a, c);

    Edge after = Edge::directedEdge(a, c);

    visited.insert(a);
    visited.insert(c);

    graph.removeEdge(before);
    graph.addEdge(after);

    oriented.insert(after);

    directStack.push(c);
}

bool MeekRules::isUnshieldedNoncollider(Node a, Node b, Node c, EdgeListGraph& graph) {
    if (!graph.isAdjacentTo(a, b)) {
        return false;
    }

    if (!graph.isAdjacentTo(c, b)) {
        return false;
    }

    if (graph.isAdjacentTo(a, c)) {
        return false;
    }

    if (graph.isAmbiguousTriple(a, b, c)) {
        return false;
    }

    return !(graph.getEndpoint(a, b) == ENDPOINT_ARROW &&
             graph.getEndpoint(c, b) == ENDPOINT_ARROW);
}

bool MeekRules::isArrowpointAllowed(Node from, Node to) {
    // Knowledge
    if (knowledge.isEmpty())
	return true;
    return !knowledge.isForbidden(from, to) && !knowledge.isRequired(to, from);
}

bool MeekRules::doesNotCreateCycle(Node from, Node to, EdgeListGraph& graph) {
    if (aggressivelyPreventCycles) {
	return !graph.isAncestorOf(to, from);
    } else {
	return true;
    }
}

void MeekRules::undirectUnforcedEdges(Node y, EdgeListGraph& graph) {
    std::unordered_set<Node> parentsToUndirect;
    std::vector<Node> parents = graph.getParents(y);

    for (Node x : parents) {
        for (Node parent : parents) {
            if (parent != x) {
                if (!graph.isAdjacentTo(parent, x)) {
                    oriented.insert(graph.getEdge(x, y));
                    goto NEXT_EDGE;
                }
            }
        }

        parentsToUndirect.insert(x);
        NEXT_EDGE:;
    }

    bool didit = false;

    for (Node x : parentsToUndirect) {
	bool mustOrient = knowledge.isRequired(x,y) || knowledge.isForbidden(y,x);

        if (!oriented.count(graph.getEdge(x, y)) && !mustOrient) {
            graph.removeEdge(x, y);
            graph.addUndirectedEdge(x, y);
            visited.insert(x);
            visited.insert(y);
            didit = true;
        }
    }

    if (didit) {
        for (Node z : graph.getAdjacentNodes(y)) {
            directStack.push(z);
        }

        directStack.push(y);
    }
}

