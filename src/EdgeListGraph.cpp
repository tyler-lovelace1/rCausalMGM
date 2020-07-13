#include "EdgeListGraph.hpp"


// Used by constructors
void EdgeListGraph::initNamesHash() {
    for (Variable* node: nodes) {
        namesHash[node->getName()] = node;
    }
}

/**
 * Constructs a new (empty) EdgeListGraph.
 */
EdgeListGraph::EdgeListGraph() {
    initNamesHash();
}

/**
 * Constructs a new graph, with no edges, using the the given variable
 * names.
 */
EdgeListGraph::EdgeListGraph(const std::vector<Variable*>& nodes) {
    for (Variable* variable : nodes) {
        if(!addNode(variable))
            throw std::invalid_argument("Issue adding variable to graph");
    }

    initNamesHash();
}

/**
 * Constructs a EdgeListGraph using the nodes and edges of the given graph.
 * If this cannot be accomplished successfully, an exception is thrown. Note
 * that any graph constraints from the given graph are forgotten in the new
 * graph.
 *
 * @param graph the graph from which nodes and edges are is to be
 *              extracted.
 * @throws std::invalid_argument if a duplicate edge is added.
 */
EdgeListGraph::EdgeListGraph(const EdgeListGraph& graph) {
    transferNodesAndEdges(graph);
    ambiguousTriples = graph.ambiguousTriples;
    underLineTriples = graph.underLineTriples;
    dottedUnderLineTriples = graph.dottedUnderLineTriples;

    for (Edge edge : graph.edgesSet) {
        if (graph.highlightedEdges.count(edge)) {
            setHighlighted(edge, true);
        }
    }

    namesHash = graph.namesHash;

    pag = graph.pag;
    pattern = graph.pattern;
}

/**
 * Transfers nodes and edges from one graph to another.  One way this is
 * used is to change graph types.  One constructs a new graph based on the
 * old graph, and this method is called to transfer the nodes and edges of
 * the old graph to the new graph.
 *
 * @param graph the graph from which nodes and edges are to be pilfered.
 * @throws IllegalArgumentException This exception is thrown if adding some
 *                                  node or edge violates one of the
 *                                  basicConstraints of this graph.
 */
void EdgeListGraph::transferNodesAndEdges(const EdgeListGraph& graph) {
    for (Variable* node : graph.nodes) {
        if (!addNode(node))
            throw std::invalid_argument("Problem copying graph nodes");
    }

    for (Edge edge : graph.edgesSet) {
        if (!addEdge(edge)) 
            throw std::invalid_argument("Problem copying edges");
    }
}

/**
 * Adds an undirected edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addUndirectedEdge(Variable* node1, Variable* node2) {
    Edge newEdge = Edge::undirectedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds an edge to the graph.
 *
 * @param edge the edge to be added
 * @return true if the edge was added, false if not.
 */
bool EdgeListGraph::addEdge(Edge& edge) {

    auto itr1 = edgeLists.find(edge.getNode1());
    auto itr2 = edgeLists.find(edge.getNode2());

    // Do not comment this out; if the user changes the names of variables, this is the
    // mechanism for adjusting the maps from nodes to edge lists to compensate.
    if (itr1 == edgeLists.end() || itr2 == edgeLists.end()) {
        initNamesHash();
        itr1 = edgeLists.find(edge.getNode1());
        itr2 = edgeLists.find(edge.getNode2());
    }

    if (itr1 == edgeLists.end() || itr2 == edgeLists.end()) {
        // Convert edge to string
        std::ostringstream ss;
        ss << edge;
        throw std::invalid_argument("Can't add an edge unless both nodes are in the graph: " + ss.str());
    }

    std::vector<Edge> edgeList1 = itr1->second;
    std::vector<Edge> edgeList2 = itr2->second;

    // If the edge is already in the graph
    if ((std::find(edgeList1.begin(), edgeList1.end(), edge) != edgeList1.end()) ||
        (std::find(edgeList2.begin(), edgeList2.end(), edge) != edgeList2.end()))
        return true;
    
    edgeList1.push_back(edge);
    edgeList2.push_back(edge);

    edgeLists[edge.getNode1()] = edgeList1;
    edgeLists[edge.getNode2()] = edgeList2;

    edgesSet.insert(edge);

    return true;
}

/**
 * Adds a node to the graph. Precondition: The proposed name of the node
 * cannot already be used by any other node in the same graph.
 *
 * @param node the node to be added.
 * @return true if the the node was added, false if not.
 */
bool EdgeListGraph::addNode(Variable* node) {
    // If nodes contains node
    if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) return true;

    if (node == NULL)
        throw std::invalid_argument("Can't add NULL node to graph");

    if (!(getNode(node->getName()) == NULL)) {
        if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) {
            namesHash[node->getName()] = node;
        }
    }

    if (edgeLists.count(node)) return false;

    edgeLists[node] = {};
    nodes.push_back(node);
    namesHash[node->getName()] = node;

    return true;
}

/**
 * Removes any relevant edge objects found in this collection. G
 *
 * @param edges the collection of edges to remove.
 * @return true if any edges in the collection were removed, false if not.
 */
bool EdgeListGraph::removeEdges(const std::vector<Edge>& edges) {
    bool change = false;

    for (Edge edge : edges) {
        bool _change = removeEdge(edge);
        change = change || _change;
    }

    return change;
}

/**
 * Removes all edges connecting node A to node B.
 *
 * @param node1 the first node.,
 * @param node2 the second node.
 * @return true if edges were removed between A and B, false if not.
 */
bool EdgeListGraph::removeEdges(Variable* node1, Variable* node2) {
    return removeEdges(getEdges(node1, node2));
}

/**
 * Removes an edge from the graph. (Note: It is dangerous to make a
 * recursive call to this method (as it stands) from a method containing
 * certain types of iterators. The problem is that if one uses an iterator
 * that iterates over the edges of node A or node B, and tries in the
 * process to remove those edges using this method, a concurrent
 * modification exception will be thrown.)
 *
 * @param edge the edge to remove.
 * @return true if the edge was removed, false if not.
 */
bool EdgeListGraph::removeEdge(Edge& edge) {
    if (!edgesSet.count(edge)) return false;

    std::vector<Edge> edgeList1 = edgeLists[edge.getNode1()];
    std::vector<Edge> edgeList2 = edgeLists[edge.getNode2()];

    edgesSet.erase(edge);
    edgeList1.erase(std::remove(edgeList1.begin(), edgeList1.end(), edge), edgeList1.end());
    edgeList2.erase(std::remove(edgeList2.begin(), edgeList2.end(), edge), edgeList2.end());

    edgeLists[edge.getNode1()] = edgeList1;
    edgeLists[edge.getNode2()] = edgeList2;

    highlightedEdges.erase(edge);
    stuffRemovedSinceLastTripleAccess = true;

    return true;
}

/**
 * Removes the edge connecting the two given nodes.
 */
bool EdgeListGraph::removeEdge(Variable* node1, Variable* node2) {
    std::vector<Edge> edges = getEdges(node1, node2);

    if (edges.size() > 1) {
        throw std::invalid_argument("There is more than one edge between " + node1->getName() + " and " + node2->getName());
    }

    return removeEdges(edges);
}

/**
 * Determines whether some edge or other exists between two nodes.
 */
bool EdgeListGraph::isAdjacentTo(Variable* node1, Variable* node2) {
    if (node1 == NULL || node2 == NULL || 
        edgeLists.find(node1) == edgeLists.end() ||
        edgeLists.find(node2) == edgeLists.end()) {
            return false;
    }

    for (Edge edge : edgeLists[node1]) {
        if (Edge::traverse(node1, edge) == node2) return true;
    }

    return false;
}

bool EdgeListGraph::isDirectedFromTo(Variable* node1, Variable* node2) {
    std::vector<Edge> edges = getEdges(node1, node2);
    if (edges.size() != 1) return false;
    Edge edge = edges[0];
    return edge.pointsTowards(node2);
}

bool EdgeListGraph::isUndirectedFromTo(Variable* node1, Variable* node2) {

    Edge edge;
    try {
        edge = getEdge(node1, node2);
    } catch (std::invalid_argument& e) {
        return false;
    }
    
    return edge.getEndpoint1() == ENDPOINT_TAIL && edge.getEndpoint2() == ENDPOINT_TAIL;
}

bool EdgeListGraph::isAmbiguousTriple(Variable* x, Variable* y, Variable* z) {
    return ambiguousTriples.count(Triple(x, y, z));
}

/**
 * @return the set of nodes adjacent to the given node. If there are multiple edges between X and Y, Y will show
 * up twice in the list of adjacencies for X, for optimality; simply create a list an and array from these to
 * eliminate the duplication.
 */
std::vector<Variable*> EdgeListGraph::getAdjacentNodes(Variable* node) {
    std::vector<Edge> edges = edgeLists[node];
    std::unordered_set<Variable*> adj;

    for (Edge edge : edges) {
        Variable* z = edge.getDistalNode(node);
        adj.insert(z);
    }

    return std::vector<Variable*>(adj.begin(), adj.end());

}

/**
 * @return the edge connecting node1 and node2, provided a unique such edge
 * exists.
 * 
 * Throws std::invalid_argument if not
 */
Edge EdgeListGraph::getEdge(Variable* node1, Variable* node2) {

    std::vector<Edge> edges = edgeLists[node1];

    if (edges.size() == 0)
        throw std::invalid_argument("No edges coming from node1");

    for (Edge edge : edges) {
        if (edge.getNode1() == node1 && edge.getNode2() == node2) {
            return edge;
        } else if (edge.getNode1() == node2 && edge.getNode2() == node1) {
            return edge;
        }
    }

    throw std::invalid_argument("node1 and node2 not connected by edge");

}

/**
 * @return the edges connecting node1 and node2.
 */
std::vector<Edge> EdgeListGraph::getEdges(Variable* node1, Variable* node2) {
    std::vector<Edge> edges = edgeLists[node1];
    std::vector<Edge> _edges;

    for (Edge edge : edges) {
        if (edge.getDistalNode(node1) == node2) {
            _edges.push_back(edge);
        }
    }

    return _edges;
}

/**
 * @return the endpoint along the edge from node to node2 at the node2 end.
 */
Endpoint EdgeListGraph::getEndpoint(Variable* node1, Variable* node2) {
    std::vector<Edge> edges = getEdges(node2);

    for (Edge edge : edges) {
        if (edge.getDistalNode(node2) == node1) return edge.getProximalEndpoint(node2);
    }

    return ENDPOINT_NULL;
}

/**
 * @return the list of parents for a node.
 */
std::vector<Variable*> EdgeListGraph::getParents(Variable* node) {
    std::vector<Variable*> parents;
    std::vector<Edge> edges = edgeLists[node];

    for (Edge edge : edges) {
        Endpoint endpoint1 = edge.getDistalEndpoint(node);
        Endpoint endpoint2 = edge.getProximalEndpoint(node);

        if (endpoint1 == ENDPOINT_TAIL && endpoint2 == ENDPOINT_ARROW) {
            parents.push_back(edge.getDistalNode(node));
        }
    }

    return parents;
}

/**
 * If there is currently an edge from node1 to node2, sets the endpoint at
 * node2 to the given endpoint; if there is no such edge, adds an edge --#
 * where # is the given endpoint. Setting an endpoint to null, provided
 * there is exactly one edge connecting the given nodes, removes the edge.
 * (If there is more than one edge, an exception is thrown.)
 *
 * @throws std::invalid_argument if the edge with the revised endpoint
 *                                  cannot be added to the graph.
 */
bool EdgeListGraph::setEndpoint(Variable* from, Variable* to, Endpoint endPoint) {
    std::vector<Edge> edges = getEdges(from, to);

    if (endPoint == ENDPOINT_NULL)
        throw std::invalid_argument("Endpoint cannot be NULL");

    if (edges.size() == 0) {
        removeEdges(from, to);
        Edge newEdge(from, to, ENDPOINT_TAIL, endPoint);
        addEdge(newEdge);
        return true;
    }

    Edge edge = edges[0];
    Edge newEdge(from, to, edge.getProximalEndpoint(from), endPoint);

    try {
        removeEdges(edge.getNode1(), edge.getNode2());
        addEdge(newEdge);
        return true;
    } catch (std::invalid_argument& e) {
        return false;
    }

    return false; // Unreachable
}

std::ostream& operator<<(std::ostream& os, EdgeListGraph& graph) {
    os << "edges(g) = \n";
    for (Edge edge : graph.edgesSet) {
        os << edge << "\n";
    }

    return os;
}
