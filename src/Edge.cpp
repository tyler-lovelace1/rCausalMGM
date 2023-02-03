#include "Edge.hpp"

// const Node& Edge::nullconst Node& = Node();

/**
 * Constructs a new edge by specifying the nodes it connects and the
 * endpoint types.
 *
 * @param node1     the first node
 * @param node2     the second node            _
 * @param endpoint1 the endpoint at the first node
 * @param endpoint2 the endpoint at the second node
 */
Edge::Edge(const Node& node1, const Node& node2, Endpoint endpoint1, Endpoint endpoint2) {
    init(node1, node2, endpoint1, endpoint2);
}

Edge::Edge(const Edge& edge) {
    init(edge.node1, edge.node2, edge.endpoint1, edge.endpoint2);
}

void Edge::init(const Node& node1, const Node& node2, Endpoint endpoint1, Endpoint endpoint2) {
    if (node1.isNull() || node2.isNull()) {
        throw std::invalid_argument("Nodes must not be null");
    }

    // Flip edges pointing left the other way.
    if(pointingLeft(endpoint1, endpoint2) || (endpoint1==endpoint2 && node2 < node1)) {
        this->node1 = node2;
        this->node2 = node1;
        this->endpoint1 = endpoint2;
        this->endpoint2 = endpoint1;
    } else {
        this->node1 = node1;
        this->node2 = node2;
        this->endpoint1 = endpoint1;
        this->endpoint2 = endpoint2;
    }
}

/**
 * @return the endpoint nearest to the given node.
 * @throws std::invalid_argument if the given node is not along the
 *                                  edge.
 */
Endpoint Edge::getProximalEndpoint(const Node& node) {
    if (node1 == node) {
        return endpoint1;
    } else if (node2 == node) {
        return endpoint2;
    }

    throw std::invalid_argument("Given node must be one along the edge");
}

/**
 * @return the endpoint furthest from the given node.
 * @throws std::invalid_argument if the given node is not along the
 *                                  edge.
 */
Endpoint Edge::getDistalEndpoint(const Node& node) {
    if (node1 == node) {
        return endpoint2;
    } else if (node2 == node) {
        return endpoint1;
    }

    throw std::invalid_argument("Given node must be one along the edge");
}

/**
 * Traverses the edge in an undirected fashion--given one node along the
 * edge, returns the node at the opposite end of the edge.
 * 
 * Returns null if the given node is not part of the edge 
 */
Node Edge::getDistalNode(const Node& node) {
    if (node1 == node) {
        return node2;
    }

    if (node2 == node) {
        return node1;
    }

    return Node();
}

/**
 * @return true just in case this edge is directed.
 */
bool Edge::isDirected() {
    return Edge::isDirectedEdge(*this);
}

/**
 * @return true just in case this edge is undirected.
 */
bool Edge::isUndirected() {
    return Edge::isUndirectedEdge(*this);
}

/**
 * @return true just in case this edge is bidirected.
 */
bool Edge::isBidirected() {
    return Edge::isBidirectionalEdge(*this);
}

/**
 * @return true just in case this edge is nondirected.
 */
bool Edge::isNondirected() {
    return Edge::isNondirectedEdge(*this);
}

/**
 * @return true just in case this edge is partially oriented.
 */
bool Edge::isPartiallyOriented() {
    return Edge::isPartiallyOrientedEdge(*this);
}

/**
 * @return true just in case the edge is pointing toward the given node--
 * that is, x --> node or x o--> node.
 */
bool Edge::pointsTowards(const Node& node) {
    Endpoint proximal = getProximalEndpoint(node);
    Endpoint distal = getDistalEndpoint(node);
    return (proximal == ENDPOINT_ARROW && (distal == ENDPOINT_TAIL || distal == ENDPOINT_CIRCLE));
}

/**
 * @return the edge with endpoints reversed.
 */
Edge Edge::reverse() {
    Edge newEdge(node2, node1, endpoint1, endpoint2);
    return newEdge;
}

bool Edge::pointingLeft(Endpoint endpoint1, Endpoint endpoint2) {
    return (endpoint1 == ENDPOINT_ARROW && (endpoint2 == ENDPOINT_TAIL || endpoint2 == ENDPOINT_CIRCLE));
}



/* STATIC FUNCTIONS (From Edges.java) */

/**
 * Constructs a new bidirected edge from nodeA to nodeB (<->).
 */
Edge Edge::bidirectedEdge(const Node& nodeA, const Node& nodeB) {
    if (nodeA.getName() < nodeB.getName()) {
	return Edge(nodeA, nodeB, ENDPOINT_ARROW, ENDPOINT_ARROW);
    } 
    return Edge(nodeB, nodeA, ENDPOINT_ARROW, ENDPOINT_ARROW);
}

/**
 * Constructs a new directed edge from nodeA to nodeB (-->).
 */
Edge Edge::directedEdge(const Node& nodeA, const Node& nodeB) {
    Edge newEdge(nodeA, nodeB, ENDPOINT_TAIL, ENDPOINT_ARROW);
    return newEdge;
}

/**
 * Constructs a new partially oriented edge from nodeA to nodeB (o->).
 */
Edge Edge::partiallyOrientedEdge(const Node& nodeA, const Node& nodeB) {
    Edge newEdge(nodeA, nodeB, ENDPOINT_CIRCLE, ENDPOINT_ARROW);
    return newEdge;
}

/**
 * Constructs a new nondirected edge from nodeA to nodeB (o-o).
 */
Edge Edge::nondirectedEdge(const Node& nodeA, const Node& nodeB) {
    if (nodeA.getName() < nodeB.getName()) {
	return Edge(nodeA, nodeB, ENDPOINT_CIRCLE, ENDPOINT_CIRCLE);
    } 
    return Edge(nodeB, nodeA, ENDPOINT_CIRCLE, ENDPOINT_CIRCLE);
}

/**
 * Constructs a new undirected edge from nodeA to nodeB (--).
 */
Edge Edge::undirectedEdge(const Node& nodeA, const Node& nodeB) {
    if (nodeA.getName() < nodeB.getName()) {
	return Edge(nodeA, nodeB, ENDPOINT_TAIL, ENDPOINT_TAIL);
    } 
    return Edge(nodeB, nodeA, ENDPOINT_TAIL, ENDPOINT_TAIL);
}

/**
 * @return true iff an edge is a bidirected edge (<->).
 */
bool Edge::isBidirectionalEdge(const Edge& edge) {
    return (edge.endpoint1 == ENDPOINT_ARROW) && (edge.endpoint2 == ENDPOINT_ARROW);
}

/**
 * @return true iff an edge is a directed edge (-->).
 */
bool Edge::isDirectedEdge(const Edge& edge) {
    return ((edge.endpoint1 == ENDPOINT_TAIL) && (edge.endpoint2 == ENDPOINT_ARROW)) ||
            (edge.endpoint2 == ENDPOINT_TAIL) && (edge.endpoint1 == ENDPOINT_ARROW);
}

/**
 * @return true iff an edge is a partially oriented edge (o->)
 */
bool Edge::isPartiallyOrientedEdge(const Edge& edge) {
    return ((edge.endpoint1 == ENDPOINT_CIRCLE) && (edge.endpoint2 == ENDPOINT_ARROW)) ||
            (edge.endpoint2 == ENDPOINT_CIRCLE) && (edge.endpoint1 == ENDPOINT_ARROW);
}

/**
 * @return true iff an edge is a nondirected edge (o-o).
 */
bool Edge::isNondirectedEdge(const Edge& edge) {
    return (edge.endpoint1 == ENDPOINT_CIRCLE) && (edge.endpoint2 == ENDPOINT_CIRCLE);
}

/**
 * @return true iff an edge is a undirected edge (-).
 */
bool Edge::isUndirectedEdge(const Edge& edge) {
    return (edge.endpoint1 == ENDPOINT_TAIL) && (edge.endpoint2 == ENDPOINT_TAIL);
}

/**
 * @return the node opposite the given node along the given edge.
 * 
 * Return null if the node is not part of the edge
 */
Node Edge::traverse(const Node& node, Edge& edge) {
    if (node.isNull()) return Node();

    if (node == edge.node1) return edge.node2;
    if (node == edge.node2) return edge.node1;

    return Node();
}

/**
 * For A -> B, given A, returns B; otherwise returns null.
 */
Node Edge::traverseDirected(const Node& node, Edge& edge) {
    if (node == edge.node1 && edge.endpoint1 == ENDPOINT_TAIL && edge.endpoint2 == ENDPOINT_ARROW)
        return edge.node2;
    if (node == edge.node2 && edge.endpoint2 == ENDPOINT_TAIL && edge.endpoint1 == ENDPOINT_ARROW)
        return edge.node1;
    
    return Node();
}

/**
 * For A -> B, given B, returns A; otherwise returns null.
 */
Node Edge::traverseReverseDirected(const Node& node, Edge& edge) {
    if (node == edge.node1 && edge.endpoint1 == ENDPOINT_ARROW && edge.endpoint2 == ENDPOINT_TAIL)
        return edge.node2;
    if (node == edge.node2 && edge.endpoint2 == ENDPOINT_ARROW && edge.endpoint1 == ENDPOINT_TAIL)
        return edge.node1;
    
    return Node();
}

Node Edge::traverseReverseSemiDirected(const Node& node, Edge& edge) {
    if (node == edge.node1 && (edge.endpoint2 == ENDPOINT_TAIL || edge.endpoint2 == ENDPOINT_CIRCLE))
        return edge.node2;
    if (node == edge.node2 && (edge.endpoint1 == ENDPOINT_TAIL || edge.endpoint1 == ENDPOINT_CIRCLE))
        return edge.node1;
    
    return Node();
}

/**
 * For A --* B or A o-* B, given A, returns B. For A <-* B, returns null.
 */
Node Edge::traverseSemiDirected(const Node& node, Edge& edge) {
    if (node == edge.node1 && (edge.endpoint1 == ENDPOINT_TAIL || edge.endpoint1 == ENDPOINT_CIRCLE))
        return edge.node2;
    if (node == edge.node2 && (edge.endpoint2 == ENDPOINT_TAIL || edge.endpoint2 == ENDPOINT_CIRCLE))
        return edge.node1;
    
    return Node();
}

/**
 * For a directed edge, returns the node adjacent to the arrow endpoint.
 *
 * @throws std::invalid_argument if the given edge is not a directed
 *                                  edge.
 */
const Node& Edge::getDirectedEdgeHead(Edge& edge) {
    if (edge.endpoint1 == ENDPOINT_ARROW && edge.endpoint2 == ENDPOINT_TAIL)
        return edge.node1;
    if (edge.endpoint2 == ENDPOINT_ARROW && edge.endpoint1 == ENDPOINT_TAIL)
        return edge.node2;

    // Convert edge to string
    std::ostringstream ss;
    ss << edge;
    throw std::invalid_argument("Not a directed edge: " + ss.str());
}

/**
 * For a directed edge, returns the node adjacent to the null endpoint.
 *
 * @throws std::invalid_argument if the given edge is not a directed
 *                                  edge.
 */
const Node& Edge::getDirectedEdgeTail(Edge& edge) {
    if (edge.endpoint2 == ENDPOINT_ARROW && edge.endpoint1 == ENDPOINT_TAIL)
        return edge.node1;
    if (edge.endpoint1 == ENDPOINT_ARROW && edge.endpoint2 == ENDPOINT_TAIL)
        return edge.node2;

    // Convert edge to string
    std::ostringstream ss;
    ss << edge;
    throw std::invalid_argument("Not a directed edge: " + ss.str());
}

std::string Edge::toString() const {
    std::ostringstream result;

    result << *this;

    return result.str();
}

std::ostream& operator<<(std::ostream& os, const Edge& edge) {
    os << edge.node1.getName() << " ";

    switch(edge.endpoint1) {
        case ENDPOINT_TAIL:
            os << "-";
            break;
        case ENDPOINT_ARROW:
            os << "<";
            break;
        case ENDPOINT_CIRCLE:
            os << "o";
            break;
    }

    os << "-";

    switch(edge.endpoint2) {
        case ENDPOINT_TAIL:
            os << "-";
            break;
        case ENDPOINT_ARROW:
            os << ">";
            break;
        case ENDPOINT_CIRCLE:
            os << "o";
            break;
    }

    os << " " << edge.node2.getName();

    return os;
}

/**
 * Two edges are equal just in case they connect the SAME nodes (not a copy) and have the
 * same endpoints proximal to each node.
 */
bool operator==(const Edge& e1, const Edge& e2) {
    bool equal;

    if (e1.node1 == e2.node1 && e1.node2 == e2.node2) {
        equal = (e1.endpoint1 == e2.endpoint1 && e1.endpoint2 == e2.endpoint2);
    } else {
        equal = (e1.node1 == e2.node2 && e1.node2 == e2.node1 && e1.endpoint1 == e2.endpoint2 && e1.endpoint2 == e2.endpoint1);
    }

    return equal;
}

bool operator!=(const Edge& e1, const Edge& e2) {
    return !(e1 == e2);
}

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

void Edge::sortEdges(std::vector<Edge>& edges) {
    // for (auto it = edges.begin(); it != edges.end(); it++) {
    // 	if (isUndirectedEdge(*it)) {
    // 	    (*it) = undirectedEdge(it->node1, it->node2);
    // 	} else if (isNondirectedEdge(*it)) {
    // 	    (*it) = nondirectedEdge(it->node1, it->node2);
    // 	} else if (isBidirectionalEdge(*it)) {
    // 	    (*it) = bidirectedEdge(it->node1, it->node2);
    // 	}
    // }
    std::sort(edges.begin(), edges.end());
}
