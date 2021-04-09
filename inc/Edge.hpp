#ifndef EDGE_HPP_
#define EDGE_HPP_

#include "Variable.hpp"

enum Endpoint {ENDPOINT_TAIL, ENDPOINT_ARROW, ENDPOINT_CIRCLE, ENDPOINT_STAR, ENDPOINT_NULL};
enum EdgeProperty {dd, nl, pd, pl};

class Edge {

private:

    Variable* node1;
    Variable* node2;
    Endpoint endpoint1;
    Endpoint endpoint2;

    //TODO - color?

    bool dashed = false;

    std::unordered_set<EdgeProperty> properties;

    bool pointingLeft(Endpoint endpoint1, Endpoint endpoint2);

public:
    Edge() {}

    /**
     * Constructs a new edge by specifying the nodes it connects and the
     * endpoint types.
     *
     * @param node1     the first node
     * @param node2     the second node            _
     * @param endpoint1 the endpoint at the first node
     * @param endpoint2 the endpoint at the second node
     */
    Edge(Variable* node1, Variable* node2, Endpoint endpoint1, Endpoint endpoint2);

    Edge(const Edge& edge);

    // Used by both constructors
    void init(Variable* node1, Variable* node2, Endpoint endpoint1, Endpoint endpoint2);

    Variable* getNode1() { return node1; }
    Variable* getNode2() { return node2; }
    Endpoint getEndpoint1() { return endpoint1; }
    Endpoint getEndpoint2() { return endpoint2; }
    void setEndpoint1(Endpoint e) { endpoint1 = e; }
    void setEndpoint2(Endpoint e) { endpoint2 = e; }

    /**
     * @return the endpoint nearest to the given node.
     * @throws std::invalid_argument if the given node is not along the
     *                                  edge.
     */
    Endpoint getProximalEndpoint(Variable* node);

    /**
     * @return the endpoint furthest from the given node.
     * @throws std::invalid_argument if the given node is not along the
     *                                  edge.
     */
    Endpoint getDistalEndpoint(Variable* node);

    /**
     * Traverses the edge in an undirected fashion--given one node along the
     * edge, returns the node at the opposite end of the edge.
     * 
     * Returns null if the given node is not part of the edge 
     */
    Variable* getDistalNode(Variable* node);

    /**
     * @return true just in case this edge is directed.
     */
    bool isDirected();

    /**
     * @return true just in case the edge is pointing toward the given node--
     * that is, x --> node or x o--> node.
     */
    bool pointsTowards(Variable* node);

    /**
     * @return the edge with endpoints reversed.
     */
    Edge reverse();

    bool isNull() { return endpoint1 == ENDPOINT_NULL && endpoint2 == ENDPOINT_NULL; }

    bool isDashed() { return dashed; }

    void addProperty(EdgeProperty property) { properties.insert(property); }
    void removeProperty(EdgeProperty property) { properties.erase(property); }
    std::unordered_set<EdgeProperty> getProperties() { return properties; }

    std::string toString();

    friend std::ostream& operator<<(std::ostream& os, Edge& Edge);
    friend bool operator==(const Edge& e1, const Edge& e2);
    friend bool operator!=(const Edge& e1, const Edge& e2);
    friend bool operator> (const Edge& e1, const Edge& e2);
    friend bool operator<= (const Edge& e1, const Edge& e2);
    friend bool operator< (const Edge& e1, const Edge& e2);
    friend bool operator>= (const Edge& e1, const Edge& e2);
    friend struct std::hash<Edge>;

    /* STATIC FUNCTIONS (From Edges.java) */

    /**
     * Constructs a new bidirected edge from nodeA to nodeB (<->).
     */
    static Edge bidirectedEdge(Variable* nodeA, Variable* nodeB);

    /**
     * Constructs a new directed edge from nodeA to nodeB (-->).
     */
    static Edge directedEdge(Variable* nodeA, Variable* nodeB);

    /**
     * Constructs a new partially oriented edge from nodeA to nodeB (o->).
     */
    static Edge partiallyOrientedEdge(Variable* nodeA, Variable* nodeB);

    /**
     * Constructs a new nondirected edge from nodeA to nodeB (o-o).
     */
    static Edge nondirectedEdge(Variable* nodeA, Variable* nodeB);

    /**
     * Constructs a new undirected edge from nodeA to nodeB (--).
     */
    static Edge undirectedEdge(Variable* nodeA, Variable* nodeB);

    /**
     * @return true iff an edge is a bidirected edge (<->).
     */
    static bool isBidirectionalEdge(const Edge& edge);

    /**
     * @return true iff an edge is a directed edge (-->).
     */
    static bool isDirectedEdge(const Edge& edge);

    /**
     * @return true iff an edge is a partially oriented edge (o->)
     */
    static bool isPartiallyOrientedEdge(const Edge& edge);

    /**
     * @return true iff an edge is a nondirected edge (o-o).
     */
    static bool isNondirectedEdge(const Edge& edge);

    /**
     * @return true iff an edge is a undirected edge (-).
     */
    static bool isUndirectedEdge(const Edge& edge);

    /**
     * @return the node opposite the given node along the given edge.
     */
    static Variable* traverse(Variable* node, Edge& edge);

    /**
     * For A -> B, given A, returns B; otherwise returns null.
     */
    static Variable* traverseDirected(Variable* node, Edge& edge);

    /**
     * For A -> B, given B, returns A; otherwise returns null.
     */
    static Variable* traverseReverseDirected(Variable* node, Edge& edge);

    static Variable* traverseReverseSemiDirected(Variable* node, Edge& edge);

    /**
     * For A --* B or A o-* B, given A, returns B. For A <-* B, returns null.
     */
    static Variable* traverseSemiDirected(Variable* node, Edge& edge);

    // Same as traverse()
    // static Variable* traverseUndirected(Variable* node, Edge& edge);

    /**
     * For a directed edge, returns the node adjacent to the arrow endpoint.
     *
     * @throws std::invalid_argument if the given edge is not a directed
     *                                  edge.
     */
    static Variable* getDirectedEdgeHead(Edge& edge);

    /**
     * For a directed edge, returns the node adjacent to the null endpoint.
     *
     * @throws std::invalid_argument if the given edge is not a directed
     *                                  edge.
     */
    static Variable* getDirectedEdgeTail(Edge& edge);

    static void sortEdges(std::vector<Edge>& edges);

};

template<> struct std::hash<Edge> {
public:
    std::size_t operator()(const Edge& k) const {
        return std::hash<Variable*>()(k.node1) + std::hash<Variable*>()(k.node2);
    }
};

#endif /* EDGE_HPP_ */
