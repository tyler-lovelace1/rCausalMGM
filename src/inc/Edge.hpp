#ifndef EDGE_HPP_
#define EDGE_HPP_

#include "Node.hpp"

enum Endpoint {ENDPOINT_TAIL, ENDPOINT_ARROW, ENDPOINT_CIRCLE, ENDPOINT_STAR, ENDPOINT_NULL};
enum EdgeProperty {dd, nl, pd, pl};

class Edge {

private:

    Node node1;
    Node node2;
    Endpoint endpoint1 = ENDPOINT_NULL;
    Endpoint endpoint2 = ENDPOINT_NULL;

    // /**
    //  * Default constructor node will serve as a null node
    //  */
    // static Node nullNode;

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
    Edge(const Node& node1, const Node& node2, Endpoint endpoint1, Endpoint endpoint2);

    Edge(const Edge& edge);

    // Used by both constructors
    void init(const Node& node1, const Node& node2, Endpoint endpoint1, Endpoint endpoint2);

    const Node& getNode1() const { return node1; }
    const Node& getNode2() const { return node2; }
    Endpoint getEndpoint1() const { return endpoint1; }
    Endpoint getEndpoint2() const { return endpoint2; }
    void setEndpoint1(Endpoint e) { endpoint1 = e; }
    void setEndpoint2(Endpoint e) { endpoint2 = e; }

    // bool isNull() { return node1.isNull() || node2.isNull(); }

    /**
     * @return the endpoint nearest to the given node.
     * @throws std::invalid_argument if the given node is not along the
     *                                  edge.
     */
    Endpoint getProximalEndpoint(const Node& node);
    // Endpoint getProximalEndpoint(Node&& node);

    /**
     * @return the endpoint furthest from the given node.
     * @throws std::invalid_argument if the given node is not along the
     *                                  edge.
     */
    Endpoint getDistalEndpoint(const Node& node);
    // Endpoint getDistalEndpoint(Node&& node);

    /**
     * Traverses the edge in an undirected fashion--given one node along the
     * edge, returns the node at the opposite end of the edge.
     * 
     * Returns null if the given node is not part of the edge 
     */
    Node getDistalNode(const Node& node);
    // const Node& getDistalNode(Node&& node);

    /**
     * @return true just in case this edge is directed.
     */
    bool isDirected();

    /**
     * @return true just in case this edge is undirected.
     */
    bool isUndirected();

    /**
     * @return true just in case this edge is bidirected.
     */
    bool isBidirected();

    /**
     * @return true just in case this edge is nondirected.
     */
    bool isNondirected();

    /**
     * @return true just in case this edge is partially oriented.
     */
    bool isPartiallyOriented();
    
    /**
     * @return true just in case the edge is pointing toward the given node--
     * that is, x --> node or x o--> node.
     */
    bool pointsTowards(const Node& node);

    /**
     * @return the edge with endpoints reversed.
     */
    Edge reverse();

    bool isNull() { return node1.isNull() || node2.isNull() || (endpoint1 == ENDPOINT_NULL && endpoint2 == ENDPOINT_NULL); }

    bool isDashed() { return dashed; }

    void addProperty(EdgeProperty property) { properties.insert(property); }
    void removeProperty(EdgeProperty property) { properties.erase(property); }
    std::unordered_set<EdgeProperty> getProperties() { return properties; }

    std::string toString() const;

    friend std::ostream& operator<<(std::ostream& os, const Edge& Edge);
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
    static Edge bidirectedEdge(const Node& nodeA, const Node& nodeB);

    /**
     * Constructs a new directed edge from nodeA to nodeB (-->).
     */
    static Edge directedEdge(const Node& nodeA, const Node& nodeB);

    /**
     * Constructs a new partially oriented edge from nodeA to nodeB (o->).
     */
    static Edge partiallyOrientedEdge(const Node& nodeA, const Node& nodeB);

    /**
     * Constructs a new nondirected edge from nodeA to nodeB (o-o).
     */
    static Edge nondirectedEdge(const Node& nodeA, const Node& nodeB);

    /**
     * Constructs a new undirected edge from nodeA to nodeB (--).
     */
    static Edge undirectedEdge(const Node& nodeA, const Node& nodeB);

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
    static Node traverse(const Node& node, Edge& edge);

    /**
     * For A -> B, given A, returns B; otherwise returns null.
     */
    static Node traverseDirected(const Node& node, Edge& edge);

    /**
     * For A -> B, given B, returns A; otherwise returns null.
     */
    static Node traverseReverseDirected(const Node& node, Edge& edge);

    static Node traverseReverseSemiDirected(const Node& node, Edge& edge);

    /**
     * For A --* B or A o-* B, given A, returns B. For A <-* B, returns null.
     */
    static Node traverseSemiDirected(const Node& node, Edge& edge);

    // Same as traverse()
    // static const Node& traverseUndirected(const Node& node, Edge& edge);

    /**
     * For a directed edge, returns the node adjacent to the arrow endpoint.
     *
     * @throws std::invalid_argument if the given edge is not a directed
     *                                  edge.
     */
    static const Node& getDirectedEdgeHead(Edge& edge);

    /**
     * For a directed edge, returns the node adjacent to the null endpoint.
     *
     * @throws std::invalid_argument if the given edge is not a directed
     *                                  edge.
     */
    static const Node& getDirectedEdgeTail(Edge& edge);

    static void sortEdges(std::vector<Edge>& edges);

};

template<> struct std::hash<Edge> {
public:
    std::size_t operator()(const Edge& k) const {
	std::size_t res = 7;
	res += 11 * (std::hash<Node>()(k.node1) + std::hash<Node>()(k.node2));
	return res;
	// return std::hash<std::string>()(k.toString());
    }
};

#endif /* EDGE_HPP_ */
