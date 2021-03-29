#ifndef EDGELISTGRAPH_HPP_
#define EDGELISTGRAPH_HPP_

#include "Variable.hpp"
#include "Edge.hpp"
#include "Triple.hpp"
#include "DataSet.hpp"


class EdgeListGraph {

private:
    /**
     * A list of the nodes in the graph, in the order in which they were added.
     */
    std::vector<Variable*> nodes;

    /**
     * The edges in the graph.
     */
    std::unordered_set<Edge> edgesSet;

    /**
     * Map from each node to the List of edges connected to that node.
     */
    std::unordered_map<Variable*, std::vector<Edge>> edgeLists;

    // TODO - PCS?
    // TODO - Triples?

    std::unordered_set<Triple> ambiguousTriples;
    std::unordered_set<Triple> underLineTriples;
    std::unordered_set<Triple> dottedUnderLineTriples;

    /**
     * Information about the algorithm used to build the graph
     */ 
    std::string algorithm;

    /**
     * undirected, true, or markov_equivalence_class
     */ 
    std::string graph_type;

    /**
     * True iff nodes were removed since the last call to an accessor for ambiguous, underline, or dotted underline
     * triples. If there are triples in the lists involving removed nodes, these need to be removed from the lists
     * first, so as not to cause confusion.
     */
    // bool stuffRemovedSinceLastTripleAccess = false;

    // std::unordered_set<Edge> highlightedEdges;

    std::unordered_map<std::string, Variable*> namesHash;

    // std::unordered_map<Variable*, std::unordered_set<Variable*>> ancestors;

    // void collectAncestorsVisit(Variable* node, std::unordered_set<Variable*>& ancestors);

    void collectDescendantsVisit(Variable* node, std::unordered_set<Variable*>& descendants);

    /**
     * closure under the child relation
     */
    void doChildClosureVisit(Variable* node, std::unordered_set<Variable*>& closure);

    /**
     * This is a simple auxiliary visit method for the isDConnectedTo() method
     * used to find the closure of a conditioning set of nodes under the parent
     * relation.
     *
     * @param node    the node in question
     * @param closure the closure of the conditioning set uner the parent
     *                relation (to be calculated recursively).
     */
    void doParentClosureVisit(Variable* node, std::unordered_set<Variable*>& closure);

    /**
     * @return true iff there is a directed path from node1 to node2.
     */
    bool existsUndirectedPathVisit(Variable* node1, Variable* node2, std::unordered_set<Variable*>& path);

    bool existsDirectedPathVisit(Variable* node1, Variable* node2, std::unordered_set<Variable*>& path);

    /**
     * @return true iff there is a semi-directed path from node1 to node2
     */
    bool existsSemiDirectedPathVisit(Variable* node1, std::unordered_set<Variable*>& nodes2, std::vector<Variable*>& path);

public:

    /**
     * Constructs a new (empty) EdgeListGraph.
     */
    EdgeListGraph();

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
    EdgeListGraph(const EdgeListGraph& graph);

    /**
     * Constructs a new graph, with no edges, using the the given variable
     * names.
     */
    EdgeListGraph(const std::vector<Variable*>& nodes);

    /**
     * Makes a graph from an R list
     */ 
    EdgeListGraph(const Rcpp::List& list, DataSet& ds);

    // Shallow copy isn't possible because Edges aren't stored by reference (is it neccesary?)

    // Used by constructors
    void initNamesHash();

    /**
     * Adds a directed edge to the graph from node A to node B.
     *
     * @param node1 the "from" node.
     * @param node2 the "to" node.
     */
    bool addDirectedEdge(Variable* node1, Variable* node2);

    /**
     * Adds an undirected edge to the graph from node A to node B.
     *
     * @param node1 the "from" node.
     * @param node2 the "to" node.
     */
    bool addUndirectedEdge(Variable* node1, Variable* node2);

    /**
     * Adds a partially oriented edge to the graph from node A to node B.
     *
     * @param node1 the "from" node.
     * @param node2 the "to" node.
     */
    bool addPartiallyOrientedEdge(Variable* node1, Variable* node2);

    /**
     * Adds a bidirected edge to the graph from node A to node B.
     *
     * @param node1 the "from" node.
     * @param node2 the "to" node.
     */
    bool addBidirectedEdge(Variable* node1, Variable* node2);

    bool existsDirectedCycle();

    bool isDirectedFromTo(Variable* node1, Variable* node2);
    bool isUndirectedFromTo(Variable* node1, Variable* node2);

    /**
     * added by ekorber, 2004/06/11
     *
     * @return true if the given edge is definitely visible (Jiji, pg 25)
     * @throws std::invalid_argument if the given edge is not a directed edge
     *                                  in the graph
     */
    bool defVisible(Edge& edge);

    /**
     * std::invalid_argument  raised (by isDirectedFromTo(getEndpoint) or by
     * getEdge) if there are multiple edges between any of the node pairs.
     */
    bool isDefNoncollider(Variable* node1, Variable* node2, Variable* node3);

    bool isDefCollider(Variable* node1, Variable* node2, Variable* node3);

    /**
     * @return true iff there is a directed path from node1 to node2.
     * a
     */
    bool existsDirectedPathFromTo(Variable* node1, Variable* node2);
    bool existsUndirectedPathFromTo(Variable* node1, Variable* node2);
    bool existsSemiDirectedPathFromTo(Variable* node1, std::unordered_set<Variable*>& nodes);

    /**
     * Determines whether a trek exists between two nodes in the graph.  A trek
     * exists if there is a directed path between the two nodes or else, for
     * some third node in the graph, there is a path to each of the two nodes in
     * question.
     */
    bool existsTrek(Variable* node1, Variable* node2);

    /**
     * @return the list of children for a node.
     */
    std::vector<Variable*> getChildren(Variable* node);

    int getConnectivity();

    std::vector<Variable*> getDescendants(std::vector<Variable*>& nodes);

    /**
     * @return the edge connecting node1 and node2, provided a unique such edge
     * exists.
     * 
     * Throws std::invalid_argument if not
     */
    Edge getEdge(Variable* node1, Variable* node2);

    Edge getDirectedEdge(Variable* node1, Variable* node2);

    /**
     * @return the list of parents for a node.
     */
    std::vector<Variable*> getParents(Variable* node);

    /**
     * @return the number of edges into the given node.
     */
    int getIndegree(Variable* node);

    int getDegree(Variable* node);

    /**
     * @return the number of edges out of the given node.
     */
    int getOutdegree(Variable* node);

    /**
     * Determines whether some edge or other exists between two nodes.
     */
    bool isAdjacentTo(Variable* node1, Variable* node2);

    /**
     * Determines whether one node is an ancestor of another.
     */
    bool isAncestorOf(Variable* node1, Variable* node2);

    bool possibleAncestor(Variable* node1, Variable* node2);

    /**
     * @return true iff node1 is a possible ancestor of at least one member of
     * nodes2
     */
    // bool possibleAncestorSet(Variable* node1, std::vector<Variable*>& nodes2);

    // std::vector<Variable*> getAncestors(std::vector<Variable*>& nodes);

    /**
     * Determines whether one node is a child of another.
     */
    bool isChildOf(Variable* node1, Variable* node2);

    /**
     * Determines whether one node is a descendent of another.
     */
    bool isDescendentOf(Variable* node1, Variable* node2);

    /**
     * @return true iff node2 is a definite nondecendent of node1
     */
    bool defNonDescendent(Variable* node1, Variable* node2);

    bool isDConnectedTo(Variable* x, Variable* y, std::vector<Variable*>& z);

    bool isDConnectedTo(std::vector<Variable*>& x, std::vector<Variable*>& y, std::vector<Variable*>& z);

    std::vector<Variable*> getSepset(Variable* x, Variable* y);

    void setNodes(std::vector<Variable*>& nodes);

    // std::unordered_set<Variable*> zAncestors(std::vector<Variable*>& z);

    bool isDSeparatedFrom(std::vector<Variable*>& x, std::vector<Variable*>& y, std::vector<Variable*>& z);

    /**
     * Determines whether one n ode is d-separated from another. According to
     * Spirtes, Richardson & Meek, two nodes are d- connected given some
     * conditioning set Z if there is an acyclic undirected path U between them,
     * such that every collider on U is an ancestor of some element in Z and
     * every non-collider on U is not in Z.  Two elements are d-separated just
     * in case they are not d-connected.  A collider is a node which two edges
     * hold in common for which the endpoints leading into the node are both
     * arrow endpoints.
     *
     * @param node1 the first node.
     * @param node2 the second node.
     * @param z     the conditioning set.
     * @return true if node1 is d-separated from node2 given set t, false if
     * not.
     * @see #isDConnectedTo
     */
    bool isDSeparatedFrom(Variable* node1, Variable* node2, std::vector<Variable*>& z);

    bool possDConnectedTo(Variable* node1, Variable* node2, std::vector<Variable*>& condNodes);

    /**
     * Determines whether an inducing path exists between node1 and node2, given
     * a set O of observed nodes and a set sem of conditioned nodes.
     *
     * @param node1 the first node.
     * @param node2 the second node.
     * @return true if an inducing path exists, false if not.
     */
    bool existsInducingPath(Variable* node1, Variable* node2);

    /**
     * Determines whether one node is a parent of another.
     *
     * @param node1 the first node.
     * @param node2 the second node.
     * @return true if node1 is a parent of node2, false if not.
     * @see #isChildOf
     * @see #getParents
     * @see #getChildren
     */
    bool isParentOf(Variable* node1, Variable* node2);

    /**
     * Determines whether one node is a proper ancestor of another.
     */
    bool isProperAncestorOf(Variable* node1, Variable* node2);

    /**
     * Determines whether one node is a proper decendent of another
     */
    bool isProperDescendentOf(Variable* node1, Variable* node2);

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
    void transferNodesAndEdges(const EdgeListGraph& graph);

    /**
     * Determines whether a node in a graph is exogenous.
     */
    bool isExogenous(Variable* node);

    /**
     * @return the set of nodes adjacent to the given node. If there are multiple edges between X and Y, Y will show
     * up twice in the list of adjacencies for X, for optimality; simply create a list an and array from these to
     * eliminate the duplication.
     */
    std::vector<Variable*> getAdjacentNodes(Variable* node);

    /**
     * Removes the edge connecting the two given nodes.
     */
    bool removeEdge(Variable* node1, Variable* node2);

    /**
     * @return the endpoint along the edge from node to node2 at the node2 end.
     */
    Endpoint getEndpoint(Variable* node1, Variable* node2);

    /**
     * If there is currently an edge from node1 to node2, sets the endpoint at
     * node2 to the given endpoint; if there is no such edge, adds an edge --#
     * where # is the given endpoint. Setting an endpoint to null, provided
     * there is exactly one edge connecting the given nodes, removes the edge.
     * (If there is more than one edge, an exception is thrown.)
     *
     * @throws IllegalArgumentException if the edge with the revised endpoint
     *                                  cannot be added to the graph.
     */
    bool setEndpoint(Variable* from, Variable* to, Endpoint endPoint);

    /**
     * Nodes adjacent to the given node with the given proximal endpoint.
     */
    std::vector<Variable*> getNodesInTo(Variable* node, Endpoint endpoint);

    /**
     * Nodes adjacent to the given node with the given distal endpoint.
     */
    std::vector<Variable*> getNodesOutTo(Variable* node, Endpoint endpoint);
    
    void collectAncestorsVisit(Variable* node, std::unordered_set<Variable*> ancestors);

    std::unordered_set<Variable*> getAncestors(std::vector<Variable*>& nodes);

    /**
     * @return a matrix of endpoints for the nodes in this graph, with nodes in
     * the same order as getNodes().
     */
    arma::Mat<Endpoint> getEndpointMatrix();

    /**
     * Adds an edge to the graph.
     *
     * @param edge the edge to be added
     * @return true if the edge was added, false if not.
     */
    bool addEdge(Edge& edge);

    /**
     * Add edge from string
     */ 
    bool addEdge(std::string edgeString);

    /**
     * Adds a node to the graph. Precondition: The proposed name of the node
     * cannot already be used by any other node in the same graph.
     *
     * @param node the node to be added.
     * @return true if the the node was added, false if not.
     */
    bool addNode(Variable* node);

    /**
     * @return the set of edges in the graph.
     */
    std::unordered_set<Edge> getEdges() { return edgesSet; }

    /**
     * @return the list of edges in the graph.  No particular ordering of the
     * edges in the list is guaranteed.
     */
    std::vector<Edge> getEdgeList() { return std::vector<Edge>(edgesSet.begin(), edgesSet.end()); }

    std::unordered_set<Triple> getAmbiguousTriples() { return ambiguousTriples; }
    std::unordered_set<Triple> getUnderLines() { return underLineTriples; }
    std::unordered_set<Triple> getDottedUnderlines() { return dottedUnderLineTriples; }

    bool isAmbiguousTriple(Variable* x, Variable* y, Variable* z) { return ambiguousTriples.count(Triple(x, y, z)); }
    bool isUnderlineTriple(Variable* x, Variable* y, Variable* z) { return underLineTriples.count(Triple(x, y, z)); }
    bool isDottedUnderlineTriple(Variable* x, Variable* y, Variable* z) { return dottedUnderLineTriples.count(Triple(x, y, z)); }

    void addAmbiguousTriple(Variable* x, Variable* y, Variable* z) { ambiguousTriples.insert(Triple(x, y, z)); } 
    void addUnderlineTriple(Variable* x, Variable* y, Variable* z) { underLineTriples.insert(Triple(x, y, z)); } 
    void addDottedUnderlineTriple(Variable* x, Variable* y, Variable* z) { dottedUnderLineTriples.insert(Triple(x, y, z)); }

    void removeAmbiguousTriple(Variable* x, Variable* y, Variable* z);
    void removeUnderlineTriple(Variable* x, Variable* y, Variable* z);
    void removeDottedUnderlineTriple(Variable* x, Variable* y, Variable* z);

    void setAmbiguousTriples(std::unordered_set<Triple>& triples);
    void setUnderLineTriples(std::unordered_set<Triple>& triples);
    void setDottedUnderLineTriples(std::unordered_set<Triple>& triples);

    void removeTriplesNotInGraph();

    std::vector<std::vector<Triple>> getTriplesLists(Variable* node);

    Triple tripleFromString(std::string tripleString);

    /**
     * Determines if the graph contains a particular edge.
     */
    bool containsEdge(Edge& edge) { return edgesSet.count(edge); }



    /**
     * @return the list of edges connected to a particular node. No particular
     * ordering of the edges in the list is guaranteed.
     */
    std::vector<Edge> getEdges(Variable* node) { 
        // Rcpp::Rcout << "Getting edges..." << std::endl;
        return edgeLists[node]; 
    }

    // TODO - hash code?

    /**
     * Resets the graph so that it is fully connects it using #-# edges, where #
     * is the given endpoint.
     */
    void fullyConnect(Endpoint endpoint);

    void reorientAllWith(Endpoint endpoint);

    /**
     * @return the node with the given name, or null if no such node exists.
     */
    Variable* getNode(std::string name) { auto itr = namesHash.find(name); return itr == namesHash.end() ? NULL : itr->second; }

    /**
     * @return the number of nodes in the graph.
     */
    int getNumNodes() { return nodes.size(); }

    /**
     * @return the number of edges in the (entire) graph.
     */
    int getNumEdges() { return edgesSet.size(); }

    /**
     * @return the number of edges connected to a particular node in the graph.
     */
    int getNumEdges(Variable* node);

    std::vector<Variable*> getNodes() { return nodes; }

    /**
     * Removes all nodes (and therefore all edges) from the graph.
     */
    void clear();

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
    bool removeEdge(Edge& edge);

    /**
     * Removes any relevant edge objects found in this collection. G
     *
     * @param edges the collection of edges to remove.
     * @return true if any edges in the collection were removed, false if not.
     */
    bool removeEdges(const std::vector<Edge>& edges);

    /**
     * Removes all edges connecting node A to node B.
     *
     * @param node1 the first node.,
     * @param node2 the second node.
     * @return true if edges were removed between A and B, false if not.
     */
    bool removeEdges(Variable* node1, Variable* node2);

    /**
     * Removes all edges
     */
    void removeEdges();

    /**
     * Removes a node from the graph.
     */
    bool removeNode(Variable* node);

    /**
     * Removes any relevant node objects found in this collection.
     *
     * @param newNodes the collection of nodes to remove.
     * @return true if nodes from the collection were removed, false if not.
     */
    bool removeNodes(std::vector<Variable*>& newNodes);

    EdgeListGraph subgraph(std::vector<Variable*>& nodes);

    /**
     * @return the edges connecting node1 and node2.
     */
    std::vector<Edge> getEdges(Variable* node1, Variable* node2);

    // TODO - Triples?

    std::vector<std::string> getNodeNames();

    std::vector<Variable*> getCausalOrdering();

    // void setHighlighted(Edge& edge, bool highlighted) { if (highlighted) highlightedEdges.insert(edge); else highlightedEdges.erase(edge); } ;

    // bool isHighlighted(Edge& edge) { return highlightedEdges.count(edge); }

    void changeName(std::string name, std::string newName);

    // Converts graph to a list in R
    Rcpp::List toList();

    // Returns true if an R list object is a valid graph
    static bool validateGraphList(const Rcpp::List& l);

    void setAlgorithm(std::string a) { algorithm = a; }
    std::string getAlgorithm() { return algorithm; } 

    void setGraphType(std::string t) { graph_type = t; }
    std::string getGraphType() { return graph_type; } 

    /**
     * @return true iff the given object is a graph that is equal to this graph,
     * in the sense that it contains the same nodes and the edges are
     * isomorphic.
     */
    friend bool operator==(const EdgeListGraph& g1, const EdgeListGraph& g2);
    friend std::ostream& operator<<(std::ostream& os, EdgeListGraph& graph);


};

#endif /* EDGELISTGRAPH_HPP_ */
