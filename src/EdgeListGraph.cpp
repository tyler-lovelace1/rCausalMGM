#include "EdgeListGraph.hpp"

#include "GraphUtils.hpp"
#include <RcppArmadillo.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <cctype>

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

    // for (Edge edge : graph.edgesSet) {
    //     if (graph.highlightedEdges.count(edge)) {
    //         setHighlighted(edge, true);
    //     }
    // }

    namesHash = graph.namesHash;
    algorithm = graph.algorithm;
    graph_type = graph.graph_type;
}

EdgeListGraph::EdgeListGraph(const Rcpp::List& list, DataSet& ds)  {
    if (!validateGraphList(list)) {
        throw std::invalid_argument("ERROR: list is not in the form of a graph");
    }

    // Nodes
    std::vector<std::string> nodeNames = list["nodes"];
    for (std::string nodeName : nodeNames) {
        try {
            if(!addNode(ds.getVariable(nodeName)))
                throw std::invalid_argument("Issue adding variable to graph");
        } catch (const std::exception& ex) {
            throw std::invalid_argument("ERROR: Could not find node " + nodeName + " in the provided data set");
        }
    }
    initNamesHash();

    // Edges
    std::vector<std::string> edgeStrings = list["edges"];
    for (std::string edgeString : edgeStrings) {
        if (!addEdge(edgeString))
            throw std::invalid_argument("Error parsing edge: " + edgeString);
    }

    //Triples
    std::vector<std::string> tripleStrings = list["ambiguous_triples"];
    for (std::string tripleString : tripleStrings) {
        ambiguousTriples.insert(tripleFromString(tripleString));
    }

    std::vector<std::string> a = list["algorithm"];
    algorithm = a[0];

    std::vector<std::string> t = list["type"];
    graph_type = t[0];
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
 * Adds a directed edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addDirectedEdge(Variable* node1, Variable* node2) {
    Edge newEdge = Edge::directedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds a bidirected edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addBidirectedEdge(Variable* node1, Variable* node2) {
    Edge newEdge = Edge::bidirectedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds a partially oriented edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addPartiallyOrientedEdge(Variable* node1, Variable* node2) {
    Edge newEdge = Edge::partiallyOrientedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds a nondirected edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addNondirectedEdge(Variable* node1, Variable* node2) {
    Edge newEdge = Edge::nondirectedEdge(node1, node2);
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

bool EdgeListGraph::addEdge(std::string edgeString) {
    std::vector<std::string> edgeComponents = GraphUtils::splitString(edgeString, " ");

    if (edgeComponents.size() != 3)
        throw std::invalid_argument("Edge from string must have 3 components (node edge node): " + edgeString);

    Variable* node1 = getNode(edgeComponents[0]);
    Variable* node2 = getNode(edgeComponents[2]);

    if (node1 == NULL)
        throw std::invalid_argument("Edge node not found in graph: " + edgeComponents[0]);

    if (node2 == NULL)
        throw std::invalid_argument("Edge node not found in graph: " + edgeComponents[2]);

    std::string edgeMid = edgeComponents[1];

    if (edgeMid.length() < 3)
        throw std::invalid_argument("Invalid edge: " + edgeString);

    char endpoint1 = edgeMid[0];
    char endpoint2 = edgeMid[edgeMid.length()-1];

    if (endpoint1 == '-') {
        if      (endpoint2 == '>') return addDirectedEdge(node1, node2);
        else if (endpoint2 == '-') return addUndirectedEdge(node1, node2);

    } else if (endpoint1 == '<') {
        if      (endpoint2 == '>') return addBidirectedEdge(node1, node2);
        else if (endpoint2 == '-') return addDirectedEdge(node2, node1);

    } else if (endpoint1 == 'o') {
        if      (endpoint2 == 'o') return addNondirectedEdge(node1, node2);
        else if (endpoint2 == '>') return addPartiallyOrientedEdge(node1, node2);
    }

    throw std::invalid_argument("Endpoints not recognized: " + edgeString);
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
 * Removes all edges
 */
void EdgeListGraph::removeEdges() {
    std::vector<Edge> edges(edgesSet.begin(), edgesSet.end());
    removeEdges(edges);
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

    // highlightedEdges.erase(edge);
    // stuffRemovedSinceLastTripleAccess = true;

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

    throw std::invalid_argument("node1 and node2 not connected by edge. node1: " + node1->getName() + " node2: " + node2->getName());

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
    if (Edge::isUndirectedEdge(newEdge)) {
	newEdge = Edge::undirectedEdge(from, to);
    } else if (Edge::isNondirectedEdge(newEdge)) {
	newEdge = Edge::nondirectedEdge(from, to);
    } else if (Edge::isBidirectionalEdge(newEdge)) {
	newEdge = Edge::bidirectedEdge(from, to);
    }

    try {
        removeEdges(edge.getNode1(), edge.getNode2());
        addEdge(newEdge);
        return true;
    } catch (std::invalid_argument& e) {
        return false;
    }

    return false; // Unreachable
}

Triple EdgeListGraph::tripleFromString(std::string tripleString) {
    tripleString = tripleString.substr(1, tripleString.size()-2); // Strip < and >

    std::vector<std::string> nodeNames = GraphUtils::splitString(tripleString, ",");

    if (nodeNames.size() != 3)
        throw std::invalid_argument("Triple must take form <X,Y,Z>: <" + tripleString + ">");

    return Triple(
        getNode(nodeNames[0]),
        getNode(nodeNames[1]),
        getNode(nodeNames[2])
    );
}

bool EdgeListGraph::validateGraphList(const Rcpp::List& l) {
    std::vector<std::string> lclass = l.attr("class");
    if (lclass[0] != "graph") return false;

    std::vector<std::string> names = l.names();

    if (names.size() != 7)               return false;
    if (names[0] != "nodes")             return false;
    if (names[1] != "edges")             return false;
    if (names[2] != "ambiguous_triples") return false;
    if (names[3] != "algorithm")         return false;
    if (names[4] != "type")              return false;
    if (names[5] != "markov.blankets")   return false;
    if (names[6] != "stabilities")       return false;

    return true;
}

/**
 * Calculate markov blankets for undirected graph
 * For a node x, all neighbors of x are in the Markov Blanket of x
 */
Rcpp::List markovBlanketUndirected(const Rcpp::List& graph) {
    std::vector<std::string> nodes = graph["nodes"];
    std::vector<std::string> edges = graph["edges"];

    std::unordered_map<std::string, std::unordered_set<std::string>> blankets;

    for (std::string n : nodes) {
        blankets[n] = std::unordered_set<std::string>();
    }

    // Get blankets
    for (std::string edgeString : edges) {
        std::vector<std::string> e = GraphUtils::splitString(edgeString, " ");

        std::string n1 = e[0];
        std::string n2 = e[2];

        blankets[n1].insert(n2);
        blankets[n2].insert(n1);
    }

    // Convert to list
    Rcpp::List result = Rcpp::List::create();

    for (std::string n : nodes) {
        result[n] = std::vector<std::string>(blankets[n].begin(), blankets[n].end());
    }

    return result;
}

/**
 * Calculate markov blankets for Partial Ancestral graphs
 *
 * Rules:
 * 1. Parents, children, and spouses (linked by fully directed edges) are treated the 
 *    same as in DAGs, and are all included in the Markov blanket.
 * 2. Unoriented, partially oriented, and bidirected edges all have to be treated like 
 *    there is latent confounding, because it hasn't been ruled out. Thus, any node 
 *    connected to the target by unoriented, partially oriented, or 
 *    bidirected edges are included in the Markov blanket.
 * 3. Unoriented, partially oriented, or bidirected edges connected to the children are 
 *    included in the Markov blanket.
 * 4. For the nodes added to the Markov blanket in rule two, add any parents of those
 *    nodes to the Markov blanket.
 * 5. Additionally, add any nodes connected to those added in rule two by unoriented, 
 *    partially oriented, or bidirected edges. We're going to cut off the Markov blanket 
 *    here arbitrarily, but in theory rules 3 and 4 should be applied recursively on 
 *    each new set of nodes connected by unoriented, partially oriented, or bidirected 
 *    edges. However, this could lead to ridiculously large Markov blankets, and it would 
 *    be unlikely the additional nodes would actually improve predictive performance.
 */
Rcpp::List markovBlanketPAG(const Rcpp::List& graph) {
    std::vector<std::string> nodes = graph["nodes"];
    std::vector<std::string> edges = graph["edges"];

    std::unordered_map<std::string, std::unordered_set<std::string>> blankets;
    std::unordered_map<std::string, std::unordered_set<std::string>> confoundingNeighbors; // Nodes connected by o-o, o->, or <->
    std::unordered_map<std::string, std::unordered_set<std::string>> parents;
    std::unordered_map<std::string, std::unordered_set<std::string>> children;

    for (std::string n : nodes) {
        blankets[n] = std::unordered_set<std::string>();
        confoundingNeighbors[n] = std::unordered_set<std::string>();
        parents[n] = std::unordered_set<std::string>();
        children[n] = std::unordered_set<std::string>();
    }

    // Get neighbors of every node
    for (std::string edgeString : edges) {
        std::vector<std::string> e = GraphUtils::splitString(edgeString, " ");

        std::string n1 = e[0];
        std::string n2 = e[2];
        std::string edge = e[1];

        if (edge == "o-o" || edge == "o->" || edge == "<->") {
            confoundingNeighbors[n1].insert(n2);
            confoundingNeighbors[n2].insert(n1);
        } else if (edge == "-->") {
            children[n1].insert(n2);
            parents[n2].insert(n1);
        } 

    }

    // Get blankets of every node
    for (std::string target : nodes) {

        std::unordered_set<std::string> spouses;
        for (std::string child : children[target]) {
            spouses.insert(parents[child].begin(), parents[child].end());
        }

        std::unordered_set<std::string> rule2(confoundingNeighbors[target]);
	std::unordered_set<std::string> rule3;
        for (std::string child : children[target]) {
            rule3.insert(confoundingNeighbors[child].begin(), confoundingNeighbors[child].end());
        }

        std::unordered_set<std::string> rule4;
        std::unordered_set<std::string> rule5;
        for (std::string Y : rule2) {
            rule4.insert(parents[Y].begin(), parents[Y].end());
            rule5.insert(confoundingNeighbors[Y].begin(), confoundingNeighbors[Y].end());
        }

        blankets[target].insert(parents[target].begin(), parents[target].end());
        blankets[target].insert(children[target].begin(), children[target].end());
        blankets[target].insert(spouses.begin(), spouses.end());
        blankets[target].insert(rule2.begin(), rule2.end());
        blankets[target].insert(rule3.begin(), rule3.end());
        blankets[target].insert(rule4.begin(), rule4.end());
	blankets[target].insert(rule5.begin(), rule5.end());
        blankets[target].erase(target);
    }

    // Convert to list
    Rcpp::List result = Rcpp::List::create();

    for (std::string n : nodes) {
        result[n] = std::vector<std::string>(blankets[n].begin(), blankets[n].end());
    }

    return result;
}


/**
 * Calculate markov blankets for Markov Equivalence Class graphs
 *
 * Updated rules for how to handle Markov blankets in PDAGs (Partially Directed Acyclic Graphs):
 * 1. Parents of target node
 * 2. Children of target node
 * 3. Spouses of target node
 * 4. If the target variable X has an undirected edge to Y, then Y and its directed parents
 *    are in the Markov blanket.
 * 5. If the target variable X has an undirected edge to Y, we do not include any node Z
 *    connected to Y by an undirected edge. This can be explained because in the Markov
 *    equivalence class, two consecutive undirected edges X --- Y --- Z can be either
 *    X --> Y --> Z, X <-- Y <-- Z, or X <-- Y --> Z. In none of these cases is Z a parent of Y.
 * 6. In the case where the target variable X is a parent of node Y that contains an
 *    undirected edge to a node Z, then Y, Z, and any directed parents W of Z are included
 *    in the Markov blanket.
 *    #### Rule 6 commented out, TODO: test whether inclusion improves MB inference
 *
 */
Rcpp::List markovBlanketMEC(const Rcpp::List& graph) {
    std::vector<std::string> nodes = graph["nodes"];
    std::vector<std::string> edges = graph["edges"];

    std::unordered_map<std::string, std::unordered_set<std::string>> blankets;
    std::unordered_map<std::string, std::unordered_set<std::string>> undirectedNeighbors;
    std::unordered_map<std::string, std::unordered_set<std::string>> parents;
    std::unordered_map<std::string, std::unordered_set<std::string>> children;

    for (std::string n : nodes) {
        blankets[n] = std::unordered_set<std::string>();
        undirectedNeighbors[n] = std::unordered_set<std::string>();
        parents[n] = std::unordered_set<std::string>();
        children[n] = std::unordered_set<std::string>();
    }

    // Get neighbors of every node
    for (std::string edgeString : edges) {
        std::vector<std::string> e = GraphUtils::splitString(edgeString, " ");

        std::string n1 = e[0];
        std::string n2 = e[2];
        std::string edge = e[1];

        if (edge == "---") {
            undirectedNeighbors[n1].insert(n2);
            undirectedNeighbors[n2].insert(n1);
        } else if (edge == "-->") {
            children[n1].insert(n2);
            parents[n2].insert(n1);
        } else if (edge == "<->") {
            children[n1].insert(n2);
            parents[n2].insert(n1);

            children[n2].insert(n1);
            parents[n1].insert(n2);
        }

    }

    // Get blankets of every node
    for (std::string target : nodes) {
        std::unordered_set<std::string> rule1(parents[target]);

        std::unordered_set<std::string> rule2(children[target]);

        std::unordered_set<std::string> rule3;
        for (std::string child : children[target]) {
            rule3.insert(parents[child].begin(), parents[child].end());
        }

        std::unordered_set<std::string> rule4(undirectedNeighbors[target]);
        for (std::string Y : undirectedNeighbors[target]) {
            rule4.insert(parents[Y].begin(), parents[Y].end());
        }

        // std::unordered_set<std::string> rule6;
        // for (std::string Y : children[target]) {
        //     for (std::string Z : undirectedNeighbors[Y]) {
        //         rule6.insert(Z);
        //         rule6.insert(parents[Z].begin(), parents[Z].end()); // W
        //     }
        // }

        blankets[target].insert(rule1.begin(), rule1.end());
        blankets[target].insert(rule2.begin(), rule2.end());
        blankets[target].insert(rule3.begin(), rule3.end());
        blankets[target].insert(rule4.begin(), rule4.end());
        // blankets[target].insert(rule6.begin(), rule6.end());
        blankets[target].erase(target);
    }

    // Convert to list
    Rcpp::List result = Rcpp::List::create();

    for (std::string n : nodes) {
        result[n] = std::vector<std::string>(blankets[n].begin(), blankets[n].end());
    }

    return result;
}

// NOTE - This function does not need to be exported (because it is called as a helper automatically every
// time a graph is returned) but it could be later by adding "//[[Rcpp::export]]" to the bottom
//' Caclulate the Markov Blanket for every node in the graph.
//' This is done by default whenever a graph is returned from an algorithm,
//' but it can also be used for graphs loaded from files or adj. mats.
//'
//' @param list The graph object
//' @export
//' @examples
//' mat <- matrix(sample(c(0,1), 16, replace=TRUE), nrow=4)
//' nodes <- c("X1", "X2", "X3", "X4")
//' g <- rCausalMGM::adjMat2Graph(mat, nodes, directed=TRUE)
//' g[["markov.blankets"]] <- rCausalMGM::calculateMarkovBlankets(g)
Rcpp::List calculateMarkovBlankets(const Rcpp::List& graph) {
    if (!EdgeListGraph::validateGraphList(graph)) {
        throw std::invalid_argument("ERROR: list is not in the form of a graph");
    }

    if (Rcpp::as<std::string>(graph["type"]) == "markov equivalence class")
      return markovBlanketMEC(graph);
    if (Rcpp::as<std::string>(graph["type"]) == "partial ancestral graph")
      return markovBlanketPAG(graph);
    else
      return markovBlanketUndirected(graph);
}

bool EdgeListGraph::isParentOf(Variable* node1, Variable* node2) {
    for (Edge edge : getEdges(node1)) {
        Variable* sub = Edge::traverseDirected(node1, edge);

        if (sub == node2) {
            return true;
        }
    }

    return false;
}

/**
 * Nodes adjacent to the given node with the given proximal endpoint.
 */
std::vector<Variable*> EdgeListGraph::getNodesInTo(Variable* node, Endpoint endpoint) {
    std::vector<Variable*> nodes;
    std::vector<Edge> edges = getEdges(node);

    for (Edge edge1 : edges) {
        if (edge1.getProximalEndpoint(node) == endpoint) {
            nodes.push_back(edge1.getDistalNode(node));
        }
    }

    return nodes;
}

/**
 * Nodes adjacent to the given node with the given distal endpoint.
 */
std::vector<Variable*> EdgeListGraph::getNodesOutTo(Variable* node, Endpoint endpoint) {
    std::vector<Variable*> nodes;
    std::vector<Edge> edges = getEdges(node);

    for (Edge edge1 : edges) {
        if (edge1.getDistalEndpoint(node) == endpoint) {
            nodes.push_back(edge1.getDistalNode(node));
        }
    }

    return nodes;
}

/**
 * Determines whether one node is an ancestor of another.
 */
bool EdgeListGraph::isAncestorOf(Variable* node1, Variable* node2) {
    std::vector<Variable*> tempList;
    tempList.push_back(node2);
    std::unordered_set<Variable*> ancestors = getAncestors(tempList);
    return (std::find(ancestors.begin(), ancestors.end(), node1) != ancestors.end());
}

std::unordered_set<Variable*> EdgeListGraph::getAncestors(std::vector<Variable*>& nodes) {
    std::unordered_set<Variable*> ancestors;

    for (Variable* node1 : nodes) {
        collectAncestorsVisit(node1, ancestors);
    }
    return ancestors;
}

void EdgeListGraph::collectAncestorsVisit(Variable* node, std::unordered_set<Variable*> &ancestors) {
    if (std::find(ancestors.begin(), ancestors.end(), node) != ancestors.end()) return;

    ancestors.insert(node);
    std::vector<Variable*> parents = getParents(node);

    if (!parents.empty()) {
        for (Variable* parent : parents) {
            collectAncestorsVisit(parent, ancestors);
        }
    }
}

bool EdgeListGraph::isDefCollider(Variable* node1, Variable* node2, Variable* node3) {
    if (!(isAdjacentTo(node1,node2) && isAdjacentTo(node2,node3))) {
        return false;
    }
    Edge edge1 = getEdge(node1, node2);
    Edge edge2 = getEdge(node2, node3);

    return edge1.getProximalEndpoint(node2) == ENDPOINT_ARROW && edge2.getProximalEndpoint(node2) == ENDPOINT_ARROW;

}

void EdgeListGraph::reorientAllWith(Endpoint endpoint) {
  for (Edge edge : getEdges()) {
      Variable* a = edge.getNode1();
      Variable* b = edge.getNode2();

      setEndpoint(a, b, endpoint);
      setEndpoint(b, a, endpoint);
  }
}

Rcpp::List EdgeListGraph::toList() {
    std::vector<std::string> nodeNames;
    for (Variable* node : nodes) {
        nodeNames.push_back(node->getName());
    }

    std::vector<std::string> edgeStrings;
    std::vector<Edge> edges = getEdgeList();
    Edge::sortEdges(edges); //TODO - commented out for testing
    for (Edge edge: edges) {
        edgeStrings.push_back(edge.toString());
    }

    std::vector<std::string> ambiguousTriplesStrings;
    for (Triple t : ambiguousTriples) {
        ambiguousTriplesStrings.push_back(t.toString());
    }
    std::sort(ambiguousTriplesStrings.begin(),
	      ambiguousTriplesStrings.end());

    Rcpp::List result = Rcpp::List::create(
        Rcpp::_["nodes"] = nodeNames,
        Rcpp::_["edges"] = edgeStrings,
        Rcpp::_["ambiguous_triples"] = ambiguousTriplesStrings,
        Rcpp::_["algorithm"] = algorithm,
        Rcpp::_["type"] = graph_type,
        Rcpp::_["markov.blankets"] = R_NilValue,
        Rcpp::_["stabilities"] = R_NilValue
    );

    result.attr("class") = "graph";
    result["markov.blankets"] = calculateMarkovBlankets(result);

    return result;

}

std::ostream& operator<<(std::ostream& os, EdgeListGraph& graph) {

    os << "Graph Nodes:\n";
    std::vector<Variable*> nodes = graph.getNodes();
    int size = nodes.size();
    int count = 0;
    for (Variable* node : nodes) {
        count++;
        os << node->getName();
        if (count < size) {
            os << ",";
        }
    }
    os << "\n\n";

    os << "Graph Edges:\n";
    std::vector<Edge> edges = graph.getEdgeList();
    Edge::sortEdges(edges);
    count = 1;

    for (Edge edge : edges) {
        os << count << ". " << edge << "\n";
        count++;
    }

    os << "\n\n";

    if (graph.ambiguousTriples.size() > 0) {
        os << "Ambiguous triples (i.e. list of triples for which there is ambiguous data about whether they are colliders or not):\n";
        for (Triple t : graph.ambiguousTriples) {
            os << t << "\n";
        }
    }


    return os;
}

// Helper
void streamGraph(const Rcpp::List& list, std::ostream& os) {
    os << "Graph Nodes:\n";

    // Nodes
    std::vector<std::string> nodeNames = list["nodes"];
    int count = 0;
    for (std::string nodeName : nodeNames) {
        count++;
        os << nodeName;
        if (count < nodeNames.size()) {
            os << ",";
        }
    }

    os << "\n\n";

    // Edges
    os << "Graph Edges:\n";
    std::vector<std::string> edgeStrings = list["edges"];
    count = 1;
    for (std::string edgeString : edgeStrings) {
        os << count << ". " << edgeString << "\n";
        count++;
    }

    os << "\n";

    // Algorithm
    std::vector<std::string> algorithms = list["algorithm"];
    os << "Algorithm: " << algorithms[0] << "\n";

    // Type
    std::vector<std::string> t = list["type"];
    os << "Type: " << t[0];

    os << "\n\n";

    //Triples
    std::vector<std::string> tripleStrings = list["ambiguous_triples"];
    if (tripleStrings.size() > 0) {
        os << "Ambiguous triples (i.e. list of triples for which there is ambiguous data about whether they are colliders or not):\n";

        for (std::string tripleString : tripleStrings) {
            os << tripleString << "\n";
        }
    }
}

//' Save a graph to a file
//'
//' @param list The graph object
//' @param filename The graph file
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::mgm(data.n100.p25)
//' rCausalMGM::saveGraph(g, "graphs/mgm_graph.txt")
// [[Rcpp::export]]
void saveGraph(const Rcpp::List& list, const std::string& filename) {
    if (!EdgeListGraph::validateGraphList(list)) {
        throw std::invalid_argument("ERROR: list is not in the form of a graph");
    }

    std::ofstream outfile;
    outfile.open(filename, std::ios::out);

    streamGraph(list, outfile);

    outfile.close();
}

//TODO - Include example graphs and example use of loadGraph()
//' Load a graph from a file
//'
//' @param filename The graph file
//' @return The graph as a List object, which can be passed into search functions
//' @export
// [[Rcpp::export]]
Rcpp::List loadGraph(const std::string& filename) {
    // Get lines from file
    std::vector<std::string> lines;
    try {
        std::ifstream f(filename);

        if(!f) {
            Rcpp::Rcout << "ERROR: Cannot open " << filename << std::endl;
            throw std::invalid_argument("Error opening file: " + filename);
        }
        std::string line;

        while (std::getline(f,line)) {
            if (line.size() > 1 && line.at(line.size()-1) == '\r')
                line = line.substr(0, line.size()-1);
            lines.push_back(line);
        }
    }
    catch(const std::exception& ex) {
        Rcpp::Rcout << "Exception: '" << ex.what() << "'!" << std::endl;
        throw std::invalid_argument("Error reading file: " + filename);
    }

    std::vector<std::string>   nodeNames = GraphUtils::splitString(lines[1], ";");
    if (nodeNames.size() <= 1) nodeNames = GraphUtils::splitString(lines[1], ","); // Check for comma delimiter
    if (nodeNames.at(nodeNames.size()-1) == "") nodeNames.pop_back();

    std::vector<std::string> edgeStrings;
    std::vector<std::string> ambiguousTriplesStrings;
    std::string algorithm = "";
    std::string graph_type = "";

    // Start the line after you read 'Graph Edges:'
    // auto edgeStart = std::find(lines.begin(), lines.end(), "Graph Edges:");
    // if (edgeStart == lines.end())
    //     throw std::invalid_argument("Unable to find 'Graph Edges:' line in " + filename);
    // int edgeStartIndex = edgeStart - lines.begin();
    // Rcpp::Rcout << "edgeStartIndex = " << edgeStartIndex << std::endl;

    // If there are edges
    if (lines.size() > 4) {

        int i = 4; // TODO: Hardcoded for now

        // return true if the string is only whitespace
        auto isWhiteSpace = [](const std::string& _s) {
            std::string s = _s;
            s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
            return s == "";
        };

        for (; i < lines.size(); i++) {
            std::string edgeString = lines[i];

            if (edgeString.find("Algorithm: ") != std::string::npos) {
                algorithm = edgeString.substr(11);
                continue;
            }

            if (edgeString.find("Type: ") != std::string::npos) {
                graph_type = edgeString.substr(6);
                continue;
            }

            if (edgeString.find("Ambiguous triples") != std::string::npos) {
                i++;
                goto TRIPLES;
            }

            if (isWhiteSpace(edgeString)) continue; // Skip empty lines

            if (edgeString.find(". ") == std::string::npos)
                throw std::invalid_argument("Error reading graph " + filename + ": edge is not formatted correctly: " + edgeString);

            edgeString = GraphUtils::splitString(edgeString, ". ")[1];

            // Returns true if the edge string is valid
            auto validateEdgeString = [](const std::string& edgeString) {
                std::vector<std::string> elements = GraphUtils::splitString(edgeString, " ");
                if (elements.size() != 3) return false;
                if (elements[1].at(0) != '<' && elements[1].at(0) != '-' && elements[1].at(0) != 'o') return false;
                if (elements[1].at(1) != '-')                                                         return false;
                if (elements[1].at(2) != '>' && elements[1].at(2) != '-' && elements[1].at(2) != 'o') return false;
                return true;
            };

            if (!validateEdgeString(edgeString))
                throw std::invalid_argument("Error reading graph " + filename + ": edge is not formatted correctly: " + edgeString);

            edgeStrings.push_back(edgeString);
        }

        TRIPLES:
        for (; i < lines.size(); i++) {
            std::string tripleString = lines[i];

            if (isWhiteSpace(tripleString)) continue; // Skip empty lines

            // tripleString = tripleString.substr(0, tripleString.size()-1); // Strip off '\r'

            auto validateTripleString = [](const std::string& tripleString) {
                if (tripleString.size() < 3) return false;
                std::string s = tripleString.substr(1, tripleString.size()-2);
                if (GraphUtils::splitString(s, ",").size() != 3) return false;
                return true;
            };

            if (!validateTripleString(tripleString))
                throw std::invalid_argument("Error reading graph " + filename + ": triple is not formatted correctly: " + tripleString);

            ambiguousTriplesStrings.push_back(tripleString);
        }
    }

    END:
    Rcpp::List result = Rcpp::List::create(
        Rcpp::_["nodes"] = nodeNames,
        Rcpp::_["edges"] = edgeStrings,
        Rcpp::_["ambiguous_triples"] = ambiguousTriplesStrings,
        Rcpp::_["algorithm"] = algorithm,
        Rcpp::_["type"] = graph_type,
        Rcpp::_["markov.blankets"] = R_NilValue,
        Rcpp::_["stabilities"] = R_NilValue
    );

    result.attr("class") = "graph";
    result["markov.blankets"] = calculateMarkovBlankets(result);

    return result;
}

//' Convert an adjacency matrix into a graph
//'
//' @param adj The adjacency matrix, NxN
//' @param nodes The names of the nodes, length N
//' @param directed TRUE if the graph should be directed, default FALSE
//' @return The graph representation of the adjacency list
//' @export
//' @examples
//' mat <- matrix(sample(c(0,1), 16, replace=TRUE), nrow=4)
//' nodes <- c("X1", "X2", "X3", "X4")
//' g <- rCausalMGM::adjMat2Graph(mat, nodes, directed=TRUE)
// [[Rcpp::export]]
Rcpp::List adjMat2Graph(arma::mat adj,
		  Rcpp::StringVector nodes,
		  Rcpp::LogicalVector directed = Rcpp::LogicalVector::create(0) // FALSE
    ) {
    std::vector<std::string> nodeNames(nodes.begin(), nodes.end());

    if (adj.n_rows != adj.n_cols || adj.n_rows <= 0) {
	throw std::invalid_argument("Input adjacency matrix is invalid");
    }

    if (adj.n_rows != nodeNames.size()) {
	throw std::invalid_argument("Input node names do not match the number of variables in the adjacency matrix");
    }

    std::vector<Variable*> _nodes;
    for (int i = 0; i < nodeNames.size(); i++) {
	_nodes.push_back(new ContinuousVariable(nodeNames[i]));
    }

    EdgeListGraph g(_nodes);

    // std::vector<std::string> edgeStrings;
    // std::vector<std::string> ambiguousTriplesStrings;

    bool dir = Rcpp::is_true(Rcpp::all(directed));
    if (dir) {
	for (arma::uword i = 0; i < adj.n_rows; i++) {
	    for (arma::uword j = 0; j < adj.n_rows; j++) {
		if (i==j) continue;
		if (adj(i,j) != 0) {
		    g.addDirectedEdge(_nodes[i], _nodes[j]);
		    // edgeStrings.push_back(nodeNames[i] + " --> " + nodeNames[j]);
		}
	    }
	}
    } else {
	for (arma::uword i = 0; i < adj.n_rows; i++) {
	    for (arma::uword j = i+1; j < adj.n_rows; j++) {
		if (adj(i,j) != 0) {
		    g.addUndirectedEdge(_nodes[i], _nodes[j]);
		    // edgeStrings.push_back(nodeNames[i] + " --- " + nodeNames[j]);
		}
	    }
	}
    }

    Rcpp::List result = g.toList();

    // Delete variables
    for (Variable * v : g.getNodes()) {
        delete v;
    }

    return result;
}

//' Display a graph object as text
//'
//' @param graph The graph object
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::mgm(data.n100.p25)
//' rCausalMGM::printGraph(g)
// [[Rcpp::export]]
void printGraph(const Rcpp::List& graph) {
    streamGraph(graph, Rcpp::Rcout);
}
