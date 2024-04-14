#include "EdgeListGraph.hpp"
#include "MeekRules.hpp"
#include "FciOrient.hpp"
#include "SepsetProducer.hpp"

// Node EdgeListGraph::nullNode = Node();

// Used by constructors
void EdgeListGraph::initNamesHash() {
    for (const Node& node: nodes) {
        namesHash[node.getName()] = node;
    }
}

// /**
//  * Constructs a new (empty) EdgeListGraph.
//  */
// EdgeListGraph::EdgeListGraph() {
//     initNamesHash();
// } 

/**
 * Constructs a new graph, with no edges, using the the given variable
 * names.
 */
EdgeListGraph::EdgeListGraph(const std::vector<Node>& nodes) {
    for (const Node& variable : nodes) {
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
// EdgeListGraph::EdgeListGraph(const EdgeListGraph& graph) {
//     transferNodesAndEdges(graph);
//     ambiguousTriples = graph.ambiguousTriples;
//     underLineTriples = graph.underLineTriples;
//     dottedUnderLineTriples = graph.dottedUnderLineTriples;

//     // for (Edge edge : graph.edgesSet) {
//     //     if (graph.highlightedEdges.count(edge)) {
//     //         setHighlighted(edge, true);
//     //     }
//     // }

//     namesHash = graph.namesHash;
//     hyperparamHash = graph.hyperparamHash;
//     algorithm = graph.algorithm;
//     graph_type = graph.graph_type;
// }

// EdgeListGraph& EdgeListGraph::operator=(const EdgeListGraph& graph) {
//     transferNodesAndEdges(graph);
//     ambiguousTriples = graph.ambiguousTriples;
//     underLineTriples = graph.underLineTriples;
//     dottedUnderLineTriples = graph.dottedUnderLineTriples;

//     // for (Edge edge : graph.edgesSet) {
//     //     if (graph.highlightedEdges.count(edge)) {
//     //         setHighlighted(edge, true);
//     //     }
//     // }

//     namesHash = graph.namesHash;
//     hyperparamHash = graph.hyperparamHash;
//     algorithm = graph.algorithm;
//     graph_type = graph.graph_type;
//     return *this;
// }


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
// EdgeListGraph::EdgeListGraph(EdgeListGraph&& graph) {
//     transferNodesAndEdges(graph);
//     ambiguousTriples = graph.ambiguousTriples;
//     underLineTriples = graph.underLineTriples;
//     dottedUnderLineTriples = graph.dottedUnderLineTriples;

//     // for (Edge edge : graph.edgesSet) {
//     //     if (graph.highlightedEdges.count(edge)) {
//     //         setHighlighted(edge, true);
//     //     }
//     // }
    
//     namesHash = graph.namesHash;
//     hyperparamHash = graph.hyperparamHash;
//     algorithm = graph.algorithm;
//     graph_type = graph.graph_type;
// }

// EdgeListGraph& EdgeListGraph::operator=(EdgeListGraph&& graph) {
//     transferNodesAndEdges(graph);
//     ambiguousTriples = graph.ambiguousTriples;
//     underLineTriples = graph.underLineTriples;
//     dottedUnderLineTriples = graph.dottedUnderLineTriples;

//     // for (Edge edge : graph.edgesSet) {
//     //     if (graph.highlightedEdges.count(edge)) {
//     //         setHighlighted(edge, true);
//     //     }
//     // }
    
//     namesHash = graph.namesHash;
//     hyperparamHash = graph.hyperparamHash;
//     algorithm = graph.algorithm;
//     graph_type = graph.graph_type;
//     return *this;
// }

EdgeListGraph::EdgeListGraph(const Rcpp::List& graph, DataSet& ds)  {
    Rcpp::List list(graph);
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

    arma::vec lambda;
    arma::vec alpha;

    if (!Rf_isNull(list["lambda"])) {
	lambda = Rcpp::as<arma::vec>(list["lambda"]);
    }

    if (!Rf_isNull(list["alpha"])) {
	alpha = Rcpp::as<arma::vec>(list["alpha"]);
    }

    hyperparamHash["lambda"] = lambda;
    hyperparamHash["alpha"] = alpha;
    // hyperparamHash["penalty"] = Rcpp::as<Rcpp::Nullable<Rcpp::NumericVector>>(list["penalty"]);
    
}


EdgeListGraph::EdgeListGraph(const Rcpp::List& graph)  {
    Rcpp::List list(graph);
    if (!validateGraphList(list)) {
        throw std::invalid_argument("ERROR: list is not in the form of a graph");
    }

    // Nodes
    std::vector<std::string> nodeNames = list["nodes"];
    for (std::string nodeName : nodeNames) {
        try {
	  if(!addNode(Node(new ContinuousVariable(nodeName))))
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

    arma::vec lambda;
    arma::vec alpha;

    if (!Rf_isNull(list["lambda"])) {
	lambda = Rcpp::as<arma::vec>(list["lambda"]);
    }

    if (!Rf_isNull(list["alpha"])) {
	alpha = Rcpp::as<arma::vec>(list["alpha"]);
    }

    hyperparamHash["lambda"] = lambda;
    hyperparamHash["alpha"] = alpha;
    
    // hyperparamHash["penalty"] = Rcpp::as<Rcpp::Nullable<Rcpp::NumericVector>>(list["penalty"]);
    
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
    for (const Node& node : graph.nodes) {
        if (!addNode(node))
            throw std::invalid_argument("Problem copying graph nodes");
    }

    for (const Edge& edge : graph.edgesSet) {
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
bool EdgeListGraph::addUndirectedEdge(const Node& node1, const Node& node2) {
    Edge newEdge = Edge::undirectedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds a directed edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addDirectedEdge(const Node& node1, const Node& node2) {
    Edge newEdge = Edge::directedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds a bidirected edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addBidirectedEdge(const Node& node1, const Node& node2) {
    Edge newEdge = Edge::bidirectedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds a partially oriented edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addPartiallyOrientedEdge(const Node& node1, const Node& node2) {
    Edge newEdge = Edge::partiallyOrientedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds a nondirected edge to the graph from node A to node B.
 *
 * @param node1 the "from" node.
 * @param node2 the "to" node.
 */
bool EdgeListGraph::addNondirectedEdge(const Node& node1, const Node& node2) {
    Edge newEdge = Edge::nondirectedEdge(node1, node2);
    return addEdge(newEdge);
}

/**
 * Adds an edge to the graph.
 *
 * @param edge the edge to be added
 * @return true if the edge was added, false if not.
 */
bool EdgeListGraph::addEdge(Edge edge) {

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

    Node node1 = getNode(edgeComponents[0]);
    Node node2 = getNode(edgeComponents[2]);

    if (node1.isNull())
        throw std::invalid_argument("Edge node not found in graph: " + edgeComponents[0]);

    if (node2.isNull())
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
bool EdgeListGraph::addNode(Node node) {
    if (node.isNull())
        throw std::invalid_argument("Can't add NULL node to graph");
  
    // If nodes contains node
    if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) return true;
    
    if (!(getNode(node.getName()).isNull())) {
        if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) {
            namesHash[node.getName()] = node;
        }
    }

    if (edgeLists.count(node)) return false;

    edgeLists[node] = {};
    nodes.push_back(node);
    namesHash[node.getName()] = node;

    return true;
}

/**
 * Adds a node to the graph. Precondition: The proposed name of the node
 * cannot already be used by any other node in the same graph.
 *
 * @param node the node to be added.
 * @return true if the the node was added, false if not.
 */
// bool EdgeListGraph::addNode(Node&& node) {
//     // If nodes contains node
//     if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) return true;

//     if (node.isNull())
//         throw std::invalid_argument("Can't add NULL node to graph");

//     if (!(getNode(node.getName()).isNull())) {
//         if (std::find(nodes.begin(), nodes.end(), node) != nodes.end()) {
//             namesHash[node.getName()] = node;
//         }
//     }

//     if (edgeLists.count(node)) return false;

//     edgeLists[node] = {};
//     nodes.push_back(node);
//     namesHash[node.getName()] = node;

//     return true;
// }

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
bool EdgeListGraph::removeEdges(const Node& node1, const Node& node2) {
    return removeEdges(getEdges(node1, node2));
}

/**
 * Removes all edges connecting node A to node B.
 *
 * @param node1 the first node.,
 * @param node2 the second node.
 * @return true if edges were removed between A and B, false if not.
 */
bool EdgeListGraph::removeEdges(Node&& node1, Node&& node2) {
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
bool EdgeListGraph::removeEdge(const Node& node1, const Node& node2) {
    std::vector<Edge> edges = getEdges(node1, node2);

    if (edges.size() > 1) {
        throw std::invalid_argument("There is more than one edge between " + node1.getName() + " and " + node2.getName());
    }

    return removeEdges(edges);
}

/**
 * Removes the edge connecting the two given nodes.
 */
bool EdgeListGraph::removeEdge(Node&& node1, Node&& node2) {
    std::vector<Edge> edges = getEdges(node1, node2);

    if (edges.size() > 1) {
        throw std::invalid_argument("There is more than one edge between " + node1.getName() + " and " + node2.getName());
    }

    return removeEdges(edges);
}

/**
 * Determines whether some edge or other exists between two nodes.
 */
bool EdgeListGraph::isAdjacentTo(const Node& node1, const Node& node2) {
    if (node1.isNull() || node2.isNull() ||
        edgeLists.find(node1) == edgeLists.end() ||
        edgeLists.find(node2) == edgeLists.end()) {
            return false;
    }

    for (Edge edge : edgeLists[node1]) {
        if (Edge::traverse(node1, edge) == node2) return true;
    }

    return false;
}

bool EdgeListGraph::isDirectedFromTo(const Node& node1, const Node& node2) {
    std::vector<Edge> edges = getEdges(node1, node2);
    if (edges.size() != 1) return false;
    Edge edge = edges[0];
    return edge.pointsTowards(node2);
}

bool EdgeListGraph::isUndirectedFromTo(const Node& node1, const Node& node2) {

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
std::vector<Node> EdgeListGraph::getAdjacentNodes(const Node& node) const {
    std::vector<Edge> edges = edgeLists.at(node);
    std::set<Node> adj;

    for (Edge edge : edges) {
        Node z = edge.getDistalNode(node);
        adj.insert(z);
    }

    return std::vector<Node>(adj.begin(), adj.end());
}

// std::vector<Node> EdgeListGraph::getMarkovBlanket(const Node& node) const {
//     std::vector<Node> mb = getAdjacentNodes(node);
//     std::set<Node> mbSet(mb.begin(), mb.end());
    
//     if (graph_type=="completed partially directed acyclic graph" ||
// 	graph_type=="directed acyclic graph") {

// 	std::vector<Node> colliderSet = getConnectedColliders(node);

// 	mbSet.insert(colliderSet.begin(), colliderSet.end());
	
// 	// for (const Node& n1 : getChildren(node)) {
// 	//     for (const Node& n2 : getParents(n1)) {
// 	// 	if (n2 != node) {
// 	// 	    mbSet.insert(n2);
// 	// 	}
// 	//     }
// 	// }
//     } if (graph_type == "partial ancestral graph" ||
// 	  graph_type == "maximal ancestral graph") {

// 	std::vector<Node> colliderSet = getConnectedColliders(node);

// 	mbSet.insert(colliderSet.begin(), colliderSet.end());
	
	
// 	// for (const Node& n1 : getPossibleChildren(node)) {
// 	//     for (const Node& n2 : getPossibleParents(n1)) {
// 	// 	if (n2 != node) {
// 	// 	    mbSet.insert(n2);
// 	// 	}
// 	//     }
// 	// }

// 	// for (const Node& n1 : getBidirectedConnectedComponent(node)) {
// 	//     if (n1 != node) {
// 	// 	mbSet.insert(n1);
// 	//     }
// 	//     for (const Edge& e : edgeLists.at(n1)) {
// 	// 	if (e.getProximalEndpoint(n1) == ENDPOINT_ARROW) {
// 	// 	    Node n2 = e.getDistalNode(n1)
// 	// 	    if (n2 != node) {
// 	// 	    	mbSet.insert(n2);
// 	// 	    }
// 	// 	}
// 	//     }
// 	// }
//     }

//     mb = std::vector<Node>(mbSet.begin(), mbSet.end());

//     return mb;
// }

std::vector<Node> EdgeListGraph::getBidirectedConnectedComponent(const Node& node) const {

    std::queue<Node> Q;

    std::set<Node> visited, C;

    Node n1, n2;

    int depth = 0, timeToInc = 1;

    // Rcpp::Rcout << "bidirCC for " << node << ":\n";
    
    Q.push(node);
    visited.insert(node);

    // Rcpp::Rcout << "  depth " << depth << "  :  ";

    while (!Q.empty()) {
	n1 = Q.front();
	Q.pop();

	// Rcpp::Rcout << n1 << " ";

	if (depth != 0) C.insert(n1);
	
	for (Edge e : edgeLists.at(n1)) {
	    if (e.isBidirected()) {
		n2 = e.getDistalNode(n1);
		if (visited.count(n2)==0) {
		    
		    visited.insert(n2);
		    Q.push(n2);

		}
	    }
	}

	timeToInc--;

	if (timeToInc==0) {
	    depth++;

	    if (depth > 2) break;

	    // Rcpp::Rcout << "\n  depth " << depth << "  :  ";

	    timeToInc = Q.size();
	}
    }

    // Rcpp::Rcout << "\n";

    std::vector<Node> bidirC(C.begin(), C.end());

    return bidirC;
}

std::vector<Node> EdgeListGraph::getMarkovBlanket(const Node& node) const {

    std::queue<Node> Q;

    std::set<Node> visited, mb;

    Node n1, n2;

    int depth = 0, timeToInc = 1;

    // Rcpp::Rcout << "MB for " << node << ":\n";
    
    Q.push(node);
    visited.insert(node);

    // Rcpp::Rcout << "  depth " << depth << "  :  ";

    while (!Q.empty()) {
	n1 = Q.front();
	Q.pop();

	// Rcpp::Rcout << "\n    " << n1 << " ";
	
	for (Edge e : edgeLists.at(n1)) {
	    if (n1 == node) {
		n2 = e.getDistalNode(n1);
		// Rcpp::Rcout << "\n      " << e << " ";
		mb.insert(n2);
		if (e.getDistalEndpoint(n1) == ENDPOINT_ARROW) {
		    visited.insert(n2);
		    Q.push(n2);
		}
	    } else if (e.getProximalEndpoint(n1) == ENDPOINT_ARROW) {
		n2 = e.getDistalNode(n1);
		if (n2 != node) mb.insert(n2);
		// Rcpp::Rcout << "\n      " << e << " ";
		if (e.isBidirected() && visited.count(n2)==0) {
		    visited.insert(n2);
		    Q.push(n2);
		}
	    }
	}

	timeToInc--;

	if (timeToInc==0) {
	    depth++;

	    if (depth > 2) break;
	    
	    // Rcpp::Rcout << "\n  depth " << depth << "  :  ";

	    timeToInc = Q.size();
	}
    }

    // Rcpp::Rcout << "\n";

    std::vector<Node> mbVec(mb.begin(), mb.end());

    return mbVec;
}


/**
 * @return the edge connecting node1 and node2, provided a unique such edge
 * exists.
 *
 * Throws std::invalid_argument if not
 */
Edge EdgeListGraph::getEdge(const Node& node1, const Node& node2) {

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

    throw std::invalid_argument("node1 and node2 not connected by edge. node1: " + node1.getName() + " node2: " + node2.getName());

}

/**
 * @return the edges connecting node1 and node2.
 */
std::vector<Edge> EdgeListGraph::getEdges(const Node& node1, const Node& node2) {
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
Endpoint EdgeListGraph::getEndpoint(const Node& node1, const Node& node2) {
    std::vector<Edge> edges = getEdges(node2);

    for (Edge edge : edges) {
        if (edge.getDistalNode(node2) == node1) return edge.getProximalEndpoint(node2);
    }

    return ENDPOINT_NULL;
}

/**
 * @return the list of parents for a node.
 */
std::vector<Node> EdgeListGraph::getParents(const Node& node) const {
    std::vector<Node> parents;
    std::vector<Edge> edges = edgeLists.at(node);

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
 * @return the list of possible parents for a node.
 */
std::vector<Node> EdgeListGraph::getPossibleParents(const Node& node) const {
    std::vector<Node> possParents;
    std::vector<Edge> edges = edgeLists.at(node);

    for (Edge edge : edges) {
        Endpoint endpoint1 = edge.getDistalEndpoint(node);
        Endpoint endpoint2 = edge.getProximalEndpoint(node);

        if (endpoint1 != ENDPOINT_ARROW && endpoint2 == ENDPOINT_ARROW) {
            possParents.push_back(edge.getDistalNode(node));
        }
    }

    return possParents;
}


/**
 * @return the list of children for a node.
 */
std::vector<Node> EdgeListGraph::getChildren(const Node& node) const {
    std::vector<Node> children;
    std::vector<Edge> edges = edgeLists.at(node);

    for (Edge edge : edges) {
        Endpoint endpoint1 = edge.getDistalEndpoint(node);
        Endpoint endpoint2 = edge.getProximalEndpoint(node);

        if (endpoint1 == ENDPOINT_ARROW && endpoint2 == ENDPOINT_TAIL) {
            children.push_back(edge.getDistalNode(node));
        }
    }

    return children;
}


/**
 * @return the list of possible children for a node.
 */
std::vector<Node> EdgeListGraph::getPossibleChildren(const Node& node) const {
    std::vector<Node> possChildren;
    std::vector<Edge> edges = edgeLists.at(node);

    for (Edge edge : edges) {
	// Endpoint endpoint1 = edge.getDistalEndpoint(node);
        Endpoint endpoint2 = edge.getProximalEndpoint(node);

        if (endpoint2 != ENDPOINT_ARROW) {
            possChildren.push_back(edge.getDistalNode(node));
        }
    }

    return possChildren;
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
bool EdgeListGraph::setEndpoint(const Node& from, const Node& to, Endpoint endPoint) {
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

    Node nodeX = getNode(nodeNames[0]);
    Node nodeY = getNode(nodeNames[1]);
    Node nodeZ = getNode(nodeNames[2]);

    if (nodeX.isNull())
        throw std::invalid_argument("Edge node not found in graph: " + nodeNames[0]);

    if (nodeY.isNull())
        throw std::invalid_argument("Edge node not found in graph: " + nodeNames[1]);

    if (nodeZ.isNull())
        throw std::invalid_argument("Edge node not found in graph: " + nodeNames[2]);
    
    return Triple(nodeX, nodeY, nodeZ);
}

bool EdgeListGraph::validateGraphList(Rcpp::List& l) {
    std::vector<std::string> lclass = l.attr("class");
    
    if (std::find(lclass.begin(), lclass.end(), "graph") == lclass.end()) return false;

    std::vector<std::string> names = l.names();

    if (names.size() < 2)                                                 return false;
    if (std::find(names.begin(), names.end(), "nodes") == names.end())    return false;
    if (std::find(names.begin(), names.end(), "edges") == names.end())    return false;
    if (std::find(names.begin(), names.end(), "type") == names.end())    return false;
    if (std::find(names.begin(), names.end(), "ambiguous_triples") == names.end())
	l["ambiguous_triples"] = std::vector<std::string>();
    if (std::find(names.begin(), names.end(), "algorithm") == names.end())
	l["algorithm"] = "";
    if (std::find(names.begin(), names.end(), "lambda") == names.end())
	l["lambda"] = R_NilValue;
    if (std::find(names.begin(), names.end(), "alpha") == names.end())
	l["alpha"] = R_NilValue;
    if (std::find(names.begin(), names.end(), "markov.blankets") == names.end())
	l["markov.blankets"] = R_NilValue; // calculateMarkovBlankets(l);
    if (std::find(names.begin(), names.end(), "neighbors") == names.end())
	l["neighbors"] = R_NilValue; // calculateMarkovBlankets(l);
    if (std::find(names.begin(), names.end(), "stabilities") == names.end())
	l["stabilities"] = R_NilValue;
    // if (names[2] != "ambiguous_triples") return false;
    // if (names[3] != "algorithm")         return false;
    // if (names[4] != "type")              return false;
    // if (names[5] != "markov.blankets")   return false;
    // if (names[6] != "stabilities")       return false;

    return true;
}

// /**
//  * Calculate markov blankets for undirected graph
//  * For a node x, all neighbors of x are in the Markov Blanket of x
//  */
// Rcpp::List markovBlanketUndirected(const Rcpp::List& graph) {
//     std::vector<std::string> nodes = graph["nodes"];
//     std::vector<std::string> edges = graph["edges"];

//     std::unordered_map<std::string, std::unordered_set<std::string>> blankets;

//     for (std::string n : nodes) {
//         blankets[n] = std::unordered_set<std::string>();
//     }

//     // Get blankets
//     for (std::string edgeString : edges) {
//         std::vector<std::string> e = GraphUtils::splitString(edgeString, " ");

//         std::string n1 = e[0];
//         std::string n2 = e[2];

//         blankets[n1].insert(n2);
//         blankets[n2].insert(n1);
//     }

//     // Convert to list
//     Rcpp::List result = Rcpp::List::create();

//     for (std::string n : nodes) {
//         result[n] = std::vector<std::string>(blankets[n].begin(), blankets[n].end());
//     }

//     return result;
// }

// /**
//  * Calculate markov blankets for Partial Ancestral graphs
//  *
//  * Rules:
//  * 1. Parents, children, and spouses (linked by fully directed edges) are treated the 
//  *    same as in DAGs, and are all included in the Markov blanket.
//  * 2. Unoriented, partially oriented, and bidirected edges all have to be treated like 
//  *    there is latent confounding, because it hasn't been ruled out. Thus, any node 
//  *    connected to the target by unoriented, partially oriented, or 
//  *    bidirected edges are included in the Markov blanket.
//  * 3. Unoriented, partially oriented, or bidirected edges connected to the children are 
//  *    included in the Markov blanket.
//  * 4. For the nodes added to the Markov blanket in rule two, add any parents of those
//  *    nodes to the Markov blanket.
//  * 5. Additionally, add any nodes connected to those added in rule two by unoriented, 
//  *    partially oriented, or bidirected edges. We're going to cut off the Markov blanket 
//  *    here arbitrarily, but in theory rules 3 and 4 should be applied recursively on 
//  *    each new set of nodes connected by unoriented, partially oriented, or bidirected 
//  *    edges. However, this could lead to ridiculously large Markov blankets, and it would 
//  *    be unlikely the additional nodes would actually improve predictive performance.
//  */
// Rcpp::List markovBlanketPAG(const Rcpp::List& graph) {
//     std::vector<std::string> nodes = graph["nodes"];
//     std::vector<std::string> edges = graph["edges"];

//     std::unordered_map<std::string, std::unordered_set<std::string>> blankets;
//     std::unordered_map<std::string, std::unordered_set<std::string>> confoundingNeighbors; // Nodes connected by o->, or <->
//     std::unordered_map<std::string, std::unordered_set<std::string>> parents;
//     std::unordered_map<std::string, std::unordered_set<std::string>> children;
//     std::unordered_map<std::string, std::unordered_set<std::string>> partialParents;
//     std::unordered_map<std::string, std::unordered_set<std::string>> partialChildren;
//     std::unordered_map<std::string, std::unordered_set<std::string>> nondirected;

//     for (std::string n : nodes) {
//         blankets[n] = std::unordered_set<std::string>();
//         partialParents[n] = std::unordered_set<std::string>();
// 	partialChildren[n] = std::unordered_set<std::string>();
// 	confoundingNeighbors[n] = std::unordered_set<std::string>();
//         parents[n] = std::unordered_set<std::string>();
//         children[n] = std::unordered_set<std::string>();
// 	nondirected[n] = std::unordered_set<std::string>();
//     }

//     // Get neighbors of every node
//     for (std::string edgeString : edges) {
//         std::vector<std::string> e = GraphUtils::splitString(edgeString, " ");

//         std::string n1 = e[0];
//         std::string n2 = e[2];
//         std::string edge = e[1];

//         if (edge == "<->") {
//             confoundingNeighbors[n1].insert(n2);
//             confoundingNeighbors[n2].insert(n1);
// 	    // confoundingChildren[n1].insert(n2);
// 	    // confoundingChildren[n2].insert(n1);
// 	} else if (edge == "o->") {
// 	    partialChildren[n1].insert(n2);
//             partialParents[n2].insert(n1);
//         } else if (edge == "-->") {
//             children[n1].insert(n2);
//             parents[n2].insert(n1);
//         } else if (edge == "o-o") {
// 	    nondirected[n1].insert(n2);
// 	    nondirected[n2].insert(n1);
// 	}

//     }

//     // Get blankets of every node
//     for (std::string target : nodes) {

//         std::unordered_set<std::string> spouses;
//         for (std::string child : children[target]) {
//             spouses.insert(parents[child].begin(), parents[child].end());
// 	    spouses.insert(partialParents[child].begin(),
// 			   partialParents[child].end());
//         }

// 	std::unordered_set<std::string> partialSpouses;
//         for (std::string child : partialChildren[target]) {
//             partialSpouses.insert(parents[child].begin(), parents[child].end());
// 	    partialSpouses.insert(partialParents[child].begin(),
// 				  partialParents[child].end());
//         }

//         std::unordered_set<std::string> rule2; // (confoundingNeighbors[target]);
// 	std::unordered_set<std::string> rule3;
//         for (std::string child : children[target]) {
//             rule2.insert(confoundingNeighbors[child].begin(), confoundingNeighbors[child].end());
//         }
	
// 	for (std::string child : partialChildren[target]) {
//             rule3.insert(confoundingNeighbors[child].begin(), confoundingNeighbors[child].end());
//         }

//         std::unordered_set<std::string> rule4;
//         std::unordered_set<std::string> rule5;
// 	std::unordered_set<std::string> rule6;
//         for (std::string Y : confoundingNeighbors[target]) {
//             rule4.insert(parents[Y].begin(), parents[Y].end());
// 	    rule5.insert(partialParents[Y].begin(), partialParents[Y].end());
//             rule6.insert(confoundingNeighbors[Y].begin(), confoundingNeighbors[Y].end());
//         }

// 	std::unordered_set<std::string> rule7;
// 	for (std::string Y : rule2) {
//             rule7.insert(parents[Y].begin(), parents[Y].end());
// 	    rule7.insert(partialParents[Y].begin(), partialParents[Y].end());
//             rule7.insert(confoundingNeighbors[Y].begin(), confoundingNeighbors[Y].end());
//         }

// 	for (std::string Y : rule3) {
//             rule7.insert(parents[Y].begin(), parents[Y].end());
// 	    rule7.insert(partialParents[Y].begin(), partialParents[Y].end());
//             rule7.insert(confoundingNeighbors[Y].begin(), confoundingNeighbors[Y].end());
//         }

// 	for (std::string Y : rule6) {
//             rule7.insert(parents[Y].begin(), parents[Y].end());
// 	    rule7.insert(partialParents[Y].begin(), partialParents[Y].end());
//             rule7.insert(confoundingNeighbors[Y].begin(), confoundingNeighbors[Y].end());
//         }

//         blankets[target].insert(parents[target].begin(), parents[target].end());
//         blankets[target].insert(children[target].begin(), children[target].end());
//         blankets[target].insert(spouses.begin(), spouses.end());
// 	blankets[target].insert(partialParents[target].begin(),
// 				partialParents[target].end());
//         blankets[target].insert(partialChildren[target].begin(),
// 				partialChildren[target].end());
//         blankets[target].insert(partialSpouses.begin(),
// 				partialSpouses.end());
// 	blankets[target].insert(confoundingNeighbors[target].begin(),
// 				confoundingNeighbors[target].end());
// 	blankets[target].insert(nondirected[target].begin(), nondirected[target].end());
	
//         blankets[target].insert(rule2.begin(), rule2.end());
//         blankets[target].insert(rule3.begin(), rule3.end());
//         blankets[target].insert(rule4.begin(), rule4.end());
// 	blankets[target].insert(rule5.begin(), rule5.end());
// 	blankets[target].insert(rule6.begin(), rule6.end());
//         blankets[target].erase(target);
//     }

//     // Convert to list
//     Rcpp::List result = Rcpp::List::create();

//     for (std::string n : nodes) {
//         result[n] = std::vector<std::string>(blankets[n].begin(), blankets[n].end());
//     }

//     return result;
// }


// /**
//  * Calculate markov blankets for Markov Equivalence Class graphs
//  *
//  * Updated rules for how to handle Markov blankets in PDAGs (Partially Directed Acyclic Graphs):
//  * 1. Parents of target node
//  * 2. Children of target node
//  * 3. Spouses of target node
//  * 4. If the target variable X has an undirected edge to Y, then Y and its directed parents
//  *    are in the Markov blanket.
//  * 5. If the target variable X has an undirected edge to Y, we do not include any node Z
//  *    connected to Y by an undirected edge. This can be explained because in the Markov
//  *    equivalence class, two consecutive undirected edges X --- Y --- Z can be either
//  *    X --> Y --> Z, X <-- Y <-- Z, or X <-- Y --> Z. In none of these cases is Z a parent of Y.
//  * 6. In the case where the target variable X is a parent of node Y that contains an
//  *    undirected edge to a node Z, then Y, Z, and any directed parents W of Z are included
//  *    in the Markov blanket.
//  *    #### Rule 6 commented out, TODO: test whether inclusion improves MB inference
//  *
//  */
// Rcpp::List markovBlanketMEC(const Rcpp::List& graph) {
//     std::vector<std::string> nodes = graph["nodes"];
//     std::vector<std::string> edges = graph["edges"];

//     std::unordered_map<std::string, std::unordered_set<std::string>> blankets;
//     std::unordered_map<std::string, std::unordered_set<std::string>> undirectedNeighbors;
//     std::unordered_map<std::string, std::unordered_set<std::string>> parents;
//     std::unordered_map<std::string, std::unordered_set<std::string>> children;

//     for (std::string n : nodes) {
//         blankets[n] = std::unordered_set<std::string>();
//         undirectedNeighbors[n] = std::unordered_set<std::string>();
//         parents[n] = std::unordered_set<std::string>();
//         children[n] = std::unordered_set<std::string>();
//     }

//     // Get neighbors of every node
//     for (std::string edgeString : edges) {
//         std::vector<std::string> e = GraphUtils::splitString(edgeString, " ");

//         std::string n1 = e[0];
//         std::string n2 = e[2];
//         std::string edge = e[1];

//         if (edge == "---") {
//             undirectedNeighbors[n1].insert(n2);
//             undirectedNeighbors[n2].insert(n1);
//         } else if (edge == "-->") {
//             children[n1].insert(n2);
//             parents[n2].insert(n1);
//         } else if (edge == "<->") {
//             children[n1].insert(n2);
//             parents[n2].insert(n1);

//             children[n2].insert(n1);
//             parents[n1].insert(n2);
//         }

//     }

//     // Get blankets of every node
//     for (std::string target : nodes) {
//         std::unordered_set<std::string> rule1(parents[target]);

//         std::unordered_set<std::string> rule2(children[target]);

//         std::unordered_set<std::string> rule3;
//         for (std::string child : children[target]) {
//             rule3.insert(parents[child].begin(), parents[child].end());
//         }

//         std::unordered_set<std::string> rule4(undirectedNeighbors[target]);
//         // for (std::string Y : undirectedNeighbors[target]) {
//         //     rule4.insert(parents[Y].begin(), parents[Y].end());
//         // }

//         // std::unordered_set<std::string> rule6;
//         // for (std::string Y : children[target]) {
//         //     for (std::string Z : undirectedNeighbors[Y]) {
//         //         rule6.insert(Z);
//         //         rule6.insert(parents[Z].begin(), parents[Z].end()); // W
//         //     }
//         // }

//         blankets[target].insert(rule1.begin(), rule1.end());
//         blankets[target].insert(rule2.begin(), rule2.end());
//         blankets[target].insert(rule3.begin(), rule3.end());
//         blankets[target].insert(rule4.begin(), rule4.end());
//         // blankets[target].insert(rule6.begin(), rule6.end());
//         blankets[target].erase(target);
//     }

//     // Convert to list
//     Rcpp::List result = Rcpp::List::create();

//     for (std::string n : nodes) {
//         result[n] = std::vector<std::string>(blankets[n].begin(), blankets[n].end());
//     }

//     return result;
// }

// // NOTE - This function does not need to be exported (because it is called as a helper automatically every
// // time a graph is returned) but it could be later by adding "//[[Rcpp::export]]" to the bottom
// //' Caclulate the Markov Blanket for every node in the graph.
// //' This is done by default whenever a graph is returned from an algorithm,
// //' but it can also be used for graphs loaded from files or adj. mats.
// //'
// //' @param list The graph object
// //' @export
// //' @examples
// //' mat <- matrix(sample(c(0,1), 16, replace=TRUE), nrow=4)
// //' nodes <- c("X1", "X2", "X3", "X4")
// //' g <- adjMat2Graph(mat, nodes, directed=TRUE)
// //' g[["markov.blankets"]] <- calculateMarkovBlankets(g)
// Rcpp::List calculateMarkovBlankets(const Rcpp::List& graph) {
//     EdgeListGraph g(graph);
//     Rcpp::List list(graph);
//     if (!EdgeListGraph::validateGraphList(list)) {
//         throw std::invalid_argument("ERROR: list is not in the form of a graph");
//     }

//     if (Rcpp::as<std::string>(list["type"])=="completed partially directed acyclic graph" ||
// 	Rcpp::as<std::string>(list["type"])=="directed acyclic graph")
// 	return markovBlanketMEC(list);
//     if (Rcpp::as<std::string>(list["type"]) == "partial ancestral graph" ||
// 	Rcpp::as<std::string>(list["type"]) == "maximal ancestral graph")
// 	return markovBlanketPAG(list);
//     else
// 	return markovBlanketUndirected(list);
// }

bool EdgeListGraph::isParentOf(const Node& node1, const Node& node2) {
    for (Edge edge : getEdges(node1)) {
        Node sub = Edge::traverseDirected(node1, edge);

        if (sub == node2) {
            return true;
        }
    }

    return false;
}


// bool EdgeListGraph::existsUndirectedPathVisit(const Node& node1, const Node& node2, std::set<Node>& path);

// bool EdgeListGraph::existsDirectedPathVisit(const Node& node1, const Node& node2,
// 					    std::set<Node>& path){
    
// }


bool EdgeListGraph::existsDirectedPathFromTo(const Node& node1, const Node& node2) {
    std::set<Node> visited;

    // std::list<Edge> path;

    // std::map<Node,Node> parentMap;

    std::stack<Node> S;

    // Rcpp::Rcout << "  Searching for directed path from " << node1 << " to " << node2 << ":\n    ";

    S.push(node1);

    while (!S.empty()) {
	Node v = S.top();
	S.pop();

	// Rcpp::Rcout << v << " ";
	
	visited.insert(v);

	for (Node ch : getChildren(v)) {
	    if (ch==node2) {
		// Rcpp::Rcout << "\n\nDIRECTED PATH FOUND\n\n";

		// parentMap[ch] = v;

		// Node n = ch;
		// do {
		//     Edge e = getEdge(n, parentMap[n]);
		//     n = parentMap[n];
		//     path.push_front(e);
		// } while (n!=node1);
		
		// Rcpp::Rcout << "  Path:\n    ";
		// for (Edge e2 : path) {
		//     Rcpp::Rcout << e2 << " ";
		// }
		// Rcpp::Rcout << "\n\n";
		return true;
	    }
	    if (visited.count(ch)==0) {
		S.push(ch);
		// parentMap[ch] = v;
	    }
	}
    }
    // Rcpp::Rcout << "\n";
    return false;
}

bool EdgeListGraph::existsSemiDirectedPathFromTo(const Node& node1, const Node& node2) {
    std::set<Node> visited;

    std::stack<Node> S;

    // Rcpp::Rcout << "  Searching for directed path from " << node1 << " to " << node2 << ":\n    ";

    S.push(node1);

    while (!S.empty()) {
	Node v = S.top();
	S.pop();

	// Rcpp::Rcout << v << " ";
	
	visited.insert(v);

	for (Node ch : getPossibleChildren(v)) {
	    if (ch==node2) {
		// Rcpp::Rcout << ch << "\n\nNODE FOUND\n\n";
		return true;
	    }
	    if (visited.count(ch)==0) {
		S.push(ch);
	    }
	}
    }
    // Rcpp::Rcout << "\n";
    return false;
}


bool EdgeListGraph::existsAlmostDirectedPathFromTo(const Node& node1, const Node& node2) {
    std::set<Node> visited;

    std::list<Edge> path;

    std::map<Node,Node> parentMap;

    bool almostFlag = false;
    bool growFlag = false;

    std::stack<Node> S;

    // Rcpp::Rcout << "  Searching for almost directed path from " << node1 << " to " << node2 << ":\n    ";

    S.push(node1);

    while (!S.empty()) {
	Node v = S.top();
	S.pop();

	// Rcpp::Rcout << v << " ";

	if (v != node1) {
	    if (!growFlag) {
		path.clear();
		almostFlag = false;

		Node n = v;
		do {
		    Edge e = getEdge(n, parentMap[n]);
		    n = parentMap[n];
		    path.push_front(e);
		    almostFlag = almostFlag || !e.isDirected();
		} while (n!=node1);
	    } else {
		Edge e = getEdge(v, parentMap[v]);
		path.push_back(e);
		almostFlag = almostFlag || !e.isDirected();
	    }
	    // if (e.isDirected()) {
	    // 	path.push_back(e);
	    // } else if (!almostFlag) {
	    // 	path.push_back(e);
	    // 	almostFlag = true;
	    // } else {
	    // 	continue;
	    // }
	}

	growFlag = false;
	
	visited.insert(v);

	for (Edge e : getEdges(v)) {
	    if (e.getDistalEndpoint(v)!=ENDPOINT_ARROW) continue;
	    if (almostFlag && !e.isDirected()) continue;

	    Node n = e.getDistalNode(v);
	    
	    if (n==node2) {		
		// Rcpp::Rcout << "\n\nALMOST DIRECTED PATH FOUND\n\n";
		// Rcpp::Rcout << "  Path:\n    ";
		// for (Edge e2 : path) {
		//     Rcpp::Rcout << e2 << " ";
		// }
		// Rcpp::Rcout << e << "\n\n";
		return true;
	    }
	    
	    if (visited.count(n)==0) {
		S.push(n);
		parentMap[n] = v;
		growFlag = true;
	    }
	}

	// if (!growFlag && !path.empty()) {
	//     if (path.back().isBidirected() || path.back().isPartiallyOriented()) {
	// 	almostFlag = false;
	//     }
	//     path.pop_back();
	// }
    }
    // Rcpp::Rcout << "\n";
    return false;
}


/**
 * Nodes adjacent to the given node with the given proximal endpoint.
 */
std::vector<Node> EdgeListGraph::getNodesInTo(const Node& node, Endpoint endpoint) {
    std::vector<Node> nodes;
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
std::vector<Node> EdgeListGraph::getNodesOutTo(const Node& node, Endpoint endpoint) {
    std::vector<Node> nodes;
    std::vector<Edge> edges = getEdges(node);

    for (Edge edge1 : edges) {
        if (edge1.getDistalEndpoint(node) == endpoint) {
            nodes.push_back(edge1.getDistalNode(node));
        }
    }

    return nodes;
}

/**
 * Determines whether node1 is an ancestor of node2.
 */
bool EdgeListGraph::isAncestorOf(const Node& node1, const Node& node2) {
    return existsDirectedPathFromTo(node1, node2);
    // std::vector<Node> tempList;
    // tempList.push_back(node2);
    // std::unordered_set<Node> ancestors = getAncestors(tempList);
    // return (std::find(ancestors.begin(), ancestors.end(), node1) != ancestors.end());
}

/**
 * Determines whether node1 is an ancestor of node2.
 */
bool EdgeListGraph::isPossibleAncestorOf(const Node& node1, const Node& node2) {
    return existsAlmostDirectedPathFromTo(node1, node2);
    // std::vector<Node> tempList;
    // tempList.push_back(node2);
    // std::unordered_set<Node> ancestors = getAncestors(tempList);
    // return (std::find(ancestors.begin(), ancestors.end(), node1) != ancestors.end());
}


std::unordered_set<Node> EdgeListGraph::getAncestors(std::vector<Node>& nodes) {
    std::unordered_set<Node> ancestors;

    for (const Node& node1 : nodes) {
        collectAncestorsVisit(node1, ancestors);
    }
    return ancestors;
}

void EdgeListGraph::collectAncestorsVisit(const Node& node, std::unordered_set<Node> &ancestors) {
    if (std::find(ancestors.begin(), ancestors.end(), node) != ancestors.end()) return;

    ancestors.insert(node);
    std::vector<Node> parents = getParents(node);

    if (!parents.empty()) {
        for (const Node& parent : parents) {
            collectAncestorsVisit(parent, ancestors);
        }
    }
}

bool EdgeListGraph::isDefCollider(const Node& node1, const Node& node2, const Node& node3) {
    if (!(isAdjacentTo(node1,node2) && isAdjacentTo(node2,node3))) {
        return false;
    }
    Edge edge1 = getEdge(node1, node2);
    Edge edge2 = getEdge(node2, node3);

    return edge1.getProximalEndpoint(node2) == ENDPOINT_ARROW && edge2.getProximalEndpoint(node2) == ENDPOINT_ARROW;

}

void EdgeListGraph::reorientAllWith(Endpoint endpoint) {    
    for (auto&& edge : getEdgeList()) {
	Node a = edge.getNode1();
	Node b = edge.getNode2();

	setEndpoint(a, b, endpoint);
	setEndpoint(b, a, endpoint);
    }
}


bool EdgeListGraph::validSink(const Node& node) {
    std::list<Node> neighbors;
    
    for (Edge edge : edgeLists[node]) {
	if (edge.getDistalEndpoint(node) == ENDPOINT_ARROW)
	    return false;
	else if (edge.getProximalEndpoint(node) == ENDPOINT_TAIL)
	    neighbors.push_back(edge.getDistalNode(node));
    }

    while (!neighbors.empty()) {
	Node y = neighbors.front();
	neighbors.pop_front();
	for (Node z : neighbors)
	    if (!isAdjacentTo(y,z))
		return false;
    }

    return true;
}


std::list<Node> EdgeListGraph::getCausalOrdering() {
    EdgeListGraph tempGraph(*this);
    std::set<Node> nodeSet(nodes.begin(), nodes.end());
    std::list<Node> order;

    while (!nodeSet.empty()) {
	for (Node n : nodeSet) {
	    if (tempGraph.validSink(n)) {
		order.push_front(n);
		tempGraph.removeNode(n);
		nodeSet.erase(n);
		break;
	    }
	}
    }

    if (order.size() != nodes.size()) {
	throw std::runtime_error("The number of nodes in causal ordering does not match number in graph");
    }

    return order;
}


Rcpp::List EdgeListGraph::toList() const {
    std::vector<std::string> nodeNames;
    Rcpp::List mbList = Rcpp::List::create();
    Rcpp::List nbList = Rcpp::List::create();
    
    for (const Node& node : nodes) {
        nodeNames.push_back(node.getName());
	
	std::vector<Node> mb = getMarkovBlanket(node); // already sorted
	std::vector<std::string> mbNames;
	for (const Node& n1 : mb) {
	    mbNames.push_back(n1.getName());
	}
	mbList[node.getName()] = mbNames;
	
	std::vector<Node> nb = getAdjacentNodes(node); // sorted
	// std::sort(nb.begin(), nb.end());
	std::vector<std::string> nbNames;
	for (const Node& n1 : nb) {
	    nbNames.push_back(n1.getName());
	}
	nbList[node.getName()] = nbNames;
    }

    std::vector<std::string> edgeStrings;
    std::vector<Edge> edges = getEdgeList();
    Edge::sortEdges(edges); //TODO - commented out for testing
    for (const Edge& edge: edges) {
        edgeStrings.push_back(edge.toString());
    }

    std::vector<std::string> ambiguousTriplesStrings;
    for (const Triple t : ambiguousTriples) {
        ambiguousTriplesStrings.push_back(t.toString());
    }
    std::sort(ambiguousTriplesStrings.begin(),
	      ambiguousTriplesStrings.end());

    Rcpp::List result = Rcpp::List::create(
        Rcpp::_["nodes"] = nodeNames,
        Rcpp::_["edges"] = edgeStrings,
        Rcpp::_["ambiguous_triples"] = ambiguousTriplesStrings,
        Rcpp::_["algorithm"] = algorithm,
	Rcpp::_["lambda"] = hyperparamHash.at("lambda").is_empty() ? R_NilValue : Rcpp::wrap(hyperparamHash.at("lambda")),
	Rcpp::_["alpha"] = hyperparamHash.at("alpha").is_empty() ? R_NilValue : Rcpp::wrap(hyperparamHash.at("alpha")),
	Rcpp::_["penalty"] = hyperparamHash.at("penalty").is_empty() ? R_NilValue : Rcpp::wrap(hyperparamHash.at("penalty")),
        Rcpp::_["type"] = graph_type,
        Rcpp::_["markov.blankets"] = mbList,
	Rcpp::_["neighbors"] = nbList,
        Rcpp::_["stabilities"] = R_NilValue
    );

    result.attr("class") = "graph";
    // result["markov.blankets"] = calculateMarkovBlankets(result);

    return result;

}

bool operator==(const EdgeListGraph& g1, const EdgeListGraph& g2) {
    if (g1.nodes.size() != g2.nodes.size()) return false;
    std::set<Node> nodeSet1(g1.nodes.begin(), g1.nodes.end());
    std::set<Node> nodeSet2(g2.nodes.begin(), g2.nodes.end());
    return nodeSet1 == nodeSet2 && g1.edgesSet == g2.edgesSet;
}

bool operator!=(const EdgeListGraph& g1, const EdgeListGraph& g2) {
    return !(g1 == g2);
}

bool operator<(const EdgeListGraph& g1, const EdgeListGraph& g2) {
    if (g1.nodes.size() != g2.nodes.size()) return g1.nodes.size() < g2.nodes.size();
    std::set<Node> nodeSet1(g1.nodes.begin(), g1.nodes.end());
    std::set<Node> nodeSet2(g2.nodes.begin(), g2.nodes.end());
    if (nodeSet1 != nodeSet2) return nodeSet1 < nodeSet2;
    return g1.edgesSet < g2.edgesSet;
}

bool operator>=(const EdgeListGraph& g1, const EdgeListGraph& g2) {
    return !(g1 < g2);
}

bool operator<=(const EdgeListGraph& g1, const EdgeListGraph& g2) {
    return (g1 < g2) || (g1 == g2);
}

bool operator>(const EdgeListGraph& g1, const EdgeListGraph& g2) {
    return !(g1 <= g2);
}


std::ostream& operator<<(std::ostream& os, EdgeListGraph graph) {

    os << "Graph Nodes:\n";
    std::vector<Node> nodes = graph.getNodes();
    int size = nodes.size();
    int count = 0;
    for (Node node : nodes) {
        count++;
        os << node.getName();
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
void streamGraph(const Rcpp::List& list, std::ostream& os, std::string ext) {

    if (ext == "sif") {

	// Edges
	std::vector<std::string> edgeStrings = list["edges"];
	for (std::string edgeString : edgeStrings) {
	    std::vector<std::string> edge;
	    size_t pos = edgeString.find(" ");
	    while (pos != std::string::npos) {
		edge.push_back(edgeString.substr(0,pos));
		edgeString.erase(0,pos+1);
		pos = edgeString.find(" ");
	    }
	    edge.push_back(edgeString);
	    
	    // Rcpp::Rcout << "edge.size() = " << edge.size() << std::endl;
	    // for (auto it = edge.begin(); it != edge.end(); it++)
	    // 	Rcpp::Rcout << *it << ", ";
	    // Rcpp::Rcout << "\n";
	    
	    if (edge.size() == 3) {
		os << edge[0] << "\t";

		if (edge[1] == "---") {
		    os << "undir\t";
		} else if (edge[1] == "o-o") {
		    os << "cc\t";
		} else if (edge[1] == "o->") {
		    os << "ca\t";
		} else if (edge[1] == "-->") {
		    os << "dir\t";
		} else if (edge[1] == "<->") {
		    os << "bidir\t";
		} else {
		    throw std::invalid_argument("ERROR: invalid edge type " +
						edge[1] + " in edge list.");
		}
		
		os << edge[2] << "\n";
	    } else {
		throw std::invalid_argument("ERROR: invalid edge not of the form var1 edgetype var2 in edge list.");
	    }
	}

    } else {
    
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

	// Lambda
	Rcpp::Nullable<Rcpp::NumericVector> lambda = list["lambda"];
	if (!lambda.isNull()) {
	    os << "Lambda: " << Rcpp::as<Rcpp::NumericVector>(lambda) << "\n";
	}

	// Alpha
	Rcpp::Nullable<Rcpp::NumericVector> alpha = list["alpha"];
	if (!alpha.isNull()) {
	    os << "Alpha: " << Rcpp::as<Rcpp::NumericVector>(alpha) << "\n";
	}

	// Penalty
	Rcpp::Nullable<Rcpp::NumericVector> penalty = list["penalty"];
	if (!penalty.isNull()) {
	    os << "Penalty: " << Rcpp::as<Rcpp::NumericVector>(penalty) << "\n";
	}
	
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
}

//' Save a graph to a file. Supported file types are ".txt" and ".sif". 
//'
//' @param graph The graph object
//' @param filename The graph filename
//' @export
// [[Rcpp::export]]
void saveGraph(const Rcpp::List& graph, const std::string& filename) {
    Rcpp::List list(graph);
    if (!EdgeListGraph::validateGraphList(list)) {
        throw std::invalid_argument("ERROR: list is not in the form of a graph");
    }

    std::string ext = filename.substr(filename.find_last_of(".") + 1);

    std::string fn(filename);

    if (ext != "txt" && ext != "sif") {
	Rcpp::Rcout << "  Unsupported file type detected. Saving graph as a txt file.\n";
	fn += ".txt";
	ext = "txt";
    }
    
    std::ofstream outfile;
    outfile.open(fn, std::ios::out);

    streamGraph(list, outfile, ext);

    outfile.close();
}

//' Load a graph from a ".txt" file
//'
//' @param filename The graph file
//' @return The graph as a graph object, which can be passed into search functions
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
	    if (line.size() > 0)
		lines.push_back(line);
        }
    }
    catch(const std::exception& ex) {
        Rcpp::Rcout << "Exception: '" << ex.what() << "'!" << std::endl;
        throw std::invalid_argument("Error reading file: " + filename);
    }

    auto nodeStart = std::find(lines.begin(), lines.end(), "Graph Nodes:");
    if (nodeStart == lines.end())
        throw std::invalid_argument("Unable to find 'Graph Nodes:' line in " + filename);
    int nodeStartIndex = nodeStart - lines.begin();
    // Rcpp::Rcout << "nodeStartIndex = " << nodeStartIndex << ": " << lines[nodeStartIndex] << std::endl;

    std::vector<std::string>   nodeNames = GraphUtils::splitString(lines[nodeStartIndex+1],
								   ";");
    if (nodeNames.size() <= 1) nodeNames = GraphUtils::splitString(lines[nodeStartIndex+1],
								   ","); // Check for comma delimiter
    if (nodeNames.at(nodeNames.size()-1) == "") nodeNames.pop_back();

    std::vector<std::string> edgeStrings;
    std::vector<std::string> ambiguousTriplesStrings;
    std::string algorithm = "";
    std::string graph_type = "";
    Rcpp::Nullable<Rcpp::NumericVector> lambda = R_NilValue;
    Rcpp::Nullable<Rcpp::NumericVector> alpha = R_NilValue;
    Rcpp::Nullable<Rcpp::NumericVector> penalty = R_NilValue;
    
    // Start the line after you read 'Graph Edges:'
    auto edgeStart = std::find(lines.begin(), lines.end(), "Graph Edges:");
    if (edgeStart == lines.end())
        throw std::invalid_argument("Unable to find 'Graph Edges:' line in " + filename);
    int edgeStartIndex = edgeStart - lines.begin();
    // Rcpp::Rcout << "edgeStartIndex = " << edgeStartIndex << ": " << lines[edgeStartIndex] << "\n  " << lines[edgeStartIndex+1] << std::endl;

    // If there are edges
    if (lines.size() > 3) {

        int i = edgeStartIndex+1;
	
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

	    if (edgeString.find("Lambda: ") != std::string::npos) {
                edgeString = edgeString.substr(8);
		std::istringstream iss(edgeString);
		double l;
		Rcpp::NumericVector temp;
		while (iss >> l)
		    temp.push_back(l);
		    
		lambda = wrap(Rcpp::clone(temp));
                continue;
            }

	    if (edgeString.find("Alpha: ") != std::string::npos) {
                edgeString = edgeString.substr(7);
		std::istringstream iss(edgeString);
		double a;
		iss >> a;
		Rcpp::NumericVector temp = {a};
		alpha = wrap(Rcpp::clone(temp));
                continue;
            }

	    if (edgeString.find("Penalty: ") != std::string::npos) {
                edgeString = edgeString.substr(9);
		std::istringstream iss(edgeString);
		double p;
		iss >> p;
		Rcpp::NumericVector temp = {p};
		penalty = wrap(Rcpp::clone(temp));
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

	    std::size_t pos = edgeString.find(". ");
            if (pos == std::string::npos)
                throw std::invalid_argument("Error reading graph " + filename + ": edge is not formatted correctly: " + edgeString);

            edgeString = edgeString.substr(pos+2, std::string::npos);

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
	Rcpp::_["lambda"] = lambda,
	Rcpp::_["alpha"] = alpha,
	Rcpp::_["penalty"] = penalty,
        Rcpp::_["type"] = graph_type,
        Rcpp::_["markov.blankets"] = R_NilValue,
	Rcpp::_["neighbors"] = R_NilValue,
        Rcpp::_["stabilities"] = R_NilValue
    );

    result.attr("class") = "graph";
    // result["markov.blankets"] = calculateMarkovBlankets(result);

    EdgeListGraph graph(result);

    return graph.toList();
}

//' Convert an adjacency matrix into a graph
//'
//' @param adj The adjacency matrix, p x p, with non-zero values indicating the presence of an adjacency.
//' @param nodes The names of the nodes, length p.
//' @param directed TRUE if the graph should be directed. This default is  FALSE.
//' @return A graph object representing the adjacency matrix.
//' @export
//' @examples
//' mat <- matrix(sample(c(0,1), 16, replace=TRUE), nrow=4)
//' mat <- mat + t(mat)
//' nodes <- c("X1", "X2", "X3", "X4")
//' g <- adjMat2Graph(mat, nodes)
// [[Rcpp::export]]
Rcpp::List adjMat2Graph(arma::mat adj, Rcpp::StringVector nodes, bool directed = false) {
    std::vector<std::string> nodeNames(nodes.begin(), nodes.end());

    if (adj.n_rows != adj.n_cols || adj.n_rows <= 0) {
	throw std::invalid_argument("Input adjacency matrix is invalid");
    }

    if (adj.n_rows != nodeNames.size()) {
	throw std::invalid_argument("Input node names do not match the number of variables in the adjacency matrix");
    }

    std::vector<Node> _nodes;
    for (int i = 0; i < nodeNames.size(); i++) {
	_nodes.push_back(Node(new ContinuousVariable(nodeNames[i])));
    }

    EdgeListGraph g(_nodes);

    // std::vector<std::string> edgeStrings;
    // std::vector<std::string> ambiguousTriplesStrings;

    // bool dir = Rcpp::is_true(Rcpp::all(directed));
    if (directed) {
	for (arma::uword i = 0; i < adj.n_rows; i++) {
	    for (arma::uword j = 0; j < adj.n_rows; j++) {
		if (i==j) continue;
		if (adj(i,j) != 0) {
		    g.addDirectedEdge(_nodes[j], _nodes[i]);
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

    return result;
}

//' Display a graph object as text.
//'
//' @description Display a graph object as text. This is the same format as written in ".txt" save files.
//'
//' @param graph The graph object
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- mgm(sim$data)
//' printGraph(g)
// [[Rcpp::export]]
void printGraph(const Rcpp::List& graph) {
    streamGraph(graph, Rcpp::Rcout);
}

bool EdgeListGraph::removeNode(const Node& node) {
    auto it = std::find(nodes.begin(), nodes.end(), node);
    if (it==nodes.end()) {
	return false;
    }

    for (Edge edge : edgeLists[node]) {
	Node node2 = edge.getDistalNode(node);

	if (node2 != node) {
	    auto jt = std::find(edgeLists[node2].begin(), edgeLists[node2].end(), edge);
	    if (jt != edgeLists[node2].end()) {
		edgeLists[node2].erase(jt);
	    }
	    edgesSet.erase(edge);
	}
    }

    std::unordered_set<Triple> ambigTriple(ambiguousTriples);

    for (Triple t : ambigTriple) {
	if (t.getX() == node || t.getY() == node || t.getZ() == node) {
	    ambiguousTriples.erase(t);
	}
    }

    edgeLists.erase(node);
    namesHash.erase(node.getName());
    nodes.erase(it);

    return true;
}

bool EdgeListGraph::updateNode(const Node& node) {
    auto it = std::find(nodes.begin(), nodes.end(), node);
    if (it==nodes.end()) {
	return false;
    }
    
    *it = node;

    for (Edge edge : edgeLists[node]) {
	Node node2 = edge.getDistalNode(node);

	if (node2 != node) {
	    Edge newEdge(node, node2,
			 edge.getProximalEndpoint(node),
			 edge.getDistalEndpoint(node));
	    
	    auto jt = std::find(edgeLists[node2].begin(), edgeLists[node2].end(), edge);
	    if (jt != edgeLists[node2].end()) {
		// edgeLists[node2].erase(jt);
		*jt = newEdge;
	    }
	    edgesSet.erase(edge);
	    edgesSet.insert(newEdge);
	}
    }

    std::unordered_set<Triple> ambigTriple(ambiguousTriples);

    for (Triple t : ambigTriple) {
	if (t.getX() == node || t.getY() == node || t.getZ() == node) {
	    ambiguousTriples.erase(t);
	    if (t.getX() == node) {
		ambiguousTriples.insert(Triple(node, t.getY(), t.getZ()));
	    } else if (t.getY() == node) {
		ambiguousTriples.insert(Triple(t.getX(), node, t.getZ()));
	    } else if (t.getZ() == node) {
		ambiguousTriples.insert(Triple(t.getX(), t.getY(), node));
	    }
	}
    }

    std::vector<Edge> edgeList1 = edgeLists[node];
    edgeLists.erase(node);
    edgeLists.insert(std::pair<Node, std::vector<Edge>>(node, edgeList1));
    namesHash.erase(node.getName());
    namesHash.insert(std::pair<std::string, Node>(node.getName(), node));
    // nodes.erase(it);

    return true;
}

bool EdgeListGraph::removeNodes(std::vector<Node>& newNodes) {
    bool success = true;
    for (const Node& node : newNodes) {
	success = success && removeNode(node);
    }
    return success;
}

EdgeListGraph EdgeListGraph::subgraph(std::vector<Node> nodes) {
    EdgeListGraph sub(*this);

    std::vector<Node> oldNodes(this->nodes);
    std::sort(oldNodes.begin(), oldNodes.end());
    std::sort(nodes.begin(), nodes.end());

    std::vector<Node> nodes2remove;

    std::set_difference(oldNodes.begin(), oldNodes.end(),
			nodes.begin(), nodes.end(),
			std::back_inserter(nodes2remove));

    sub.removeNodes(nodes2remove);

    return sub;
}

EdgeListGraph EdgeListGraph::getCPDAG() {
    EdgeListGraph cpdag(*this);
	    
    cpdag.reorientAllWith(ENDPOINT_TAIL);

    for (const Node& b: cpdag.getNodes()) {
	
	std::vector<Node> adjacentNodes = cpdag.getAdjacentNodes(b);

	if (adjacentNodes.size() < 2) {
	    continue;
	}

	// std::sort(adjacentNodes.begin(), adjacentNodes.end());

	ChoiceGenerator cg(adjacentNodes.size(), 2);
	std::vector<int> *choice;

	for (choice = cg.next(); choice != NULL; choice = cg.next()) {
	    const Node& a = adjacentNodes[(*choice)[0]];
	    const Node& c = adjacentNodes[(*choice)[1]];

	    // Skip triples that are shielded.
	    if (cpdag.isAdjacentTo(a, c)) {
		continue;
	    }

	    if (isDefCollider(a, b, c)) {
		cpdag.setEndpoint(a, b, ENDPOINT_ARROW);
    		cpdag.setEndpoint(c, b, ENDPOINT_ARROW);
	    }
	}
    }

    MeekRules rules;
    rules.setAggressivelyPreventCycles(true);
    rules.orientImplied(cpdag);

    cpdag.setGraphType("completed partially directed acyclic graph");

    return cpdag;
}


EdgeListGraph EdgeListGraph::getMoral() {
    EdgeListGraph moral(*this);

    moral.reorientAllWith(ENDPOINT_TAIL);

    for (const Node& b: moral.getNodes()) {
	
	std::vector<Node> adjacentNodes = moral.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            continue;
        }

	// std::sort(adjacentNodes.begin(), adjacentNodes.end());

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
	    const Node& a = adjacentNodes[(*choice)[0]];
            const Node& c = adjacentNodes[(*choice)[1]];

            // Skip triples that are shielded.
            if (moral.isAdjacentTo(a, c)) {
                continue;
            }

	    if (isDefCollider(a, b, c)) {
		moral.addUndirectedEdge(a, c);
	    }
        }
    }
    
    moral.setGraphType("undirected");

    return moral;
}


EdgeListGraph EdgeListGraph::getPAG(std::vector<Node>& latent) {

    // std::vector<Node> oldNodes(this->nodes);
    // std::sort(oldNodes.begin(), oldNodes.end());
    // std::sort(latent.begin(), latent.end());

    // std::vector<Node> nodes2keep;

    // std::set_difference(oldNodes.begin(), oldNodes.end(),
    // 			latent.begin(), latent.end(),
    // 			std::back_inserter(nodes2keep));

    std::vector<Node> observedNodes;

    for (const Node& node : this->nodes) {
	auto it = std::find(latent.begin(), latent.end(), node);
	
	if (it == latent.end()) {
	    observedNodes.push_back(node);
	} else {
	    it->setObserved(false);
	    if (!updateNode(*it))
		throw std::runtime_error("Failed to update Node " + it->getName());
	}
    }
    
    // EdgeListGraph mag = subgraph(nodes2keep);
    
    EdgeListGraph pag(observedNodes);

    for (int i = 0; i < observedNodes.size()-1; i++) {
    	for (int j = i+1; j < observedNodes.size(); j++) {
    	    // if (i==j) continue;
    	    std::vector<Node> inducingPath = GraphUtils::getInducingPath(observedNodes[i],
    									 observedNodes[j],
    									 *this);
    	    // Rcpp::Rcout << "Inducing path from " << observedNodes[i] << " to " << observedNodes[j] << ":\n  ";
    	    // for (const Node& node : inducingPath) {
    	    // 	Rcpp::Rcout << node.getName() << " ";
    	    // }
    	    // Rcpp::Rcout << std::endl;

    	    if (inducingPath.size() > 1) {
    		pag.addNondirectedEdge(observedNodes[i], observedNodes[j]);
    	    }
    	}
    }

    SepsetMap sepsets;

    // for (int i = 0; i < observedNodes.size()-1; i++) {
    // 	for (int j = i+1; j < observedNodes.size(); j++) {
    // 	    if (!pag.isAdjacentTo(observedNodes[i], observedNodes[j])) {
    // 		std::vector<Node> sepset = GraphUtils::getSepset(observedNodes[i],
    // 								 observedNodes[j],
    // 								 *this);
		
    // 		if (sepset.size() == 1 && sepset[0].isNull()) {
    // 		    Rcpp::Rcout << "No sepset found for " << observedNodes[i] << " and " << observedNodes[j] << "\n";
    // 		} else{
		
    // 		    Rcpp::Rcout << "Sepset for " << observedNodes[i] << " and " << observedNodes[j] << ":\n  ";
    // 		    for (const Node& node : sepset) {
    // 		        Rcpp::Rcout << node.getName() << " ";
    // 		    }
    // 		    Rcpp::Rcout << std::endl;
		
    // 		    sepsets.set(observedNodes[i], observedNodes[j], sepset);
    // 		}
    // 	    }
    // 	}
    // }

    std::set<NodePair> checkedSepsets;

    for (const Node& b: pag.getNodes()) {
	
	std::vector<Node> adjacentNodes = pag.getAdjacentNodes(b);

	if (adjacentNodes.size() < 2) {
	    continue;
	}

	// std::sort(adjacentNodes.begin(), adjacentNodes.end());

	ChoiceGenerator cg(adjacentNodes.size(), 2);
	std::vector<int> *choice;

	for (choice = cg.next(); choice != NULL; choice = cg.next()) {
	    const Node& a = adjacentNodes[(*choice)[0]];
	    const Node& c = adjacentNodes[(*choice)[1]];

	    // Skip triples that are shielded.
	    if (pag.isAdjacentTo(a, c)) {
		continue;
	    }

	    if (checkedSepsets.count(std::minmax(a, c))>0) continue;

	    checkedSepsets.insert(std::minmax(a, c));

	    std::vector<Node> sepset = GraphUtils::getSepset(a, c, *this);
		
	    if (sepset.size() == 1 && sepset[0].isNull()) {
		// Rcpp::Rcout << "No sepset found for " << a << " and " << c << "\n";
		continue;
	    } else{
		
		// Rcpp::Rcout << "Sepset for " << a << " and " << c << ":\n  ";
		// for (const Node& node : sepset) {
		//     Rcpp::Rcout << node.getName() << " ";
		// }
		// Rcpp::Rcout << std::endl;
		
		sepsets.set(a, c, sepset);
	    }
	}
    }
    
    SepsetProducer sepsetsset(pag, sepsets, nullptr);
    sepsetsset.setOrientRule(ORIENT_SEPSETS);

    FciOrient rules(sepsetsset);
    rules.setCompleteRuleSetUsed(true);
    rules.ruleR0(pag);
    rules.doFinalOrientation(pag);

    pag.setGraphType("partial ancestral graph");

    return pag;
}

//' Calculate the CPDAG for a given DAG
//'
//' @description Create the completed partially directed acyclic graph (CPDAG) for the input directed acyclic graph (DAG). The CPDAG represents the Markov equivalence class of the true cauasl DAG. The PC algorithms are only identifiable up to the Markov equivalence class, so assessments of causal structure recovery should be compared to the CPDAG rather than the causal DAG.
//'
//' @param graph The graph object used to generate the CPDAG. Should be the ground-truth causal DAG
//' @return The CPDAG corresponding to the input DAG
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' sim$cpdag <- cpdag(sim$graph)
//' print(sim$cpdag)
// [[Rcpp::export]]
Rcpp::List cpdag(const Rcpp::List& graph) {
    EdgeListGraph dag(graph);

    Rcpp::List result = dag.getCPDAG().toList();

    return result;
}


//' Calculate the moral graph for a given DAG
//'
//' @description Create the moral graph for the input directed acyclic graph (DAG). The moral graph is the undirected graphical model that is equivalent to the input DAG.
//'
//' @param graph The graph object used to generate the moral graph. Should be the ground-truth causal DAG
//' @return The moral graph corresponding to the input DAG
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' sim$moral <- moral(sim$graph)
//' print(sim$moral)
// [[Rcpp::export]]
Rcpp::List moral(const Rcpp::List& graph) {
    EdgeListGraph dag(graph);

    Rcpp::List result = dag.getMoral().toList();

    return result;
}


//' Calculate the undirected skeleton for a given DAG
//'
//' @description Create the skeleton graph for the input directed acyclic graph (DAG). The skeleton graph is the undirected graph that contains the same adjacencies as the input DAG.
//'
//' @param graph The graph object used to generate the skeleton graph. Should be the ground-truth causal DAG
//' @return The skeleton graph corresponding to the input DAG
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' sim$skeleton <- skeleton(sim$graph)
//' print(sim$skeleton)
// [[Rcpp::export]]
Rcpp::List skeleton(const Rcpp::List& graph) {
    EdgeListGraph skeleton(graph);

    skeleton.reorientAllWith(ENDPOINT_TAIL);

    skeleton.setGraphType("undirected");

    Rcpp::List result = skeleton.toList();    

    return result;
}


//' Calculate the PAG for a given DAG and set of latent variables
//'
//' @description Create the partial ancestral graph (PAG) for the input directed acyclic graph (DAG). The PAG represents the Markov equivalence class of the true cauasl MAG. The FCI algorithms are only identifiable up to the Markov equivalence class, so assessments of causal structure recovery should be compared to the PAG rather than the causal MAG.
//'
//' @param graph The graph object used to generate the PAG. Should be the ground-truth causal DAG
//' @param latent The names of latent (unobserved) variables in the causal DAG. The default is NULL.
//' @return The PAG corresponding to the input DAG
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' sim$pag <- pag(sim$graph)
//' print(sim$pag)
// [[Rcpp::export]]
Rcpp::List pag(const Rcpp::List& graph,
	       Rcpp::Nullable<Rcpp::StringVector> latent = R_NilValue) {
    EdgeListGraph dag(graph);

    std::vector<Node> latentNodes;

    if (latent.isNotNull()) {
	for (std::string nodeName : Rcpp::as<std::vector<std::string>>(latent)) {
	    latentNodes.push_back(Node(new ContinuousVariable(nodeName)));
	}
    }

    Rcpp::List result = dag.getPAG(latentNodes).toList();

    return result;
}


// //' Calculate the skeleton Structural Hamming Distance (SHD) between two graphs. This only counts the missing and added edges, and does not consider edge orientation
// //'
// //' @param graph1 A graph object
// //' @param graph2 A graph object
// //' @return The skeleton SHD between the two graph objects
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' g <- mgm(data.n100.p25)
// //' printGraph(g)
// // [[Rcpp::export]]
// double skeletonSHD(const Rcpp::List& graph1, const Rcpp::List& graph2) {
//     double shd = 0.0;

//     EdgeListGraph g1(graph1);
//     EdgeListGraph g2(graph2);

//     std::vector<Node> nodesG1(g1.getNodes());
//     std::vector<Node> nodesG2(g2.getNodes());
//     std::set<Node> nodeSetG1(nodesG1.begin(), nodesG1.end());
//     std::set<Node> nodeSetG2(nodesG2.begin(), nodesG2.end());
//     if (nodeSetG1 != nodeSetG2) {
// 	throw std::invalid_argument("Both graphs must be constructed over the same set of nodes.");
//     }

//     if (g1.getGraphType() != g2.getGraphType()) {
// 	Rcpp::warning("The two graphs are not of the same type (undirected, completed partially directed acyclic graph, or partial ancestral graph).");
//     }

//     for (const Edge& edge : g1.getEdges()) {
// 	if (!g2.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
// 	    shd++;
// 	}
//     }

//     for (const Edge& edge : g2.getEdges()) {
// 	if (!g1.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
// 	    shd++;
// 	}
//     }
    
//     return shd;
// }


// //' Calculate the orientation Structural Hamming Distance (SHD) between two graphs. This only counts the incorrect edge endpoints for edges present in both graphs, and does not consider differences in the graph skeleton. Each different endpoint adds 0.5 to the orientation SHD (i.e. A o-> B vs. A --> B). Thus, a completely misoriented edge adds 1 to the orientation SHD (i.e. A o-o B vs. A <-- B).
// //'
// //' @param graph1 A graph object
// //' @param graph2 A graph object
// //' @return The skeleton SHD btween the two graph objects
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' g <- mgm(data.n100.p25)
// //' printGraph(g)
// // [[Rcpp::export]]
// double orientationSHD(const Rcpp::List& graph1, const Rcpp::List& graph2) {
//     double shd = 0.0;

//     EdgeListGraph g1(graph1);
//     EdgeListGraph g2(graph2);

//     std::vector<Node> nodesG1(g1.getNodes());
//     std::vector<Node> nodesG2(g2.getNodes());
//     std::set<Node> nodeSetG1(nodesG1.begin(), nodesG1.end());
//     std::set<Node> nodeSetG2(nodesG2.begin(), nodesG2.end());
//     if (nodeSetG1 != nodeSetG2) {
// 	throw std::invalid_argument("Both graphs must be constructed over the same set of nodes.");
//     }

//     if (g1.getGraphType() != g2.getGraphType()) {
// 	throw std::invalid_argument("Both graphs must be of the same type (undirected, completed partially directed acyclic graph, or partial ancestral graph).");
//     }

//     for (Edge edge1 : g1.getEdges()) {
// 	if (g2.isAdjacentTo(edge1.getNode1(), edge1.getNode2())) {
// 	    Edge edge2 = g2.getEdge(edge1.getNode1(), edge1.getNode2());
	    
// 	    if (edge1.getDistalEndpoint(edge1.getNode1()) !=
// 		edge2.getDistalEndpoint(edge1.getNode1())) {
// 		shd += 0.5;
// 	    }

// 	    if (edge1.getProximalEndpoint(edge1.getNode1()) !=
// 		edge2.getProximalEndpoint(edge1.getNode1())) {
// 		shd += 0.5;
// 	    }
// 	}
//     }
    
//     return shd;
// }

//' Structural Hamming Distance (SHD)
//'
//' @description Calculate the Structural Hamming Distance (SHD) between two graphs.
//'
//' @param graph1 A graph object
//' @param graph2 A graph object
//' @return The SHD btween the two graph objects
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' SHD(g, cpdag(sim$graph))
// [[Rcpp::export]]
Rcpp::NumericVector SHD(const Rcpp::List& graph1, const Rcpp::List& graph2) {
    double SHD=0.0, orientSHD=0.0, skelSHD=0.0;

    EdgeListGraph g1(graph1);
    EdgeListGraph g2(graph2);

    std::vector<Node> nodesG1(g1.getNodes());
    std::vector<Node> nodesG2(g2.getNodes());
    std::set<Node> nodeSetG1(nodesG1.begin(), nodesG1.end());
    std::set<Node> nodeSetG2(nodesG2.begin(), nodesG2.end());
    if (nodeSetG1 != nodeSetG2) {
	throw std::invalid_argument("Both graphs must be constructed over the same set of nodes.");
    }

    if (g1.getGraphType() != g2.getGraphType()) {
	throw std::invalid_argument("Both graphs must be of the same type (undirected, completed partially directed acyclic graph, or partial ancestral graph).");
    }

    for (Edge edge1 : g1.getEdges()) {
	if (g2.isAdjacentTo(edge1.getNode1(), edge1.getNode2())) {
	    Edge edge2 = g2.getEdge(edge1.getNode1(), edge1.getNode2());
	    
	    if (edge1 != edge2) {
		orientSHD += 1;
		if (g1.getGraphType()=="completed partially directed acyclic graph") {
		    if (!edge1.isUndirected() && !edge2.isUndirected()) {
			orientSHD += 1;
		    }
		}
		if (g1.getGraphType()=="partial ancestral graph") {
		    if (!edge1.isNondirected() && !edge1.isNondirected()) {
			orientSHD += 1;
		    }
		}
	    }
	    
	} else {
	    skelSHD += 1;
	    if (g1.getGraphType()=="completed partially directed acyclic graph" && !edge1.isUndirected()) {
		orientSHD += 1;
	    }
	    if (g1.getGraphType()=="partial ancestral graph" && !edge1.isNondirected()) {
		orientSHD += 1;
	    }
	}
    }

    for (Edge edge2 : g2.getEdges()) {
	if (!g1.isAdjacentTo(edge2.getNode1(), edge2.getNode2())) {
	    skelSHD += 1;
	    if (g2.getGraphType()=="completed partially directed acyclic graph" && !edge2.isUndirected()) {
		orientSHD += 1;
	    }
	    if (g2.getGraphType()=="partial ancestral graph" && !edge2.isNondirected()) {
		orientSHD += 1;
	    }
	}
    }

    // orientSHD = orientationSHD(graph1, graph2);
    // skelSHD = skeletonSHD(graph1, graph2);
    SHD = orientSHD + skelSHD;
    
    return Rcpp::NumericVector::create(Rcpp::_["SHD"] = SHD);
}

//' Adjacency Precision-Recall Metrics
//'
//' @description Calculate the skeleton precision, recall, F1, and Matthew's Correlation Coefficient (MCC) between an estimated and ground truth graph.
//'
//' @param estimate An estimated graph object
//' @param groundTruth A ground truth graph object
//' @return The skeleton precision, recall, F1, and MCC, between the two graph objects
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' prMetricsAdjacency(g, cpdag(sim$graph))
// [[Rcpp::export]]
Rcpp::NumericVector prMetricsAdjacency(const Rcpp::List& estimate,
				       const Rcpp::List& groundTruth) {
    EdgeListGraph est(estimate);
    EdgeListGraph truth(groundTruth);
    
    std::vector<Node> nodesEst(est.getNodes());
    std::vector<Node> nodesTruth(truth.getNodes());
    std::set<Node> nodeSetEst(nodesEst.begin(), nodesEst.end());
    std::set<Node> nodeSetTruth(nodesTruth.begin(), nodesTruth.end());
    if (nodeSetEst != nodeSetTruth) {
	throw std::invalid_argument("Estimate and ground truth graphs must be constructed over the same set of nodes.");
    }

    double tp = 0, fp = 0, fn = 0, tn = nodeSetEst.size() * (nodeSetEst.size()-1) / 2;

    

    for (const Edge& edge : est.getEdges()) {
	if (truth.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
	    tp++;
	} else {
	    fp++;
	}
    }

    for (const Edge& edge : truth.getEdges()) {
	if (!est.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
	    fn++;
	}
    }

    tn = tn - tp - fp - fn;

    double precision = tp / ((double) (tp + fp));
    double recall = tp / ((double) (tp + fn));
    double f1 = 2 * precision * recall / (precision + recall);
    double mcc = (tp*tn - fp*fn) / std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));

    return Rcpp::NumericVector::create(
	Rcpp::_["adjPrecision"] = precision,
        Rcpp::_["adjRecall"] = recall,
        Rcpp::_["adjF1"] = f1,
        Rcpp::_["adjMCC"] = mcc);
}


//' Orientation Precision-Recall Metrics
//'
//' @description Calculate the orientation precision, recall, F1, and Matthew's Correlation Coefficient (MCC) between an estimated and ground truth graph.
//'
//' @param estimate An estimated graph object
//' @param groundTruth A ground truth graph object of the same type as the estimated graph object
//' @param groundTruthDAG A ground truth graph object containing the true causal DAG. Only necessary for calculating the or precision, recall, F1, and MCC for partial ancestral graphs (PAGs)
//' @return The orientation precision, recall, F1, and MCC, between the two graph objects
//' @export
//' @examples
//' data("train_n10000_p10")
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' prMetricsOrientation(g, cpdag(sim$graph))
// [[Rcpp::export]]
Rcpp::NumericVector prMetricsOrientation(const Rcpp::List& estimate,
					 const Rcpp::List& groundTruth,
					 const Rcpp::Nullable<Rcpp::List>& groundTruthDAG = R_NilValue) {
    EdgeListGraph est(estimate);
    EdgeListGraph truth(groundTruth);
    
    std::vector<Node> nodesEst(est.getNodes());
    std::vector<Node> nodesTruth(truth.getNodes());
    std::set<Node> nodeSetEst(nodesEst.begin(), nodesEst.end());
    std::set<Node> nodeSetTruth(nodesTruth.begin(), nodesTruth.end());
    if (nodeSetEst != nodeSetTruth) {
	throw std::invalid_argument("Estimate and ground truth graphs must be constructed over the same set of nodes.");
    }

    if (est.getGraphType() != truth.getGraphType()) {
	throw std::invalid_argument("Estimate and ground truth graphs must be of the same type (undirected, completed partially directed acyclic graph, or partial ancestral graph).");
    }

    if (est.getGraphType() == "undirected") {
	if (groundTruthDAG.isNotNull())
	    Rcpp::warning("The ground truth DAG is only needed for assessing the orientation accuracy of partial ancestral graphs.");
	
	return Rcpp::NumericVector::create(
	    Rcpp::_["orientPrecision"] = NA_REAL,
	    Rcpp::_["orientRecall"] = NA_REAL,
	    Rcpp::_["orientF1"] = NA_REAL,
	    Rcpp::_["orientMCC"] = NA_REAL);
	
    } else if (est.getGraphType() == "completed partially directed acyclic graph") {
	if (groundTruthDAG.isNotNull())
	    Rcpp::warning("The ground truth DAG is only needed for assessing the orientation accuracy of partial ancestral graphs.");
	
	double tp = 0, fp = 0, fn = 0, tn = 0;

	for (Edge edge1 : est.getEdges()) {
	    if (truth.isAdjacentTo(edge1.getNode1(), edge1.getNode2())) {
	        Edge edge2 = truth.getEdge(edge1.getNode1(), edge1.getNode2());
		
		if (edge1.getDistalEndpoint(edge1.getNode1()) ==
		    edge2.getDistalEndpoint(edge1.getNode1())) {
		    if (edge2.getDistalEndpoint(edge1.getNode1()) == ENDPOINT_ARROW) {
			tp++;
		    } else if (edge2.getDistalEndpoint(edge1.getNode1()) == ENDPOINT_TAIL) {
			tn++;
		    } else {
			throw std::runtime_error("Invalid edge endpoint for completed partially directed acyclic graphs");
		    }
		} else {
		    if (edge2.getDistalEndpoint(edge1.getNode1()) == ENDPOINT_ARROW) {
			fn++;
		    } else if (edge2.getDistalEndpoint(edge1.getNode1()) == ENDPOINT_TAIL) {
			fp++;
		    } else {
			throw std::runtime_error("Invalid edge endpoint for completed partially directed acyclic graphs");
		    }
		}

	        if (edge1.getProximalEndpoint(edge1.getNode1()) ==
		    edge2.getProximalEndpoint(edge1.getNode1())) {
		    if (edge2.getProximalEndpoint(edge1.getNode1()) == ENDPOINT_ARROW) {
			tp++;
		    } else if (edge2.getProximalEndpoint(edge1.getNode1()) == ENDPOINT_TAIL) {
			tn++;
		    } else {
			throw std::runtime_error("Invalid edge endpoint for completed partially directed acyclic graphs");
		    }
		} else {
		    if (edge2.getProximalEndpoint(edge1.getNode1()) == ENDPOINT_ARROW) {
			fn++;
		    } else if (edge2.getProximalEndpoint(edge1.getNode1()) == ENDPOINT_TAIL) {
			fp++;
		    } else {
			throw std::runtime_error("Invalid edge endpoint for completed partially directed acyclic graphs");
		    }
		}
	    }
	}

	// Rcpp::Rcout << "TP = " << tp << std::endl;
	// Rcpp::Rcout << "FP = " << fp << std::endl;
	// Rcpp::Rcout << "FN = " << fn << std::endl;
	// Rcpp::Rcout << "TN = " << tn << std::endl;
	// Rcpp::Rcout << "MCC numerator = " << (tp*tn - fp*fn) << std::endl;
	// Rcpp::Rcout << "MCC denominator = " << std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)) << std::endl;
	
	double precision = tp / ((double) (tp + fp));
	double recall = tp / ((double) (tp + fn));
	double f1 = 2 * precision * recall / (precision + recall);
	double mcc = (tp*tn - fp*fn) / std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));
	
	return Rcpp::NumericVector::create(
	    Rcpp::_["orientPrecision"] = precision,
	    Rcpp::_["orientRecall"] = recall,
	    Rcpp::_["orientF1"] = f1,
	    Rcpp::_["orientMCC"] = mcc);
	
    } else if (est.getGraphType() == "partial ancestral graph") {
	
	if (groundTruthDAG.isNull())
	    throw std::invalid_argument("The ground truth DAG is needed for assessing the orientation accuracy of partial ancestral graphs.");
	
	EdgeListGraph dag(Rcpp::as<Rcpp::List>(groundTruthDAG));

	double tp = 0, fp = 0, fn = 0, tn = 0;

	for (Edge edge1 : est.getEdges()) {
	    Node node1 = edge1.getNode1();
	    Node node2 = edge1.getNode2();
	    if (truth.isAdjacentTo(node1, node2)) {
	        Edge edge2 = truth.getEdge(node1, node2);

		if (edge1.getDistalEndpoint(node1) ==
		    edge2.getDistalEndpoint(node1)) {
		    if (edge1.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
		    	tp += 1;
		    } else if (edge1.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
		    	tn += 1;
		    } else if (edge1.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
		    	tp += 0.5;
		    	tn += 0.5;
		    } else {
		    	throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
		    }
		    // tp++;
		} else {
		    if (edge2.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
			// check if node2 is ancestor of node1
			bool ancestor = dag.isAncestorOf(node2, node1);
			if (edge1.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
			    // node2 is estimated to be an ancestor of node1
			    tn += ancestor;
			    fn += !ancestor;
			} else if (edge1.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
			    // node2 is estimated to be a descendent of node1
			    tp += !ancestor;
			    fp += ancestor;
			} else {
			    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
			}
		    } else if (edge2.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
			if (edge1.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
			    fp += 1;
			} else if (edge1.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
			    fp += 0.5;
			} else {
			    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
			}
		    } else if (edge2.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
			if (edge1.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
			    fn += 1;
			} else if (edge1.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
			    fn += 0.5;
			} else {
			    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
			}
		    } else {
			throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
		    }
		}

		
		if (edge1.getProximalEndpoint(node1) ==
		    edge2.getProximalEndpoint(node1)) {
		    if (edge1.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
		    	tp += 1;
		    } else if (edge1.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
		    	tn += 1;
		    } else if (edge1.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
		    	tp += 0.5;
		    	tn += 0.5;
		    } else {
		    	throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
		    }
		    // tp++;
		} else {
		    if (edge2.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
			// check if node1 is ancestor of node2
			bool ancestor = dag.isAncestorOf(node1, node2);
			if (edge1.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
			    // node1 is estimated to be an ancestor of node2
			    tn += ancestor;
			    fn += !ancestor;
			} else if (edge1.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
			    // node1 is estimated to be a descendent of node2
			    tp += !ancestor;
			    fp += ancestor;
			} else {
			    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
			}
		    } else if (edge2.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
			if (edge1.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
			    fp += 1;
			} else if (edge1.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
			    fp += 0.5;
			} else {
			    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
			}
		    } else if (edge2.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
			if (edge1.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
			    fn += 1;
			} else if (edge1.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
			    fn += 0.5;
			} else {
			    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
			}
		    } else {
			throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
		    }
		}
	        
	    }
	}

	// Rcpp::Rcout << "TP = " << tp << std::endl;
	// Rcpp::Rcout << "FP = " << fp << std::endl;
	// Rcpp::Rcout << "FN = " << fn << std::endl;
	// Rcpp::Rcout << "TN = " << tn << std::endl;
	// Rcpp::Rcout << "MCC numerator = " << (tp*tn - fp*fn) << std::endl;
	// Rcpp::Rcout << "MCC denominator = " << std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)) << std::endl;
	
	double precision = tp / ((double) (tp + fp));
	double recall = tp / ((double) (tp + fn));
	double f1 = 2 * precision * recall / (precision + recall);
	double mcc = (tp*tn - fp*fn) / std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));
	
	return Rcpp::NumericVector::create(
	    Rcpp::_["orientPrecision"] = precision,
	    Rcpp::_["orientRecall"] = recall,
	    Rcpp::_["orientF1"] = f1,
	    Rcpp::_["orientMCC"] = mcc);
	
    } else {
	throw std::invalid_argument("Unrecognized graph type");
    }
}


//' Causal Orientaion Precision-Recall Metrics for CPDAGs
//'
//' @description Calculate the causal orientation precision, recall, and F1 between an estimated CPDAG and ground truth graph causal DAG.
//'
//' @param estimate An estimated graph object.
//' @param groundTruthDAG A ground truth graph object of the type "directed acyclic graph".
//' @return The causal orientation precision, recall, and F1 between the two graph objects
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' prMetricsCausal(g, sim$graph)
// [[Rcpp::export]]
Rcpp::NumericVector prMetricsCausal(const Rcpp::List& estimate,
				    const Rcpp::List& groundTruthDAG) {
    EdgeListGraph est(estimate);
    EdgeListGraph truth(groundTruthDAG);
    
    std::vector<Node> nodesEst(est.getNodes());
    std::vector<Node> nodesTruth(truth.getNodes());
    std::set<Node> nodeSetEst(nodesEst.begin(), nodesEst.end());
    std::set<Node> nodeSetTruth(nodesTruth.begin(), nodesTruth.end());
    if (nodeSetEst != nodeSetTruth) {
	throw std::invalid_argument("Estimate and ground truth graphs must be constructed over the same set of nodes.");
    }

    if (est.getGraphType() != "completed partially directed acyclic graph") {
	throw std::invalid_argument("This causal orientation precision and recall is only defined for completed partially directed acyclic graphs");
    }

    if (truth.getGraphType() != "directed acyclic graph") {
	throw std::invalid_argument("Ground truth graph must be a directed acyclic graph.");
    }

    // if (est.getGraphType() != truth.getGraphType()) {
    // 	throw std::invalid_argument("Estimate and ground truth graphs must be of the same type (undirected, completed partially directed acyclic graph, or partial ancestral graph).");
    // }

    if (est.getGraphType() == "undirected") {
	// if (groundTruthDAG.isNotNull())
	//     Rcpp::warning("The ground truth DAG is only needed for assessing the orientation accuracy of partial ancestral graphs.");
	
	return Rcpp::NumericVector::create(
	    Rcpp::_["causalPrecision"] = NA_REAL,
	    Rcpp::_["causalRecall"] = NA_REAL,
	    Rcpp::_["causalF1"] = NA_REAL);
	
    } else if (est.getGraphType() == "completed partially directed acyclic graph") {
	// if (groundTruthDAG.isNotNull())
	//     Rcpp::warning("The ground truth DAG is only needed for assessing the causalation accuracy of partial ancestral graphs.");
	
	double tp = 0, fp = 0, fn = 0, tn = 0;

	for (Edge edge1 : est.getEdges()) {
	    if (truth.isAdjacentTo(edge1.getNode1(), edge1.getNode2())) {
	        Edge edge2 = truth.getEdge(edge1.getNode1(), edge1.getNode2());

		if (edge1.isUndirected()) {
		    fn++;
		} else if (edge1.isDirected()) {
		    if (edge1 == edge2) {
			tp++;
			tn++;
		    } else {
			fp++;
			fn++;
		    }
		} else {
		    throw std::runtime_error("Invalid edge endpoint for completed partially directed acyclic graphs");
		}
	    } else {
		if (edge1.isDirected()) {
		    fp++;
		}
	    }
	}

	for (Edge edge1 : truth.getEdges()) {
	    if (!est.isAdjacentTo(edge1.getNode1(), edge1.getNode2())) {
		fn++;
	    }
	}

	// Rcpp::Rcout << "TP = " << tp << std::endl;
	// Rcpp::Rcout << "FP = " << fp << std::endl;
	// Rcpp::Rcout << "FN = " << fn << std::endl;
	// Rcpp::Rcout << "TN = " << tn << std::endl;
	// Rcpp::Rcout << "MCC numerator = " << (tp*tn - fp*fn) << std::endl;
	// Rcpp::Rcout << "MCC denominator = " << std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)) << std::endl;
	
	double precision = tp / ((double) (tp + fp));
	double recall = tp / ((double) (tp + fn));
	double f1 = 2 * precision * recall / (precision + recall);
	// double mcc = (tp*tn - fp*fn) / std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));
	
	return Rcpp::NumericVector::create(
	    Rcpp::_["causalPrecision"] = precision,
	    Rcpp::_["causalRecall"] = recall,
	    Rcpp::_["causalF1"] = f1);
	
    } else if (est.getGraphType() == "partial ancestral graph") {

	return Rcpp::NumericVector::create(
	    Rcpp::_["causalPrecision"] = NA_REAL,
	    Rcpp::_["causalRecall"] = NA_REAL,
	    Rcpp::_["causalF1"] = NA_REAL);
	
	// if (groundTruthDAG.isNull())
	//     throw std::invalid_argument("The ground truth DAG is needed for assessing the orientation accuracy of partial ancestral graphs.");
	
	// EdgeListGraph dag(Rcpp::as<Rcpp::List>(groundTruthDAG));

	// double tp = 0, fp = 0, fn = 0, tn = 0;

	// for (Edge edge1 : est.getEdges()) {
	//     Node node1 = edge1.getNode1();
	//     Node node2 = edge1.getNode2();
	//     if (truth.isAdjacentTo(node1, node2)) {
	//         Edge edge2 = truth.getEdge(node1, node2);

	// 	if (edge1.getDistalEndpoint(node1) ==
	// 	    edge2.getDistalEndpoint(node1)) {
	// 	    if (edge1.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
	// 	    	tp += 1;
	// 	    } else if (edge1.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
	// 	    	tn += 1;
	// 	    } else if (edge1.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 	    	tp += 0.5;
	// 	    	tn += 0.5;
	// 	    } else {
	// 	    	throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 	    }
	// 	    // tp++;
	// 	} else {
	// 	    if (edge2.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 		// check if node2 is ancestor of node1
	// 		bool ancestor = dag.isAncestorOf(node2, node1);
	// 		if (edge1.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
	// 		    // node2 is estimated to be an ancestor of node1
	// 		    tn += ancestor;
	// 		    fn += !ancestor;
	// 		} else if (edge1.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
	// 		    // node2 is estimated to be a descendent of node1
	// 		    tp += !ancestor;
	// 		    fp += ancestor;
	// 		} else {
	// 		    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 		}
	// 	    } else if (edge2.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
	// 		if (edge1.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
	// 		    fp += 1;
	// 		} else if (edge1.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 		    fp += 0.5;
	// 		} else {
	// 		    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 		}
	// 	    } else if (edge2.getDistalEndpoint(node1) == ENDPOINT_ARROW) {
	// 		if (edge1.getDistalEndpoint(node1) == ENDPOINT_TAIL) {
	// 		    fn += 1;
	// 		} else if (edge1.getDistalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 		    fn += 0.5;
	// 		} else {
	// 		    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 		}
	// 	    } else {
	// 		throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 	    }
	// 	}

		
	// 	if (edge1.getProximalEndpoint(node1) ==
	// 	    edge2.getProximalEndpoint(node1)) {
	// 	    if (edge1.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
	// 	    	tp += 1;
	// 	    } else if (edge1.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
	// 	    	tn += 1;
	// 	    } else if (edge1.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 	    	tp += 0.5;
	// 	    	tn += 0.5;
	// 	    } else {
	// 	    	throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 	    }
	// 	    // tp++;
	// 	} else {
	// 	    if (edge2.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 		// check if node1 is ancestor of node2
	// 		bool ancestor = dag.isAncestorOf(node1, node2);
	// 		if (edge1.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
	// 		    // node1 is estimated to be an ancestor of node2
	// 		    tn += ancestor;
	// 		    fn += !ancestor;
	// 		} else if (edge1.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
	// 		    // node1 is estimated to be a descendent of node2
	// 		    tp += !ancestor;
	// 		    fp += ancestor;
	// 		} else {
	// 		    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 		}
	// 	    } else if (edge2.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
	// 		if (edge1.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
	// 		    fp += 1;
	// 		} else if (edge1.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 		    fp += 0.5;
	// 		} else {
	// 		    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 		}
	// 	    } else if (edge2.getProximalEndpoint(node1) == ENDPOINT_ARROW) {
	// 		if (edge1.getProximalEndpoint(node1) == ENDPOINT_TAIL) {
	// 		    fn += 1;
	// 		} else if (edge1.getProximalEndpoint(node1) == ENDPOINT_CIRCLE) {
	// 		    fn += 0.5;
	// 		} else {
	// 		    throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 		}
	// 	    } else {
	// 		throw std::runtime_error("Invalid edge endpoint for partial ancestral graphs");
	// 	    }
	// 	}
	        
	//     }
	// }

	// // Rcpp::Rcout << "TP = " << tp << std::endl;
	// // Rcpp::Rcout << "FP = " << fp << std::endl;
	// // Rcpp::Rcout << "FN = " << fn << std::endl;
	// // Rcpp::Rcout << "TN = " << tn << std::endl;
	// // Rcpp::Rcout << "MCC numerator = " << (tp*tn - fp*fn) << std::endl;
	// // Rcpp::Rcout << "MCC denominator = " << std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)) << std::endl;
	
	// double precision = tp / ((double) (tp + fp));
	// double recall = tp / ((double) (tp + fn));
	// double f1 = 2 * precision * recall / (precision + recall);
	// double mcc = (tp*tn - fp*fn) / std::sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));
	
	// return Rcpp::NumericVector::create(
	//     Rcpp::_["orientPrecision"] = precision,
	//     Rcpp::_["orientRecall"] = recall,
	//     Rcpp::_["orientF1"] = f1,
	//     Rcpp::_["orientMCC"] = mcc);
	
    } else {
	throw std::invalid_argument("Unrecognized graph type");
    }
}


//' Combined adjaceny and orientation precision-recall metrics
//'
//' @description Calculate the precision, recall, F1, and Matthew's Correlation Coefficient (MCC) for the adjacencies and orientations of an estimated graph compared to the ground truth. This is the concatenated output of the adjacency PR metrics and the orientation PR metrics.
//' @param estimate An estimated graph object
//' @param groundTruth A ground truth graph object of the same type as the estimated graph object
//' @param groundTruthDAG A ground truth graph object containing the true causal DAG. Only necessary for calculating the or precision, recall, F1, and MCC for partial ancestral graphs (PAGs)
//' @return The orientation precision, recall, F1, and MCC, between the two graph objects
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' prMetrics(g, cpdag(sim$graph))
// [[Rcpp::export]]
Rcpp::NumericVector prMetrics(const Rcpp::List& estimate,
			      const Rcpp::List& groundTruth,
			      const Rcpp::Nullable<Rcpp::List>& groundTruthDAG = R_NilValue) {

    Rcpp::NumericVector orientationPR = prMetricsOrientation(estimate, groundTruth,
							     groundTruthDAG);

    Rcpp::NumericVector adjPR = prMetricsAdjacency(estimate, groundTruth);
    
    return Rcpp::NumericVector::create(
	Rcpp::_["adjPrecision"] = adjPR["adjPrecision"],
	Rcpp::_["adjRecall"] = adjPR["adjRecall"],
	Rcpp::_["adjF1"] = adjPR["adjF1"],
	Rcpp::_["adjMCC"] = adjPR["adjMCC"],
	Rcpp::_["orientPrecision"] = orientationPR["orientPrecision"],
	Rcpp::_["orientRecall"] = orientationPR["orientRecall"],
	Rcpp::_["orientF1"] = orientationPR["orientF1"],
	Rcpp::_["orientMCC"] = orientationPR["orientMCC"]);
}


//' Combined graph recovery metrics
//'
//' @description Calculate the SHD, precision, recall, F1, and Matthew's Correlation Coefficient (MCC) for the adjacencies and orientations of an estimated graph compared to the ground truth. This is the concatenated output of the SHD, adjacency PR metrics, and the orientation PR metrics.
//' @param estimate An estimated graph object
//' @param groundTruth A ground truth graph object of the same type as the estimated graph object
//' @param groundTruthDAG A ground truth graph object containing the true causal DAG. Only necessary for calculating the or precision, recall, F1, and MCC for partial ancestral graphs (PAGs)
//' @return The orientation precision, recall, F1, and MCC, between the two graph objects
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' allMetrics(g, cpdag(sim$graph))
// [[Rcpp::export]]
Rcpp::NumericVector allMetrics(const Rcpp::List& estimate,
			       const Rcpp::List& groundTruth,
			       const Rcpp::Nullable<Rcpp::List>& groundTruthDAG = R_NilValue) {

    Rcpp::NumericVector pr = prMetrics(estimate, groundTruth, groundTruthDAG);

    Rcpp::NumericVector shd = SHD(estimate, groundTruth);
    
    return Rcpp::NumericVector::create(
	Rcpp::_["SHD"] = shd["SHD"],
	Rcpp::_["adjPrecision"] = pr["adjPrecision"],
	Rcpp::_["adjRecall"] = pr["adjRecall"],
	Rcpp::_["adjF1"] = pr["adjF1"],
	Rcpp::_["adjMCC"] = pr["adjMCC"],
	Rcpp::_["orientPrecision"] = pr["orientPrecision"],
	Rcpp::_["orientRecall"] = pr["orientRecall"],
	Rcpp::_["orientF1"] = pr["orientF1"],
	Rcpp::_["orientMCC"] = pr["orientMCC"]);
}
