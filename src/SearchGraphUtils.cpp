#include "SearchGraphUtils.hpp"

/**
 * Orients edges in the graph based on prior knowledge
 */
void SearchGraphUtils::pcOrientbk(Knowledge& knowledge, EdgeListGraph& graph) {
    
    std::set<Edge> edgeSet = graph.getEdges();
    for (Edge edge : edgeSet) {
	Node node1 = edge.getNode1();
	Node node2 = edge.getNode2();

	if (knowledge.isForbidden(node1, node2) ||
	    knowledge.isRequired(node2, node1)) {
	    
	    graph.removeEdges(node1, node2);
	    
	    if (!knowledge.isForbidden(node2, node1))
		graph.addDirectedEdge(node2, node1);
	    
	} else if (knowledge.isForbidden(node2, node1) ||
		   knowledge.isRequired(node1, node2)) {
	    
	    graph.removeEdges(node1, node2);
	    
	    if (!knowledge.isForbidden(node1, node2))
		graph.addDirectedEdge(node1, node2);
	}
	
    }
}

bool SearchGraphUtils::isArrowheadAllowed(const Node& from, const Node& to,
					  Knowledge& knowledge) {
    return !knowledge.isRequired(to, from) && !knowledge.isForbidden(from, to);
}

/**
 * Step C of PC; orients colliders using specified sepset. That is, orients x *-* y *-* z as x *-> y <-* z just in
 * case y is in Sepset({x, z}).
 */
std::vector<Triple> SearchGraphUtils::orientCollidersUsingSepsets(SepsetMap& set,
								  EdgeListGraph& graph,
								  Knowledge& knowledge,
								  bool verbose) {

    // if (verbose) Rcpp::Rcout << "Starting Collider Orientation" << std::endl;

    std::vector<Triple> colliders;
    std::vector<Triple> colliderList;
    // std::vector<Triple> tripleList;
    std::unordered_map<Triple, double> scores;

    // int numTriples = 0;

    std::vector<Node> nodes = graph.getNodes();

    for (const Node& b : nodes) {
        std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) continue;

	// std::sort(adjacentNodes.begin(),
	// 	  adjacentNodes.end(),
	// 	  [] (const Node& a, const Node& b) {return a < b; }
	//     );

        ChoiceGenerator cg(adjacentNodes.size(), 2);

        std::vector<int> *combination;

        for (combination = cg.next(); combination != NULL; combination = cg.next()) {
            Node a = adjacentNodes[(*combination)[0]];
            Node c = adjacentNodes[(*combination)[1]];

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c)) continue;

	    
            std::vector<Node> sepset;
	    if (set.isInSepsetMap(a, c)) {
		sepset = set.get(a, c);
		
		if (std::find(sepset.begin(), sepset.end(), b) == sepset.end() &&
		    isArrowheadAllowed(a, b, knowledge) &&
		    isArrowheadAllowed(c, b, knowledge)) { 
                    colliderList.push_back(Triple(a, b, c));
		    scores[Triple(a, b, c)] = set.getPValue(a, c);
                }
            }
        }
    }

    // Most independent ones first.
    std::sort(colliderList.begin(), colliderList.end(),
    	      [&](Triple& t1, Triple& t2) {
		  if (scores[t1] == scores[t2])
		      return t1 < t2;
		  return scores[t1] < scores[t2];
	      });

    for (Triple triple : colliderList) {
        Node a = triple.getX();
        Node b = triple.getY();
        Node c = triple.getZ();

        if (!knowledge.isForbidden(b, a) && !knowledge.isForbidden(b, c)) {
            if (orientCollider(set, graph, a, b, c)) {
    		colliders.push_back(triple);
    	    }
        }
    }

    return colliders;
}

void SearchGraphUtils::orientCollidersUsingOrderedColliders(
    std::vector<Triple>& orderedColliders,
    EdgeListGraph& graph,
    Knowledge& knowledge,
    bool verbose) {
    
    graph.reorientAllWith(ENDPOINT_TAIL);
    pcOrientbk(knowledge, graph);

    for (Triple triple : orderedColliders) {
        Node a = triple.getX();
        Node b = triple.getY();
        Node c = triple.getZ();

        // if (!(graph.getEndpoint(b, a) == ENDPOINT_ARROW
	//       || graph.getEndpoint(b, c) == ENDPOINT_ARROW)) {
	//     orientCollider(orderedColliders, graph, a, b, c);
        // }
	
	if (!knowledge.isForbidden(b, a) && !knowledge.isForbidden(b, c)) {
	    orientCollider(orderedColliders, graph, a, b, c);
	}
    }
}


bool SearchGraphUtils::orientCollider(SepsetMap& set, EdgeListGraph& graph, const Node& a, const Node& b, const Node& c) {
    // if (wouldCreateBadCollider(set, graph, a, b)) return false;
    // if (wouldCreateBadCollider(set, graph, c, b)) return false;
    if (graph.getEdges(a, b).size() > 1) return false;
    if (graph.getEdges(b, c).size() > 1) return false;
    graph.removeEdge(a, b);
    graph.removeEdge(c, b);
    graph.addDirectedEdge(a, b);
    graph.addDirectedEdge(c, b);
    // Rcpp::Rcout << "ORIENTED SUCCESSFULLY" << std::endl;
    return true;
}

bool SearchGraphUtils::wouldCreateBadCollider(SepsetMap& set, EdgeListGraph& graph, const Node& x, const Node& y) {

    for (Node z : graph.getAdjacentNodes(y)) {
        if (x == z) continue;

	std::vector<Node> sepset;
	if (set.isInSepsetMap(x, z)) { 
	    sepset = set.get(x, z);
	} else {
	    continue;
	}	
	
	if (!graph.isAdjacentTo(x, z) &&
	    graph.getEndpoint(z, y) == ENDPOINT_ARROW &&
	    std::find(sepset.begin(), sepset.end(), y) != sepset.end()) {
	    return true;
	}
    }

    return false;
}


bool SearchGraphUtils::orientCollider(std::vector<Triple>& colliders, EdgeListGraph& graph,
				      const Node& a, const Node& b, const Node& c) {
    // if (wouldCreateBadCollider(colliders, graph, a, b)) return false;
    // if (wouldCreateBadCollider(colliders, graph, c, b)) return false;
    if (graph.getEdges(a, b).size() > 1) return false;
    if (graph.getEdges(b, c).size() > 1) return false;
    graph.removeEdge(a, b);
    graph.removeEdge(c, b);
    graph.addDirectedEdge(a, b);
    graph.addDirectedEdge(c, b);
    // Rcpp::Rcout << "ORIENTED SUCCESSFULLY" << std::endl;
    // graph.setEndpoint(a, b, ENDPOINT_ARROW);
    // graph.setEndpoint(c, b, ENDPOINT_ARROW);
    return true;
}

bool SearchGraphUtils::wouldCreateBadCollider(std::vector<Triple>& colliders,
					      EdgeListGraph& graph,
					      const Node& x, const Node& y) {

    for (Node z : graph.getAdjacentNodes(y)) {
        if (x == z) continue;

	if (!graph.isAdjacentTo(x, z) &&
	    graph.getEndpoint(z, y) == ENDPOINT_ARROW &&
	    std::find(colliders.begin(), colliders.end(), Triple(x,y,z)) != colliders.end()) {
	    return true;
	}
    }

    return false;
}
