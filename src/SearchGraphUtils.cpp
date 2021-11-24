#include "SearchGraphUtils.hpp"


/**
 * Step C of PC; orients colliders using specified sepset. That is, orients x *-* y *-* z as x *-> y <-* z just in
 * case y is in Sepset({x, z}).
 */
std::vector<Triple> SearchGraphUtils::orientCollidersUsingSepsets(SepsetMap& set,
								  EdgeListGraph& graph,
								  bool verbose) {

    if (verbose) Rcpp::Rcout << "Starting Collider Orientation" << std::endl;

    std::vector<Triple> colliders;
    std::vector<Triple> colliderList;
    // std::vector<Triple> tripleList;
    std::unordered_map<Triple, double> scores;

    // int numTriples = 0;

    std::vector<Variable*> nodes = graph.getNodes();

    for (Variable* b : nodes) {
        std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) continue;

	std::sort(adjacentNodes.begin(),
		  adjacentNodes.end(),
		  [] (Variable* a, Variable* b) {return a->getName() < b->getName(); }
	    );

        ChoiceGenerator cg(adjacentNodes.size(), 2);

        std::vector<int> *combination;

        for (combination = cg.next(); combination != NULL; combination = cg.next()) {
            Variable* a = adjacentNodes[(*combination)[0]];
            Variable* c = adjacentNodes[(*combination)[1]];

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c)) continue;

            std::vector<Variable*>* sepset = set.get(a, c);

            if (sepset != NULL) {

		// numTriples++;
		// tripleList.push_back(Triple(a, b, c));
		// scores[Triple(a, b, c)] = set.getPValue(a, c);
		
                if (std::find(sepset->begin(), sepset->end(), b) == sepset->end()) { // !sepset.contains(b)
                    // graph.setEndpoint(a, b, ENDPOINT_ARROW);
                    // graph.setEndpoint(c, b, ENDPOINT_ARROW);
		    // colliders.push_back(Triple(a, b, c));
                    colliderList.push_back(Triple(a, b, c));
		    scores[Triple(a, b, c)] = set.getPValue(a, c);
                }
            }
        }
    }

    // std::sort(tripleList.begin(), tripleList.end(),
    // 	      [&](const Triple& t1, const Triple& t2) {
    // 		  return scores[t1] < scores[t2];
    // 	      });

    // int i = 1;
    // for (Triple triple : tripleList) {
    //     Variable* a = triple.getX();
    //     Variable* b = triple.getY();
    //     Variable* c = triple.getZ();

    // 	double fdrpval = numTriples / ((double) i) * scores[triple];

    // 	if (fdrpval > 0.05) {

    // 	    std::vector<Variable*>* sepset = set.get(a, c);

    // 	    if (std::find(sepset->begin(), sepset->end(), b) == sepset->end()) {
    // 		graph.setEndpoint(a, b, ENDPOINT_ARROW);
    // 		graph.setEndpoint(c, b, ENDPOINT_ARROW);
    // 		colliders.push_back(Triple(a, b, c));
    // 	    }
    // 	}
	
    // 	i++;
    // }

    // Most independent ones first.
    std::sort(colliderList.begin(), colliderList.end(),
    	      [&](Triple& t1, Triple& t2) {
		  int s1 = set.get(t1.getX(), t1.getZ())->size();
		  int s2 = set.get(t2.getX(), t2.getZ())->size();
		  if (s1 == s2) {
		      return scores[t1] > scores[t2];
		  }
		  return s1 < s2;
    	      });

    for (Triple triple : colliderList) {
        Variable* a = triple.getX();
        Variable* b = triple.getY();
        Variable* c = triple.getZ();

        if (!(graph.getEndpoint(b, a) == ENDPOINT_ARROW ||
    	      graph.getEndpoint(b, c) == ENDPOINT_ARROW)) {
            if (orientCollider(set, graph, a, b, c)) {
    		colliders.push_back(triple);
    	    }
        }
    }

    // for (Edge edge : graph.getEdges()) {
    // 	if ((edge.getEndpoint1() == ENDPOINT_ARROW) && (edge.getEndpoint2() == ENDPOINT_ARROW)) {
    // 	    Variable* x = edge.getNode1();
    // 	    Variable* y = edge.getNode2();

    // 	    graph.setEndpoint(x, y, ENDPOINT_TAIL);
    // 	    graph.setEndpoint(y, x, ENDPOINT_TAIL);
    // 	}
    // }

    if (verbose) Rcpp::Rcout << "Finishing Collider Orientation" << std::endl;

    return colliders;
}

bool SearchGraphUtils::orientCollider(SepsetMap& set, EdgeListGraph& graph, Variable* a, Variable* b, Variable* c) {
    if (wouldCreateBadCollider(set, graph, a, b)) return false;
    if (wouldCreateBadCollider(set, graph, c, b)) return false;
    if (graph.getEdges(a, b).size() > 1) return false;
    if (graph.getEdges(b, c).size() > 1) return false;
    graph.removeEdge(a, b);
    graph.removeEdge(c, b);
    graph.addDirectedEdge(a, b);
    graph.addDirectedEdge(c, b);
    // Rcpp::Rcout << "ORIENTED SUCCESSFULLY" << std::endl;
    return true;
}

bool SearchGraphUtils::wouldCreateBadCollider(SepsetMap& set, EdgeListGraph& graph, Variable* x, Variable* y) {
    std::unordered_set<Variable*> empty = {};
    std::unordered_set<Variable*> ySet = {y};

    for (Variable* z : graph.getAdjacentNodes(y)) {
        if (x == z) continue;

        // if (!graph.isAdjacentTo(x, z) &&
	//     graph.getEndpoint(z, y) == ENDPOINT_ARROW &&
	//     !sepset(x, z, empty, ySet)) {
	//     return true;
	// }
	std::vector<Variable*>* sepset = set.get(x, z);
	if (!graph.isAdjacentTo(x, z) &&
	    graph.getEndpoint(z, y) == ENDPOINT_ARROW &&
	    std::find(sepset->begin(), sepset->end(), y) != sepset->end()) {
	    return true;
	}
    }

    return false;
}
