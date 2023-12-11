#include "FciOrient.hpp"
#include <algorithm>
#include <queue>


//============================CONSTRUCTORS============================//

/**
 * Constructs a new FCI search for the given independence test and
 * background knowledge.
 */
FciOrient::FciOrient(SepsetProducer& sepsets) {
    this->sepsets = sepsets;
    this->whyOrient = std::unordered_map<std::string,std::string>();

    //Different sepsetProducers
    // if (sepsets instanceof SepsetsGreedy) {
    //     SepsetsGreedy _sepsets = (SepsetsGreedy) sepsets;
    //     this->dag = _sepsets.getDag();
    // } else if (sepsets instanceof DagSepsets) {
    //     DagSepsets _sepsets = (DagSepsets) sepsets;
    //     this->edag = _sepsets.getDag();
    // }
}

FciOrient::FciOrient(SepsetProducer& sepsets, std::unordered_map<std::string,std::string> whyOrient) {
    this->sepsets = sepsets;
    this->whyOrient = whyOrient;
//   if (sepsets instanceof SepsetsGreedy) {
//       SepsetsGreedy _sepsets = (SepsetsGreedy) sepsets;
//       this.dag = _sepsets.getDag();
//   } else if (sepsets instanceof  DagSepsets) {
//       DagSepsets _sepsets = (DagSepsets) sepsets;
//       this.dag = _sepsets.getDag();
//   }
}

//============================PUBLIC METHODS============================//

EdgeListGraph FciOrient::orient(EdgeListGraph& graph) {
    ruleR0(graph);

    // Step CI D. (Zhang's step F4.)
    doFinalOrientation(graph);

    // graph.closeInducingPaths();   //to make sure it's a legal PAG

    return graph;
}

/**
 * Orients colliders in the graph. (FCI Step C)
 * <p>
 * Zhang's step F3, rule R0.
 */
void FciOrient::ruleR0(EdgeListGraph& graph) {
    graph.reorientAllWith(ENDPOINT_CIRCLE);
    fciOrientbk(knowledge, graph);

    // Rcpp::Rcout << "Executing rule 0:\n" << graph << std::endl;

    std::vector<Node> nodes = graph.getNodes();

    for (const Node& b : nodes) {
	std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);

	if (adjacentNodes.size() < 2) {
	    continue;
	}

	ChoiceGenerator cg(adjacentNodes.size(), 2);
	std::vector<int>  *combination;

	for (combination = cg.next(); combination != NULL; combination = cg.next()) {
	    Node a = adjacentNodes.at(combination->at(0));
	    Node c = adjacentNodes.at(combination->at(1));

	    // Skip triples that are shielded or where a and c are censored.
	    if (graph.isAdjacentTo(a, c) || (a.isCensored() && c.isCensored())) {
		continue;
	    }

	    if (graph.isDefCollider(a, b, c)) {
		continue;
	    }

	    if (sepsets.isCollider(a, b, c)) {
		if (!isArrowpointAllowed(a, b, graph)) {
		    continue;
		}

		if (!isArrowpointAllowed(c, b, graph)) {
		    continue;
		}

		graph.setEndpoint(a, b, ENDPOINT_ARROW);
		graph.setEndpoint(c, b, ENDPOINT_ARROW);

		/**FOR ANALYZING CONSISTENCY**/
		// std::vector<Node> temp = sepsets.getSepset(a,c);
		whyOrient.insert(std::pair<std::string, std::string>(a.getName() + "," + b.getName(), "0"));
		whyOrient.insert(std::pair<std::string, std::string>(c.getName() + "," + b.getName(), "0"));
		// if(temp.empty()) {
		//     whyOrient.insert(std::pair<std::string, std::string>(a->getName() + "," + b->getName(), "0,null"));
		//     whyOrient.insert(std::pair<std::string, std::string>(c->getName() + "," + b->getName(), "0,null"));
		// }
		// else {
		//     std::string x = "0";
		//     for(Node p:temp) {
		//         x+=(","+p->getName());
		//     }
		//     whyOrient.insert(std::pair<std::string, std::string>(a->getName()+","+b->getName(),x));
		//     whyOrient.insert(std::pair<std::string, std::string>(c->getName()+","+b->getName(),x));
		// }
	    }
	}
    }
}

/**
 * Orients the graph according to rules in the graph (FCI step D).
 * <p>
 * Zhang's step F4, rules R1-R10.
 */
void FciOrient::doFinalOrientation(EdgeListGraph& graph) {
    if (completeRuleSetUsed) {
	// Rcpp::Rcout << "  doing complete rule set\n";
	zhangFinalOrientation(graph);
    } else {
	// Rcpp::Rcout << "  doing spirtes rule set\n";
	spirtesFinalOrientation(graph);
    }
}


//Does all 3 of these rules at once instead of going through all
// triples multiple times per iteration of doFinalOrientation.
void FciOrient::rulesR1R2cycle(EdgeListGraph& graph) {
    std::vector<Node> nodes = graph.getNodes();

    for (const Node& B : nodes) {
	std::vector<Node> adj = graph.getAdjacentNodes(B);

	if (adj.size() < 2) {
	    continue;
	}

	ChoiceGenerator cg(adj.size(), 2);
	std::vector<int> *combination;

	for (combination = cg.next(); combination != NULL; combination = cg.next()) {
	    Node A = adj.at(combination->at(0));
	    Node C = adj.at(combination->at(1));

	    //choice gen doesnt do diff orders, so must switch A & C around.
	    ruleR1(A, B, C, graph);
	    ruleR1(C, B, A, graph);
	    ruleR2(A, B, C, graph);
	    ruleR2(C, B, A, graph);
	}
    }
}


/**
 * Implements the double-triangle orientation rule, which states that if
 * D*-oB, A*->B<-*C and A*-oDo-*C, then
 * D*->B.
 * <p>
 * This is Zhang's rule R3.
 */
void FciOrient::ruleR3(EdgeListGraph& graph) {
    std::vector<Node> nodes = graph.getNodes();


    for (const Node& B : nodes) {

	std::vector<Node> intoBArrows = graph.getNodesInTo(B, ENDPOINT_ARROW);
	std::vector<Node> intoBCircles = graph.getNodesInTo(B, ENDPOINT_CIRCLE);

	for (const Node& D : intoBCircles) {
	    if (intoBArrows.size() < 2) {
		continue;
	    }

	    ChoiceGenerator gen(intoBArrows.size(), 2);
	    std::vector<int> *combination;

	    for (combination = gen.next(); combination != NULL; combination = gen.next()) {
		Node A = intoBArrows.at(combination->at(0));
		Node C = intoBArrows.at(combination->at(1));

		if (graph.isAdjacentTo(A, C)) {
		    continue;
		}

		if (!graph.isAdjacentTo(A, D)
		    || !graph.isAdjacentTo(C, D)) {
		    continue;
		}

		if (!sepsets.isNoncollider(A, D, C)) {
		    continue;
		}

		if (graph.getEndpoint(A, D) != ENDPOINT_CIRCLE) {
		    continue;
		}

		if (graph.getEndpoint(C, D) != ENDPOINT_CIRCLE) {
		    continue;
		}

		if (!isArrowpointAllowed(D, B, graph)) {
		    continue;
		}

		// check does not create cycle
		if (!doesNotCreateAlmostCycle(D, B, graph)) {
		    return;
		}

		// Rcpp::Rcout << "\n    R3: " << graph.getEdge(D, B) << " to ";

		graph.setEndpoint(D, B, ENDPOINT_ARROW);

		// Rcpp::Rcout << graph.getEdge(D, B) << "\n";

		whyOrient.insert(std::pair<std::string, std::string>(D.getName() + "," + B.getName(),"3,"+A.getName()+","+C.getName()));

		changeFlag = true;
	    }
	}
    }
}


/**
 * The triangles that must be oriented this way (won't be done by another
 * rule) all look like the ones below, where the dots are a collider path
 * from L to A with each node on the path (except L) a parent of C.
 * <pre>
 *          B
 *         xo           x is either an arrowhead or a circle
 *        /  \
 *       v    v
 * L....A --> C
 * </pre>
 * <p>
 * This is Zhang's rule R4, discriminating undirectedPaths.
 */
void FciOrient::ruleR4B(EdgeListGraph& graph) {
    if (skipDiscriminatingPathRule) {
	return;
    }

    std::vector<Node> nodes = graph.getNodes();

    for (const Node& b : nodes) {

	//potential A and C candidate pairs are only those
	// that look like this:   A<-*Bo-*C
	std::vector<Node> possA = graph.getNodesOutTo(b, ENDPOINT_ARROW);
	std::vector<Node> possC = graph.getNodesInTo(b, ENDPOINT_CIRCLE);

	for (const Node& a : possA) {
	    for (const Node& c : possC) {
		if (!graph.isParentOf(a, c)) {
		    continue;
		}

		if (graph.getEndpoint(b, c) != ENDPOINT_ARROW) {
		    continue;
		}

		ddpOrient(a, b, c, graph);
	    }
	}
    }
}

/**
 * a method to search "back from a" to find a DDP. It is called with a
 * reachability list (first consisting only of a). This is breadth-first,
 * utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. </p>
 * The body of a DDP consists of colliders that are parents of c.
 */
void FciOrient::ddpOrient(const Node& a,
			  const Node& b,
			  const Node& c,
			  EdgeListGraph& graph) {
    std::queue<Node> Q;

    std::unordered_set<Node> V;

    Node e;
    int distance = 0;

    std::unordered_map<Node, Node> previous;

    std::vector<Node> cParents = graph.getParents(c);

    Q.push(a);
    V.insert(a);
    V.insert(b);
    previous.insert(std::pair<Node,Node>(a, b));

    while (!Q.empty()) {
	Node t = Q.front();
	Q.pop();

	if (e.isNull() || e == t) {
	    e = t;
	    distance++;
	    if (distance > 0 && (distance > (maxPathLength == -1) ? 1000 : maxPathLength)) {
		return;
	    }
	}

	const std::vector<Node> nodesInTo = graph.getNodesInTo(t, ENDPOINT_ARROW);

	for (const Node& d : nodesInTo) {
	    if (std::count(V.begin(), V.end(), d) != 0) {
		continue;
	    }

	    previous.insert(std::pair<Node, Node>(d, t));
	    Node p = previous.at(t);

	    if (!graph.isDefCollider(d, t, p)) {
		continue;
	    }

	    previous.insert(std::pair<Node, Node>(d, t));

	    if (!graph.isAdjacentTo(d, c)) {
		if (doDdpOrientation(d, a, b, c, previous, graph)) {
		    return;
		}
	    }
	    if (std::count(cParents.begin(), cParents.end(), d) != 0) {
		Q.push(d);
		V.insert(d);
	    }
	}
    }
}

/**
 * Implements Zhang's rule R5, orient circle undirectedPaths: for any Ao-oB,
 * if there is an uncovered circle path u =
 * <A,C,...,D,B> such that A,D nonadjacent and B,C nonadjacent, then A---B
 * and orient every edge on u undirected.
 */
void FciOrient::ruleR5(EdgeListGraph& graph) {
    std::vector<Node> nodes = graph.getNodes();

    for (const Node& a : nodes) {
	std::vector<Node> adjacents = graph.getNodesInTo(a, ENDPOINT_CIRCLE);

	for (const Node& b : adjacents) {
	    if (!(graph.getEndpoint(a, b) == ENDPOINT_CIRCLE)) {
		continue;
	    }
	    // We know Ao-oB.

	    std::vector<std::vector<Node>> ucCirclePaths = getUcCirclePaths(a, b, graph);

	    for (const std::vector<Node>& u : ucCirclePaths) {
		if (u.size() < 3) {
		    continue;
		}

		Node c = u.at(1);
		Node d = u.at(u.size() - 2);

		if (graph.isAdjacentTo(a, d)) {
		    continue;
		}
		if (graph.isAdjacentTo(b, c)) {
		    continue;
		}
		// We know u is as required: R5 applies!


		graph.setEndpoint(a, b, ENDPOINT_TAIL);
		graph.setEndpoint(b, a, ENDPOINT_TAIL);
		orientTailPath(u, graph);
		changeFlag = true;
	    }
	}
    }
}

/**
 * Implements Zhang's rules R6 and R7, applies them over the graph once.
 * Orient single tails. R6: If A---Bo-*C then A---B--*C. R7: If A--oBo-*C
 * and A,C nonadjacent, then A--oB--*C
 */
void FciOrient::ruleR6R7(EdgeListGraph& graph) {
    std::vector<Node> nodes = graph.getNodes();

    for (const Node& b : nodes) {
	std::vector<Node> adjacents = graph.getAdjacentNodes(b);

	if (adjacents.size() < 2) {
	    continue;
	}

	ChoiceGenerator cg(adjacents.size(), 2);
	std::vector<int>  *combination;

	for (combination = cg.next(); combination != NULL; combination = cg.next()) {
	    Node a = adjacents.at(combination->at(0));
	    Node c = adjacents.at(combination->at(1));

	    if (graph.isAdjacentTo(a, c)) {
		continue;
	    }

	    if (!(graph.getEndpoint(b, a) == ENDPOINT_TAIL)) {
		continue;
	    }
	    if (!(graph.getEndpoint(c, b) == ENDPOINT_CIRCLE)) {
		continue;
	    }
	    // We know A--*Bo-*C.

	    if (graph.getEndpoint(a, b) == ENDPOINT_TAIL) {

		// We know A---Bo-*C: R6 applies!
		graph.setEndpoint(c, b, ENDPOINT_TAIL);

		changeFlag = true;
	    }

	    if (graph.getEndpoint(a, b) == ENDPOINT_CIRCLE) {
//                    if (graph.isAdjacentTo(a, c)) continue;


		// We know A--oBo-*C and A,C nonadjacent: R7 applies!
		graph.setEndpoint(c, b, ENDPOINT_TAIL);
		changeFlag = true;
	    }

	}
    }
}

/**
 * Implements Zhang's rules R8, R9, R10, applies them over the graph once.
 * Orient arrow tails. I.e., tries R8, R9, and R10 in that sequence on each
 * Ao->C in the graph.
 */
void FciOrient::rulesR8R9R10(EdgeListGraph& graph) {
    std::vector<Node> nodes = graph.getNodes();

    for (const Node& c : nodes) {
	std::vector<Node> intoCArrows = graph.getNodesInTo(c, ENDPOINT_ARROW);

	for (const Node& a : intoCArrows) {
	    if (!(graph.getEndpoint(c, a) == ENDPOINT_CIRCLE)) {
		continue;
	    }
	    // We know Ao->C.

	    // Try each of R8, R9, R10 in that order, stopping ASAP.
	    if (!ruleR8(a, c, graph)) {
		bool b = ruleR9(a, c, graph);

		if (!b) {
		    ruleR10(a, c, graph);
		}
	    }
	}
    }
}

/**
 * @param maxPathLength the maximum length of any discriminating path, or -1
 * if unlimited.
 */
void FciOrient::setMaxPathLength(int maxPathLength) {
    if (maxPathLength < -1) {
	throw std::invalid_argument("Max path length must be -1 (unlimited) or >= 0: " + std::to_string(maxPathLength));
    }

    this->maxPathLength = maxPathLength;
}

//============================PRIVATE METHODS============================//

void FciOrient::printWrongColliderMessage(const Node& a, const Node& b, const Node& c, std::string location, EdgeListGraph& graph) {
    // if (&truePag != NULL && graph.isDefCollider(a, b, c) && !truePag.isDefCollider(a, b, c)) {
    if (truePag.getNumNodes() >= 1 && graph.isDefCollider(a, b, c) && !truePag.isDefCollider(a, b, c)) {
	Rcpp::Rcout << location << ": Orienting collider by mistake: " << a << "*->" << b << "<-*" << c;
    }
}

void FciOrient::spirtesFinalOrientation(EdgeListGraph& graph) {
    changeFlag = true;
    bool firstTime = true;

    while (changeFlag) {
	changeFlag = false;
	rulesR1R2cycle(graph);
	ruleR3(graph);

	// R4 requires an arrow orientation.
	if (changeFlag || firstTime && !knowledge.isEmpty()) {
	    ruleR4B(graph);
	    firstTime = false;
	}
    }
}

void FciOrient::zhangFinalOrientation(EdgeListGraph& graph) {
    changeFlag = true;
    bool firstTime = true;

    while (changeFlag) {
	changeFlag = false;
	rulesR1R2cycle(graph);
	ruleR3(graph);

	// R4 requires an arrow orientation.
	if (changeFlag || firstTime && !knowledge.isEmpty()) {
	    ruleR4B(graph);
	    firstTime = false;
	}
    }

    if (isCompleteRuleSetUsed()) {
	// Now, by a remark on page 100 of Zhang's dissertation, we apply rule
	// R5 once.
	ruleR5(graph);

	// Now, by a further remark on page 102, we apply R6,R7 as many times
	// as possible.
	changeFlag = true;

	while (changeFlag) {
	    changeFlag = false;
	    ruleR6R7(graph);
	}

	// Finally, we apply R8-R10 as many times as possible.
	changeFlag = true;

	while (changeFlag) {
	    changeFlag = false;
	    rulesR8R9R10(graph);
	}

    }
}

/// R1, away from collider
// If a*->bo-*c and a, c not adjacent then a*->b->c
void FciOrient::ruleR1(const Node& a, const Node& b, const Node& c, EdgeListGraph& graph) {
    if (graph.isAdjacentTo(a, c)) {
	return;
    }

    if (graph.getEndpoint(a, b) == ENDPOINT_ARROW && graph.getEndpoint(c, b) == ENDPOINT_CIRCLE) {
	if (!isArrowpointAllowed(b, c, graph)) {
	    return;
	}

	// check does not create cycle
	if (!doesNotCreateAlmostCycle(b, c, graph)) {
	    return;
	}

	// Rcpp::Rcout << "\n    R1: " << graph.getEdge(b,c) << " to ";

	graph.setEndpoint(c, b, ENDPOINT_TAIL);
	graph.setEndpoint(b, c, ENDPOINT_ARROW);

	// Rcpp::Rcout << graph.getEdge(b,c) << "\n";
	
	whyOrient.insert(std::pair<std::string, std::string>(b.getName() + "," + c.getName(),"1"));
	changeFlag = true;
    }
}


//if a*-oc and either a-->b*->c or a*->b-->c, then a*->c
// This is Zhang's rule R2.
void FciOrient::ruleR2(const Node& a, const Node& b, const Node& c, EdgeListGraph& graph) {
    if ((graph.isAdjacentTo(a, c))
	&& (graph.getEndpoint(a, c) == ENDPOINT_CIRCLE)) {

	if ((graph.getEndpoint(a, b) == ENDPOINT_ARROW)
	    && (graph.getEndpoint(b, c) == ENDPOINT_ARROW) && ((graph.getEndpoint(b, a) == ENDPOINT_TAIL)
							       || (graph.getEndpoint(c, b) == ENDPOINT_TAIL))) {

	    if (!isArrowpointAllowed(a, c, graph)) {
		return;
	    }

	    // check does not create cycle
	    if (!doesNotCreateAlmostCycle(a, c, graph)) {
	    	return;
	    }

	    // Rcpp::Rcout << "\n    R2: " << graph.getEdge(a,c) << " to ";

	    graph.setEndpoint(a, c, ENDPOINT_ARROW);

	    // Rcpp::Rcout << graph.getEdge(a,c) << "\n";
	    
	    whyOrient.insert(std::pair<std::string, std::string>(a.getName() + "," + c.getName(), "2," + b.getName()));

	    changeFlag = true;
	}
    }
}

/**
 * The triangles that must be oriented this way (won't be done by another
 * rule) all look like the ones below, where the dots are a collider path
 * from L to A with each node on the path (except L) a parent of C.
 * <pre>
 *          B
 *         xo           x is either an arrowhead or a circle
 *        /  \
 *       v    v
 * L....A --> C
 * </pre>
 * <p>
 * This is Zhang's rule R4, discriminating undirectedPaths.
 */
void FciOrient::ruleR4A(EdgeListGraph& graph) {
    std::vector<Node> nodes = graph.getNodes();

    for (const Node& b : nodes) {

	//potential A and C candidate pairs are only those
	// that look like this:   A<-*Bo-*C
	std::vector<Node> possA = graph.getNodesOutTo(b, ENDPOINT_ARROW);
	std::vector<Node> possC = graph.getNodesInTo(b, ENDPOINT_CIRCLE);

	for (const Node& a : possA) {
	    for (const Node& c : possC) {
		if (!graph.isParentOf(a, c)) {
		    continue;
		}

		if (graph.getEndpoint(b, c) != ENDPOINT_ARROW) {
		    continue;
		}

		std::list<Node> reachable;
		reachable.push_back(a);
	    }
	}
    }
}

/**
 * a method to search "back from a" to find a DDP. It is called with a
 * reachability list (first consisting only of a). This is breadth-first,
 * utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. </p>
 * The body of a DDP consists of colliders that are parents of c.
 */
void FciOrient::reachablePathFind(const Node& a, const Node& b, const Node& c,
                                  std::list<Node> reachable, EdgeListGraph& graph) {
//        Map<Node, Node> next = new HashMap<Node, Node>();
// RFCI: stores the next node in the disciminating path
//        // path containing the nodes in the traiangle
//        next.put(a, b);
//        next.put(b, c);

    std::vector<Node> v = graph.getParents(c);
    std::unordered_set<Node> cParents(v.begin(), v.end());


    // Needed to avoid cycles in failure case.
    std::unordered_set<Node> visited;
    visited.insert(b);
    visited.insert(c);

    Node e = reachable.front();
    int distance = 0;

    // We don't want to include a,b,or c on the path, so they are added to
    // the "visited" set.  b and c are added explicitly here; a will be
    // added in the first while iteration.
    while (reachable.size() > 0) {
	Node x = reachable.front();
	reachable.pop_front();
	visited.insert(x);

	if (e == x) {
	    e = x;
	    distance++;
	    const int maxPathLength_ = (maxPathLength == -1) ? 1000 : maxPathLength;

	    if (distance > 0 && distance > maxPathLength_) {
		continue;
	    }
	}

	// Possible DDP path endpoints.
	std::vector<Node> pathExtensions = graph.getNodesInTo(x, ENDPOINT_ARROW);
	for (const Node& var: visited) {
	    std::remove(pathExtensions.begin(), pathExtensions.end(), var);
	}

	for (const Node& d : pathExtensions) {
	    // If d is reachable and not adjacent to c, its a DDP
	    // endpoint, so do DDP orientation. Otherwise, if d <-> c,
	    // add d to the list of reachable nodes.
	    if (!graph.isAdjacentTo(d, c)) {
		// Check whether <a, b, c> should be reoriented given
		// that d is not adjacent to c; if so, orient and stop.
		doDdpOrientation(d, a, b, c, graph);
		return;
	    } else if (cParents.count(d)) {
		if (graph.getEndpoint(x, d) == ENDPOINT_ARROW) {
		    reachable.push_back(d);
		}
	    }
	}
    }
}

/**
 * Orients the edges inside the definte discriminating path triangle. Takes
 * the left endpoint, and a,b,c as arguments.
 */
void FciOrient::doDdpOrientation(const Node& d, const Node& a, const Node& b, const Node& c, EdgeListGraph& graph) {

    std::vector<Node> sepset = getSepset(d, c);

    if (sepset.empty()) {
	return;
    }

    if (std::count(sepset.begin(), sepset.end(), b) != 0) {
	graph.setEndpoint(c, b, ENDPOINT_TAIL);
	changeFlag = true;

    } else {
	if (!isArrowpointAllowed(a, b, graph)) {
	    return;
	}

	if (!isArrowpointAllowed(c, b, graph)) {
	    return;
	}

	// check does not create cycle
	if (!doesNotCreateAlmostCycle(a, b, graph) || !doesNotCreateAlmostCycle(c, b, graph)) {
	    return;
	}

	graph.setEndpoint(a, b, ENDPOINT_ARROW);
	graph.setEndpoint(c, b, ENDPOINT_ARROW);
	changeFlag = true;
    }
}

/**
 * Orients the edges inside the definte discriminating path triangle. Takes
 * the left endpoint, and a,b,c as arguments.
 */
bool FciOrient::doDdpOrientation(const Node& d, const Node& a, const Node& b, const Node& c,
                                 std::unordered_map<Node, Node> previous, EdgeListGraph& graph) {

    if (dag.getNumNodes() != 0) {
	if (dag.isAncestorOf(b, c)) {
	    graph.setEndpoint(c, b, ENDPOINT_TAIL);
	    whyOrient.insert(std::pair<std::string, std::string>(c.getName() + "," + b.getName(),"4"));
	    changeFlag = true;
	} else {
	    if (!isArrowpointAllowed(a, b, graph)) {
		return false;
	    }

	    if (!isArrowpointAllowed(c, b, graph)) {
		return false;
	    }

	    // check does not create cycle
	    if (!doesNotCreateAlmostCycle(a, b, graph) || !doesNotCreateAlmostCycle(c, b, graph)) {
	    	return false;
	    }

	    graph.setEndpoint(a, b, ENDPOINT_ARROW);
	    graph.setEndpoint(c, b, ENDPOINT_ARROW);
	    whyOrient.insert(std::pair<std::string, std::string>(a.getName() + "," + b.getName(),"4"));
	    whyOrient.insert(std::pair<std::string, std::string>(c.getName() + "," + b.getName(),"4"));

	    changeFlag = true;
	}

	return true;
    }

    if (graph.isAdjacentTo(d, c)) {
	throw std::invalid_argument("Illegal Argument");
    }

    std::vector<Node> path = getPath(d, previous);

    bool ind = getSepsets().isIndependent(d, c, path);

    // std::vector<Node> path2 (path);
    // path2.remove(b);
    std::vector<Node> path2(path);
    if (std::find(path2.begin(), path2.end(), b) != path2.end()) {
	path2.erase(std::find(path2.begin(), path2.end(), b));
    }
    for(const Node& var: path) {
	if (var != b) {
	    path2.push_back(var);
	}
    }


    bool ind2 = getSepsets().isIndependent(d, c, path2);

    if (!ind && !ind2) {
	std::vector<Node> sepset = getSepsets().getSepset(d, c);

	if (sepset.empty()) {
	    return false;
	}

	ind = (std::count(sepset.begin(), sepset.end(), b));
    }

    //        printDdp(d, path, a, b, c, graph);
    if (ind) {
	//            if (sepset.contains(b)) {
	graph.setEndpoint(c, b, ENDPOINT_TAIL);

	changeFlag = true;
	return true;
    } else {
	if (!isArrowpointAllowed(a, b, graph)) {
	    return false;
	}

	if (!isArrowpointAllowed(c, b, graph)) {
	    return false;
	}

	// check does not create cycle
	if (!doesNotCreateAlmostCycle(a, b, graph) || !doesNotCreateAlmostCycle(c, b, graph)) {
	    return false;
	}

	graph.setEndpoint(a, b, ENDPOINT_ARROW);
	graph.setEndpoint(c, b, ENDPOINT_ARROW);

	changeFlag = true;
	return true;
    }


}

// graph is only used in the deleted verbose section, making it unnecessary here.
void FciOrient::printDdp(const Node& d, std::vector<Node> path,
			 const Node& a, const Node& b, const Node& c,
			 EdgeListGraph& graph) {
    std::vector<Node> nodes;
    nodes.push_back(d);
    for(const Node& var: path) {
	nodes.push_back(var);
    }
    nodes.push_back(a);
    nodes.push_back(b);
    nodes.push_back(c);
}


std::vector<Node> FciOrient::getPath(const Node& c,
				     std::unordered_map<Node, Node> previous) {
    std::vector<Node> l;

    Node p = c;

    do {
	p = previous.at(p);

	if (!p.isNull()) {
	    l.push_back(p);
	    // Rcpp::Rcout << "The loop that never ends\n";
	} else {
	    // Rcpp::Rcout << "The loop actually ends\n";
	}
    } while (!p.isNull());

    return l;
}

/**
 * Orients every edge on a path as undirected (i.e. A---B).
 * <p>
 * DOES NOT CHECK IF SUCH EDGES ACTUALLY EXIST: MAY DO WEIRD THINGS IF
 * PASSED AN ARBITRARY LIST OF NODES THAT IS NOT A PATH.
 *
 * @param path The path to orient as all tails.
 */
void FciOrient::orientTailPath(std::vector<Node> path, EdgeListGraph& graph) {
    for (int i = 0; i < path.size() - 1; i++) {
	Node n1 = path.at(i);
	Node n2 = path.at(i + 1);

	graph.setEndpoint(n1, n2, ENDPOINT_TAIL);
	graph.setEndpoint(n2, n1, ENDPOINT_TAIL);
	changeFlag = true;
    }
}

/**
 * Gets a list of every uncovered partially directed path between two nodes
 * in the graph.
 * <p>
 * Probably extremely slow.
 *
 * @param n1 The beginning node of the undirectedPaths.
 * @param n2 The ending node of the undirectedPaths.
 * @return A list of uncovered partially directed undirectedPaths from n1 to
 * n2.
 */
std::vector<std::vector<Node>> FciOrient::getUcPdPaths(const Node& n1, const Node& n2, EdgeListGraph& graph) {
    std::vector<std::vector<Node>> ucPdPaths;

    std::vector<Node> soFar;

    soFar.push_back(n1);

    std::vector<Node> adjacencies = graph.getAdjacentNodes(n1);
    for (const Node& curr : adjacencies) {
	getUcPdPsHelper(curr, soFar, n2, ucPdPaths, graph);
    }

    return ucPdPaths;
}

/**
 * Used in getUcPdPaths(n1,n2) to perform a breadth-first search on the
 * graph.
 * <p>
 * ASSUMES soFar CONTAINS AT LEAST ONE NODE!
 * <p>
 * Probably extremely slow.
 *
 * @param curr The getModel node to test for addition.
 * @param soFar The getModel partially built-up path.
 * @param end The node to finish the undirectedPaths at.
 * @param ucPdPaths The getModel list of uncovered p.d. undirectedPaths.
 */
void FciOrient::getUcPdPsHelper(const Node& curr, std::vector<Node> soFar, const Node& end,
                                std::vector<std::vector<Node>> ucPdPaths, EdgeListGraph& graph) {

    if (std::count(soFar.begin(), soFar.end(), curr) != 0) {
	return;
    }

    Node prev = soFar.at(soFar.size() - 1);
    if (graph.getEndpoint(prev, curr) == ENDPOINT_TAIL
	|| graph.getEndpoint(curr, prev) == ENDPOINT_ARROW) {
	return; // Adding curr would make soFar not p.d.
    } else if (soFar.size() >= 2) {
	Node prev2 = soFar.at(soFar.size() - 2);
	if (graph.isAdjacentTo(prev2, curr)) {
	    return; // Adding curr would make soFar not uncovered.
	}
    }

    soFar.push_back(curr); // Adding curr is OK, so let's do it.

    if (curr == end) {
	// We've reached the goal! Save soFar as a path.
	std::vector<Node> soFar_copy(soFar);
	ucPdPaths.push_back(soFar_copy);
    } else {
	// Otherwise, try each node adjacent to the getModel one.
	std::vector<Node> adjacents = graph.getAdjacentNodes(curr);
	for (const Node& next : adjacents) {
	    getUcPdPsHelper(next, soFar, end, ucPdPaths, graph);
	}
    }
    soFar.pop_back();
}


/**
 * Gets a list of every uncovered circle path between two nodes in the graph
 * by iterating through the uncovered partially directed undirectedPaths and
 * only keeping the circle undirectedPaths.
 * <p>
 * Probably extremely slow.
 *
 * @param n1 The beginning node of the undirectedPaths.
 * @param n2 The ending node of the undirectedPaths.
 * @return A list of uncovered circle undirectedPaths between n1 and n2.
 */
std::vector<std::vector<Node>> FciOrient::getUcCirclePaths(const Node& n1, const Node& n2, EdgeListGraph& graph) {
    std::vector<std::vector<Node>> ucCirclePaths;
    std::vector<std::vector<Node>> ucPdPaths = getUcPdPaths(n1, n2, graph);

    for (const std::vector<Node>& path : ucPdPaths) {
	for (int i = 0; i < path.size() - 1; i++) {
	    Node j = path.at(i);
	    Node sj = path.at(i + 1);

	    if (!(graph.getEndpoint(j, sj) == ENDPOINT_CIRCLE)) {
		break;
	    }
	    if (!(graph.getEndpoint(sj, j) == ENDPOINT_CIRCLE)) {
		break;
	    }
	    // This edge is OK, it's all circles.

	    if (i == path.size() - 2) {
		// We're at the last edge, so this is a circle path.
		ucCirclePaths.push_back(path);
	    }
	}
    }

    return ucCirclePaths;
}

/**
 * Tries to apply Zhang's rule R8 to a pair of nodes A and C which are
 * assumed to be such that Ao->C.
 * <p>
 * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
 * <p>
 * R8: If Ao->C and A-->B-->C or A--oB-->C, then A-->C.
 *
 * @param a The node A.
 * @param c The node C.
 * @return Whether or not R8 was successfully applied.
 */
bool FciOrient::ruleR8(const Node& a, const Node& c, EdgeListGraph& graph) {

    std::vector<Node> intoCArrows = graph.getNodesInTo(c, ENDPOINT_ARROW);

    for (const Node& b : intoCArrows) {
	// We have B*->C.
	if (!graph.isAdjacentTo(a, b)) {
	    continue;
	}
	if (!graph.isAdjacentTo(b, c)) {
	    continue;
	}

	// We have A*-*B*->C.
	if (!(graph.getEndpoint(b, a) == ENDPOINT_TAIL)) {
	    continue;
	}
	if (!(graph.getEndpoint(c, b) == ENDPOINT_TAIL)) {
	    continue;
	}
	// We have A--*B-->C.

	if (graph.getEndpoint(a, b) == ENDPOINT_TAIL) {
	    continue;
	}
	// We have A-->B-->C or A--oB-->C: R8 applies!


	graph.setEndpoint(c, a, ENDPOINT_TAIL);
	changeFlag = true;
	return true;
    }

    return false;
}

/**
 * Tries to apply Zhang's rule R9 to a pair of nodes A and C which are
 * assumed to be such that Ao->C.
 * <p>
 * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
 * <p>
 * R9: If Ao->C and there is an uncovered p.d. path u=<A,B,..,C> such that
 * C,B nonadjacent, then A-->C.
 *
 * @param a The node A.
 * @param c The node C.
 * @return Whether or not R9 was succesfully applied.
 */
bool FciOrient::ruleR9(const Node& a, const Node& c, EdgeListGraph& graph) {
    std::vector<std::vector<Node>> ucPdPsToC = getUcPdPaths(a, c, graph);

    for (const std::vector<Node>& u : ucPdPsToC) {
	Node b = u.at(1);
	if (graph.isAdjacentTo(b, c)) {
	    continue;
	}
	if (b == c) {
	    continue;
	}
	// We know u is as required: R9 applies!


	graph.setEndpoint(c, a, ENDPOINT_TAIL);
	changeFlag = true;
	return true;
    }

    return false;
}

/**
 * Tries to apply Zhang's rule R10 to a pair of nodes A and C which are
 * assumed to be such that Ao->C.
 * <p>
 * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
 * <p>
 * R10: If Ao->C, B-->C<--D, there is an uncovered p.d. path u1=<A,M,...,B>
 * and an uncovered p.d. path u2=
 * <A,N,...,D> with M != N and M,N nonadjacent then A-->C.
 *
 * @param a The node A.
 * @param c The node C.
 * @return Whether or not R10 was successfully applied.
 */
bool FciOrient::ruleR10(const Node& a, const Node& c, EdgeListGraph& graph) {
    std::vector<Node> intoCArrows = graph.getNodesInTo(c, ENDPOINT_ARROW);

    for (const Node& b : intoCArrows) {
	if (b == a) {
	    continue;
	}

	if (!(graph.getEndpoint(c, b) == ENDPOINT_TAIL)) {
	    continue;
	}
	// We know Ao->C and B-->C.

	for (const Node& d : intoCArrows) {
	    if (d == a || d == b) {
		continue;
	    }

	    if (!(graph.getEndpoint(d, c) == ENDPOINT_TAIL)) {
		continue;
	    }
	    // We know Ao->C and B-->C<--D.

	    std::vector<std::vector<Node>> ucPdPsToB = getUcPdPaths(a, b, graph);
	    std::vector<std::vector<Node>> ucPdPsToD = getUcPdPaths(a, d, graph);
	    for (const std::vector<Node>& u1 : ucPdPsToB) {
		Node m = u1.at(1);
		for (const std::vector<Node>& u2 : ucPdPsToD) {
		    Node n = u2.at(1);

		    if (m == n) {
			continue;
		    }
		    if (graph.isAdjacentTo(m, n)) {
			continue;
		    }
		    // We know B,D,u1,u2 as required: R10 applies!


		    graph.setEndpoint(c, a, ENDPOINT_TAIL);
		    changeFlag = true;
		    return true;
		}
	    }
	}
    }

    return false;
}

void FciOrient::fciOrientbk(Knowledge& knowledge, EdgeListGraph& graph) {
    
    std::set<Edge> edgeSet = graph.getEdges();
    for (Edge edge : edgeSet) {
	Node node1 = edge.getNode1();
	Node node2 = edge.getNode2();

	// If n1 --> n2 forbidden, then set n1 <-* n2
	if (knowledge.isForbidden(node1, node2)) {
	    if (isArrowpointAllowed(node2, node1, graph)) {
		graph.setEndpoint(node2, node1, ENDPOINT_ARROW);
	    }
	}

	// If n2 --> n1 forbidden, then set n2 <-* n1
	if (knowledge.isForbidden(node2, node1)) {
	    if (isArrowpointAllowed(node1, node2, graph)) {
		graph.setEndpoint(node1, node2, ENDPOINT_ARROW);
	    }
	}

	// If n1 --> n2 is required, then set n1 --> n2
	if (knowledge.isRequired(node1, node2)) {
	    if (isArrowpointAllowed(node1, node2, graph)) {
		graph.setEndpoint(node1, node2, ENDPOINT_ARROW);
		graph.setEndpoint(node2, node1, ENDPOINT_TAIL);
	    }
	}

	// If n2 --> n1 is required, then set n2 --> n1
	if (knowledge.isRequired(node2, node1)) {
	    if (isArrowpointAllowed(node2, node1, graph)) {
		graph.setEndpoint(node2, node1, ENDPOINT_ARROW);
		graph.setEndpoint(node1, node2, ENDPOINT_TAIL);
	    }
	}
    }
}


/**
 * Helper method. Appears to check if an arrowpoint is permitted by
 * background knowledge.
 *
 * @param x The possible other node.
 * @param y The possible point node.
 * @return Whether the arrowpoint is allowed.
 */
bool FciOrient::isArrowpointAllowed(const Node& x, const Node& y, EdgeListGraph& graph) {
    if (!graph.isAdjacentTo(x, y)) return false;
    
    if (graph.getEndpoint(x, y) == ENDPOINT_ARROW) {
	return true;
    }

    if (graph.getEndpoint(x, y) == ENDPOINT_TAIL) {
	return false;
    }

    if (graph.getEndpoint(y, x) == ENDPOINT_ARROW &&
	graph.getEndpoint(x, y) == ENDPOINT_CIRCLE) {
	// switched the commented sections
	// return true;
	if (knowledge.isForbidden(x, y)) {
	    return true;
	}
    }

    if (graph.getEndpoint(y, x) == ENDPOINT_TAIL &&
	graph.getEndpoint(x, y) == ENDPOINT_CIRCLE) {
	// added return true as did above in place of knowledge
	// return true;
	if (knowledge.isForbidden(x, y)) {
	    return false;
	}

	if (!doesNotCreateAlmostCycle(x, y, graph)) {
	    return false;
	}
    }

    return graph.getEndpoint(x, y) == ENDPOINT_CIRCLE;
}


bool FciOrient::doesNotCreateCycle(Node from, Node to, EdgeListGraph& graph) {
    if (aggressivelyPreventCycles) {
	return !graph.isAncestorOf(to, from);
    } else {
	return true;
    }
}

bool FciOrient::doesNotCreateAlmostCycle(Node from, Node to, EdgeListGraph& graph) {
    if (aggressivelyPreventCycles) {
	return !graph.existsAlmostDirectedPathFromTo(to, from);
    } else {
	return true;
    }
}
