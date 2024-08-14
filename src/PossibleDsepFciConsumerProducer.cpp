#include "PossibleDsepFciConsumerProducer.hpp"

//============================CONSTRUCTORS============================//

/**
 * Creates a new SepSet and assumes that none of the variables have yet been checked.
 *
 * @param graph The GaSearchGraph on which to work
 * @param test  The IndependenceChecker to use as an oracle
 */
PossibleDsepFciConsumerProducer::PossibleDsepFciConsumerProducer(EdgeListGraph& graph, IndependenceTest *test, int threads) : PossibleDsepFciConsumerProducer(test, threads)
{
    // Unable to compare EdgeListGraph to NULL
    // if (graph == NULL) {
    //     throw std::invalid_argument("null GaSearchGraph passed in PossibleDSepSearch constructor!");
    // }
    if (test == NULL) {
        throw std::invalid_argument("null IndependenceChecker passed in PossibleDSepSearch constructor!");
    }
    this->graph = graph;
}

PossibleDsepFciConsumerProducer::PossibleDsepFciConsumerProducer(IndependenceTest *test, int threads) : taskQueue(MAX_QUEUE_SIZE)
{
    this->test = test;
    this->sepset =  SepsetMap();

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }

    setMaxPathLength(maxReachablePathLength);
    // this->poisonEdge = Edge(new ContinuousVariable("EJWMX3RCpPi0qbp"),
    // 			    new ContinuousVariable("nLtWU7DmeZyYPZs"),
    // 			    ENDPOINT_NULL, ENDPOINT_NULL);
}

//========================PRIVATE METHODS==========================//


//This method is based on background knowledge aka IKnowledge which is not yet implemented
/**
 * Removes from the list of nodes any that cannot be parents of x given the background knowledge.
 */
std::vector<Node> PossibleDsepFciConsumerProducer::possibleParents(const Node& x, std::vector<Node>& nodes) {

    std::vector<Node> possibleParents;

    for (const Node& z : nodes) {

        if (possibleParentOf(z, x)) {
            possibleParents.push_back(z);
        }
    }
    return possibleParents;
}

bool PossibleDsepFciConsumerProducer::possibleParentOf(const Node& x, const Node& z) {
    return !knowledge.isForbidden(z, x) && !knowledge.isRequired(x, z);
}

/**
 * A variable v is in Possible-D-Sep(A,B) iff
 * <pre>
 * 	(i) v != A & v != B
 * 	(ii) there is an undirected path U between A and v such that for every
 * 		 subpath <X,Y,Z> of U either:
 * 		(a) Y is a collider on the subpath, or
 * 		(b) X is adjacent to Z.
 * </pre>
 */
std::unordered_set<Node> PossibleDsepFciConsumerProducer::getPossibleDsep(const Node& node1, const Node& node2, int maxPathLength) {

    std::unordered_set<Node>  dsep = GraphUtils::possibleDsep(node1, node2, graph, maxPathLength);

    dsep.erase(node1);
    dsep.erase(node2);

    return dsep;
}

//========================PUBLIC METHODS==========================//

/**
 * Performs pairwise comparisons of each variable in the graph with the variables that have not already been
 * checked. We get the Possible-D-Sep sets for the pair of variables, and we check to see if they are independent
 * conditional on some subset of the union of Possible-D-Sep sets. This method returns the SepSet passed in the
 * constructor (if any), possibly augmented by some edge removals in this step. The GaSearchGraph passed in the
 * constructor is directly changed.
 */
SepsetMap& PossibleDsepFciConsumerProducer::search() {
    std::unordered_map<Edge, std::vector<Node>> edgeCondsetMap;

    concurrentSearch(graph, edgeCondsetMap);

    for (const std::pair<Edge, std::vector<Node>>& entry : edgeCondsetMap) {
        Edge edge = entry.first;
        std::vector<Node> condSet = entry.second;
        Node x = edge.getNode1();
        Node y = edge.getNode2();

        graph.removeEdge(x, y);
        if (verbose) Rcpp::Rcout << "      Removed " << x.getName() << " --- " << y.getName() << "\n";
        sepset.set(x, y, condSet);
    }
    return sepset;
}

void PossibleDsepFciConsumerProducer::concurrentSearch(EdgeListGraph& graph, std::unordered_map<Edge, std::vector<Node>>& edgeCondsetMap) {
    const std::set<Edge> edges(graph.getEdges());
    std::vector<RcppThread::Thread> threads;

    threads.push_back(RcppThread::Thread( [&] { PossibleDsepProducer(edges); } ));


    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( [&] { PossibleDsepConsumer(edgeCondsetMap); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }
    return;
}

void PossibleDsepFciConsumerProducer::setDepth(int depth) {
    if (depth < -1) {
        throw  std::invalid_argument("Depth must be -1 (unlimited) or >= 0");
    }

    this->depth = depth;
}

void PossibleDsepFciConsumerProducer::setMaxPathLength(int maxReachablePathLength) {
    if (maxReachablePathLength < -1) {
	throw std::invalid_argument("Max path length must be -1 (unlimited) or >= 0");
    }
    this->maxReachablePathLength = maxReachablePathLength == std::numeric_limits<int>::max() ? -1 : maxReachablePathLength;
}

void PossibleDsepFciConsumerProducer::PossibleDsepProducer(std::set<Edge> edges) {
    // PossibleDsepTask poisonPill(Edge(), std::vector<Node>());

    int maxPossDsepSize = 0;

    bool printFlag = false;

    for (const Edge& edge : edges) {
        Node x = edge.getNode1();
        Node y = edge.getNode2();

        // std::unordered_set<Node> possibleDsepSet = getPossibleDsep(x, y, maxReachablePathLength);
	std::set<Node> possibleDsepSet = GraphUtils::possibleDsep2(x, y, graph, maxReachablePathLength);
        std::vector<Node> possibleDsep;

        possibleDsep.insert(possibleDsep.end(), possibleDsepSet.begin(), possibleDsepSet.end());
	std::sort(possibleDsep.begin(), possibleDsep.end());

	std::vector<Node> adjx = graph.getAdjacentNodes(x);
	adjx.erase(std::remove(adjx.begin(), adjx.end(), y), adjx.end());


	// RcppThread::Rcout << ss.str();

        bool noEdgeRequired = knowledge.noEdgeRequired(x, y);
        // bool noEdgeRequired = true;

        if (!noEdgeRequired)
            continue;

	int depth_ = getDepth() == -1 ? 1000 : getDepth();
	std::vector<Node> possParents;

	if (adjx != possibleDsep) {

	    if (possibleDsep.size() > maxPossDsepSize) {
		maxPossDsepSize = possibleDsep.size();
		if (verbose) {
		    printFlag = true;
		    RcppThread::Rcout << "\r      Largest Encountered Possible-D-Sep: "
				      << maxPossDsepSize;
		}
	    }

	    // RcppThread::Rcout << "      Possible-D-Sep(" << x << "," << y << ")  :  { "; // != Adj(" << x <<") \\ { " << y <<" }\n";
	
	    // // RcppThread::Rcout << "    Edge: " << edge << "  :  { ";
	    // for (Node node : possibleDsep) {
	    // 	RcppThread::Rcout << node << " ";
	    // }
	    // RcppThread::Rcout << "}\n";

	    // std::set<Node> adjSet(adjx.begin(), adjx.end());

	    //possible parents is not fully implemented without background knowledge
	    possParents = possibleParents(x, possibleDsep);

	    for (int d = 0; d <= std::min((std::size_t) depth_, possParents.size()); d++) {
		ChoiceGenerator cg(possParents.size(), d);
		std::vector<int> *combination;

		for (combination = cg.next(); combination != NULL; combination = cg.next()) {
		    std::vector<Node> condSet = GraphUtils::asList(*combination, possParents);
		    // bool allAdj = true;
		    // for (const Node& n : condSet) {
		    // 	if (adjSet.count(n)==0) {
		    // 	    allAdj = false;
		    // 	    break;
		    // 	}
		    // }
		    
		    // if (!allAdj) {
		    PossibleDsepTask newTask(edge, condSet);
		    taskQueue.push(newTask);
		    // }
		}
	    }
	}

        // possibleDsepSet = getPossibleDsep(y, x, maxReachablePathLength);	
	possibleDsepSet = GraphUtils::possibleDsep2(y, x, graph, maxReachablePathLength);
        possibleDsep.clear();
        possibleDsep.insert(possibleDsep.end(), possibleDsepSet.begin(), possibleDsepSet.end());

	std::sort(possibleDsep.begin(), possibleDsep.end());

	std::vector<Node> adjy = graph.getAdjacentNodes(y);
	adjy.erase(std::remove(adjy.begin(), adjy.end(), x), adjy.end());

	if (adjy != possibleDsep) {

	    if (possibleDsep.size() > maxPossDsepSize) {
	        maxPossDsepSize = possibleDsep.size();
		if (verbose) {
		    printFlag = true;
		    RcppThread::Rcout << "\r      Largest Encountered Possible-D-Sep: "
				      << maxPossDsepSize;
		}
	    }

	    // RcppThread::Rcout << "      Possible-D-Sep(" << y << "," << x << ")  :  { "; // != Adj(" << y <<") \\ { " << x <<" }\n";
	
	    // // RcppThread::Rcout << "    Edge: " << edge << "  :  { ";
	    // for (Node node : possibleDsep) {
	    // 	RcppThread::Rcout << node << " ";
	    // }
	    // RcppThread::Rcout << "}\n";
	    
	    // std::set<Node> adjSet(adjy.begin(), adjy.end());

	    //possible parents is not fully implemented without background knowledge
	    possParents = possibleParents(y, possibleDsep);

	    for (int d = 0; d <= std::min((std::size_t) depth_, possParents.size()); d++) {
		ChoiceGenerator cg (possParents.size(), d);
		std::vector<int> *combination;

		for (combination = cg.next(); combination != NULL; combination = cg.next()) {
		    std::vector<Node> condSet = GraphUtils::asList(*combination, possParents);
		    // bool allAdj = true;
		    // for (const Node& n : condSet) {
		    // 	if (adjSet.count(n)==0) {
		    // 	    allAdj = false;
		    // 	    break;
		    // 	}
		    // }
		    
		    // if (!allAdj) {
		    PossibleDsepTask newTask(edge, condSet);
		    taskQueue.push(newTask);
		    // }
		}
	    }
	}

	// RcppThread::Rcout << "      Largest Possible-D-Sep: " << maxPossDsepSize << std::endl;

	if (RcppThread::isInterrupted()) {
	    break;
	}
    }

    if (verbose && printFlag) {
	RcppThread::Rcout << std::endl;
    }

    for (int i = 0; i < parallelism; i++) {
	taskQueue.push(PossibleDsepTask());
    }
    return;
}

void PossibleDsepFciConsumerProducer::PossibleDsepConsumer(std::unordered_map<Edge, std::vector<Node>>& edgeCondsetMap) {
    PossibleDsepTask task = taskQueue.pop();
    while (!task.edge.isNull()) {
	
        if (edgeCondsetMap.count(task.edge) == 0) {
	    Node x = task.edge.getNode1();
	    Node y = task.edge.getNode2();
	    // if (x > y) {
	    // 	x = task.edge.getNode2();
	    // 	y = task.edge.getNode1();
	    // }
	    
            if (test->isIndependent(x, y, task.condSet)) {
                std::lock_guard<std::mutex> edgeLock(edgeMutex);
		edgeCondsetMap[task.edge] = task.condSet;
                // edgeCondsetMap.insert(std::pair<Edge,
		// 		      std::vector<Node>>(task.edge, task.condSet));
	    }
        }
	
        task = taskQueue.pop();

	if (RcppThread::isInterrupted()) {
	    break;
	}
    }
    return;
}
