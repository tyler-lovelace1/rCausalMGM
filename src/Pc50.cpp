#include "Pc50.hpp"

#include "GraphUtils.hpp"
#include "BlockingQueue.hpp"

bool Pc50::isCollider(Triple& t) {
    return score[t] > 0.5;
}

bool Pc50::isNonCollider(Triple& t) {
    return score[t] <= 0.5;
}

void Pc50::orientCollider(const Node& a, const Node& b, const Node& c) {
    if (wouldCreateBadCollider(a, b)) return;
    if (wouldCreateBadCollider(c, b)) return;
    if (graph.getEdges(a, b).size() > 1) return;
    if (graph.getEdges(b, c).size() > 1) return;
    graph.removeEdge(a, b);
    graph.removeEdge(c, b);
    graph.addDirectedEdge(a, b);
    graph.addDirectedEdge(c, b);
    // Rcpp::Rcout << "ORIENTED SUCCESSFULLY" << std::endl;
}

bool Pc50::wouldCreateBadCollider(const Node& x, const Node& y) {
    std::unordered_set<Node> empty = {};
    std::unordered_set<Node> ySet = {y};

    for (const Node& z : graph.getAdjacentNodes(y)) {
        if (x == z) continue;

        // if (!graph.isAdjacentTo(x, z) &&
	//     graph.getEndpoint(z, y) == ENDPOINT_ARROW &&
	//     !sepset(x, z, empty, ySet)) {
	//     return true;
	// }
	if (!graph.isAdjacentTo(x, z) &&
	    graph.getEndpoint(z, y) == ENDPOINT_ARROW &&
	    !colliders[Triple(x, y, z)]) {
	    return true;
	}
    }

    return false;
}

void Pc50::orientUnshieldedTriples() {

    BlockingQueue<ColliderTask> taskQueue(100);

    sepsetCount.clear();
    score.clear();
    colliders.clear();

    auto producer = [&]() {
        std::vector<Node> nodes = graph.getNodes();

	int maxDepth;

        for (const Node& y : nodes) {
            std::vector<Node> adjacentNodes = graph.getAdjacentNodes(y);

            if (adjacentNodes.size() < 2)
                continue;

	    std::sort(adjacentNodes.begin(),
		      adjacentNodes.end(),
		      [] (const Node& a, const Node& b) { return a < b; }
		);
	    
            ChoiceGenerator cg(adjacentNodes.size(), 2);
            std::vector<int> *combination;
            for (combination = cg.next(); combination != NULL; combination = cg.next()) {
                Node x = adjacentNodes[(*combination)[0]];
                Node z = adjacentNodes[(*combination)[1]];

                if (graph.isAdjacentTo(x, z))
                    continue;

                std::vector<Node> adjx = graph.getAdjacentNodes(x);
                std::vector<Node> adjz = graph.getAdjacentNodes(z);

		taskQueue.push(ColliderTask(Triple(x, y, z), {}));

		maxDepth = std::max(adjx.size(), adjz.size());

		if (depth != -1) maxDepth = std::min(depth, maxDepth);

                for (int d = 1; d <= maxDepth; d++) {
                    if (adjx.size() >= 2 && d <= adjx.size()) {
                        ChoiceGenerator gen(adjx.size(), d);

                        std::vector<int> *choice;
                        for (choice = gen.next(); choice != NULL; choice = gen.next()) {
                            std::vector<Node> v = GraphUtils::asList(*choice, adjx);
                            taskQueue.push(ColliderTask(Triple(x, y, z), v));
                        }
                    }

                    if (adjz.size() >= 2 && d <= adjz.size()) {
                        ChoiceGenerator gen(adjz.size(), d);

                        std::vector<int> *choice;
                        for (choice = gen.next(); choice != NULL; choice = gen.next()) {
                            std::vector<Node> v = GraphUtils::asList(*choice, adjz);
                            taskQueue.push(ColliderTask(Triple(x, y, z), v));
                        }
                    }
                }
            }
	    
	    if (RcppThread::isInterrupted()) {
	      break;
	    }
        }

        //Poison Pill
        for (int i = 0; i < parallelism; i++) {
            taskQueue.push(ColliderTask());
        }
    };

    auto consumer = [&]() {
        while(true) {
            ColliderTask ct = taskQueue.pop();
            Triple t = ct.t;
            std::vector<Node> sepset = ct.sepset;

            // Poison pill
            if (t.x.isNull() && t.y.isNull() && t.z.isNull()) {
                return;
            }

	    bool indep = independenceTest->isIndependent(t.x, t.z, sepset);

	    {
		std::lock_guard<std::mutex> mapLock(mapMutex);
		
		if (sepsetCount.count(Triple(t.x, t.y, t.z)) == 0) {
		    sepsetCount[Triple(t.x, t.y, t.z)] = {0, 0};
		}
		
		if (indep) {
		    if (std::find(sepset.begin(), sepset.end(), t.y) != sepset.end()) {
			std::pair<int, int> current = sepsetCount[Triple(t.x, t.y, t.z)];
			sepsetCount[Triple(t.x, t.y, t.z)] = {current.first + 1, current.second};
		    } else {
			std::pair<int, int> current = sepsetCount[Triple(t.x, t.y, t.z)];
			sepsetCount[Triple(t.x, t.y, t.z)] = {current.first, current.second + 1};
		    }
		}
	    }
        }
    };

    if (parallelism <= 0) {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }

    std::vector<RcppThread::Thread> threads;

    threads.push_back(RcppThread::Thread( producer ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( consumer ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    std::vector<Triple> colliderList;
    // colliderList.reserve(sepsetCount.size());
    
    // std::ofstream logfile;
    // logfile.open("pc50_debug.log", std::ios_base::app);
    for (auto element : sepsetCount) {
        Triple t = element.first;
	std::pair<int, int> count = element.second;
	colliders[t] = false;

	score[t] = count.second / (double) std::max(count.first + count.second, 1);

	// logfile << t << ":   {" << count.first << "," << count.second
	// 	<< "}   :   " << score[t] << std::endl;

        if (isCollider(t)) {
            // colliderAllowed (knowledge)
            if (true) {
		colliders[t] = true;
		colliderList.push_back(t);
                // graph.setEndpoint(t.x, t.y, ENDPOINT_ARROW);
                // graph.setEndpoint(t.z, t.y, ENDPOINT_ARROW);
            }
        }
    }
    // logfile.close();

    // Most independent ones first.
    std::sort(colliderList.begin(), colliderList.end(),
	      [&](const Triple& t1, const Triple& t2) {
		  if (score[t1] == score[t2])
		      return t1 < t2;
		  return score[t1] > score[t2];
	      });

    for (Triple triple : colliderList) {
        Node a = triple.getX();
        Node b = triple.getY();
        Node c = triple.getZ();

        if (!(graph.getEndpoint(b, a) == ENDPOINT_ARROW
	      || graph.getEndpoint(b, c) == ENDPOINT_ARROW)) {
	    // Rcpp::Rcout << "orienting collider " << triple << "  :  " << score[triple] << std::endl;
            orientCollider(a, b, c);
        }
    }
}

/**
 * Constructs a PC50 algorithm that uses the given independence test as oracle. This does not make a copy of the
 * independence test, for fear of duplicating the data set!
 */
Pc50::Pc50(IndependenceTest *independenceTest) {
    if (independenceTest == NULL) {
        throw std::invalid_argument("independenceTest may not be NULL.");
    } 

    this->independenceTest = independenceTest;

}

/**
 * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
 * checked.
 *
 * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
 *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
 *              machines.
 */
void Pc50::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

std::unordered_set<Edge> Pc50::getAdjacencies() {
    return graph.getEdges();
}

std::unordered_set<Edge> Pc50::getNonadjacencies() {
    EdgeListGraph complete = GraphUtils::completeGraph(graph);
    std::unordered_set<Edge> nonAdjacencies = complete.getEdges();
    EdgeListGraph undirected = GraphUtils::undirectedGraph(graph);
    for (Edge edge : undirected.getEdges()) {
        nonAdjacencies.erase(edge);
    }
    return nonAdjacencies;
}

/**
 * Runs PC starting with a fully connected graph over all of the variables in the domain of the independence test.
 * See PC for caveats. The number of possible cycles and bidirected edges is far less with PC50 than with PC.
 */
EdgeListGraph Pc50::search() {
    return search(independenceTest->getVariables());
}

EdgeListGraph Pc50::search(const std::vector<Node>& nodes) {
    FasStableProducerConsumer fas(initialGraph, independenceTest, parallelism);
    return search(fas, nodes);
}

EdgeListGraph Pc50::search(FasStableProducerConsumer& fas, const std::vector<Node>& nodes) {
    if (verbose) Rcpp::Rcout << "Starting PC50 algorithm..." << std::endl;

    // allTriples = {};

    if (independenceTest == NULL)
        throw std::invalid_argument("independenceTest of Pc50 may not be NULL.");

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<Node> allNodes = independenceTest->getVariables();

    for (const Node& node : nodes) {
        if (std::find(allNodes.begin(), allNodes.end(), node) == allNodes.end())
            throw std::invalid_argument("All of the given nodes must be in the domain of the independence test provided.");
    }

    fas.setDepth(depth);
    fas.setVerbose(verbose);
    fas.setFDR(fdr);

    // Note that we are ignoring the sepset map returned by this method
    // on purpose; it is not used in this search.
    graph = fas.search();
    sepsets = fas.getSepsets();

    if (verbose) Rcpp::Rcout << "  Orienting edges..." << std::endl;

    orientUnshieldedTriples();

    // Rcpp::Rcout << graph << std::endl;

    MeekRules meekRules;
    meekRules.orientImplied(graph);

    // // Set algorithm and type
    // std::ostringstream alg;
    // alg << "Pc50: alpha = " << independenceTest->getAlpha();
    // graph.setAlgorithm(alg.str());
    // graph.setGraphType("markov equivalence class"); 

    // Set algorithm and type
    std::ostringstream alg;
    if (initialGraph==NULL) {
	alg << "PC50";
    } else {
	alg << initialGraph->getAlgorithm() << "-" << "PC50";
    }
    graph.setAlgorithm(alg.str());
    graph.setGraphType("completed partially directed acyclic graph");
    
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    if (verbose) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "PC50 Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

    return graph;
}
