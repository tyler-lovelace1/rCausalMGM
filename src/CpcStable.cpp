#include "CpcStable.hpp"

#include "GraphUtils.hpp"
#include "BlockingQueue.hpp"

bool CpcStable::isCollider(Triple t) {
    std::pair<int, int> count = sepsetCount[t];

    return (count.first == 0) && (count.second > 0);
}

bool CpcStable::isNonCollider(Triple t) {
    std::pair<int, int> count = sepsetCount[t];

    return (count.first > 0) && (count.second == 0);
}

void CpcStable::orientUnshieldedTriples() {

    BlockingQueue<ColliderTask> taskQueue(10000);

    auto producer = [&]() {
        std::vector<Variable*> nodes = graph.getNodes();

        for (Variable* y : nodes) {
            std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(y);

            if (adjacentNodes.size() < 2)
                continue;

	    std::sort(adjacentNodes.begin(),
		      adjacentNodes.end(),
		      [] (Variable* a, Variable* b) {return a->getName() < b->getName(); }
		);
	    
            ChoiceGenerator cg(adjacentNodes.size(), 2);
            std::vector<int> *combination;
            for (combination = cg.next(); combination != NULL; combination = cg.next()) {
                Variable* x = adjacentNodes[(*combination)[0]];
                Variable* z = adjacentNodes[(*combination)[1]];
		// if (adjacentNodes[(*combination)[0]]->getName()
		//     < adjacentNodes[(*combination)[1]]->getName()) {
		//     x = adjacentNodes[(*combination)[0]];
		//     z = adjacentNodes[(*combination)[1]];
		// } else {
		//     x = adjacentNodes[(*combination)[1]];
		//     z = adjacentNodes[(*combination)[0]];
		// }

                if (graph.isAdjacentTo(x, z))
                    continue;

                // std::unique_lock<std::mutex> mapLock(mapMutex);
                // sepsetCount[Triple(x, y, z)] = {0, 0};
                // mapLock.unlock();

                std::vector<Variable*> adjx = graph.getAdjacentNodes(x);
                std::vector<Variable*> adjz = graph.getAdjacentNodes(z);

		taskQueue.push(ColliderTask(Triple(x, y, z), {}));
		// taskQueue.push(ColliderTask(Triple(z, y, x), {}));

                for (int d = 1; d <= std::max(adjx.size(), adjz.size()); d++) {
                    if (adjx.size() >= 2 && d <= adjx.size()) {
                        ChoiceGenerator gen(adjx.size(), d);

                        std::vector<int> *choice;
                        for (choice = gen.next(); choice != NULL; choice = gen.next()) {
                            std::vector<Variable*> v = GraphUtils::asList(*choice, adjx);
                            taskQueue.push(ColliderTask(Triple(x, y, z), v));
			    // taskQueue.push(ColliderTask(Triple(z, y, x), v));
                        }
                    }

                    if (adjz.size() >= 2 && d <= adjz.size()) {
                        ChoiceGenerator gen(adjz.size(), d);

                        std::vector<int> *choice;
                        for (choice = gen.next(); choice != NULL; choice = gen.next()) {
                            std::vector<Variable*> v = GraphUtils::asList(*choice, adjz);
                            taskQueue.push(ColliderTask(Triple(x, y, z), v));
			    // taskQueue.push(ColliderTask(Triple(z, y, x), v));
                        }
                    }
                }
            }
        }

        //Poison Pill
        for (int i = 0; i < parallelism; i++) {
            taskQueue.push(ColliderTask(Triple(NULL, NULL, NULL), {}));
        }
    };

    auto consumer = [&]() {
	// std::ofstream logfile;
        while(true) {
            ColliderTask ct = taskQueue.pop();
            Triple t = ct.t;
            std::vector<Variable*> sepset = ct.sepset;

            // Poison pill
            if (t.x == NULL && t.y == NULL && t.z == NULL) {
                return;
            }

	    bool indep = independenceTest->isIndependent(t.x, t.z, sepset);

	    {
		std::unique_lock<std::mutex> mapLock(mapMutex);
		mapCondition.wait(mapLock,
			      [this] {
				  if (!mapModifying) {
				      return mapModifying = true;
				  }
				  return false;
			      });
		if (sepsetCount.count(Triple(t.x, t.y, t.z)) == 0) {
		    sepsetCount[Triple(t.x, t.y, t.z)] = {0, 0};
		}
		// std::sort(sepset.begin(),
		// 	  sepset.end(),
		// 	  [] (Variable* a, Variable* b) {return a->getName() < b->getName(); }
		//     );
		// logfile.open("../cpc_debug.log", std::ios_base::app);
		// logfile << mapLock.owns_lock() << "\t" << t.x->getName() << " ? " << t.z->getName() << " | [";
		// for (Variable* n : sepset) logfile << n->getName() << ",";
		// logfile << "]\t" << indep << "\n";
		// logfile.close();
		if (indep) {
		    // if (sepset.contains(t.y))
		    if (std::find(sepset.begin(), sepset.end(), t.y) != sepset.end()) {
			std::pair<int, int> current = sepsetCount[Triple(t.x, t.y, t.z)];
			sepsetCount[Triple(t.x, t.y, t.z)] = {current.first + 1, current.second};
		    } else {
			// std::lock_guard<std::mutex> mapLock(mapMutex);
			std::pair<int, int> current = sepsetCount[Triple(t.x, t.y, t.z)];
			sepsetCount[Triple(t.x, t.y, t.z)] = {current.first, current.second + 1};
		    }
		}
		mapModifying = false;
	    }
	    mapCondition.notify_one();
        }
    };

    std::vector<std::thread> threads;

    threads.push_back(std::thread( producer ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( consumer ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    // std::ofstream logfile;
    // logfile.open("../cpc_debug.log", std::ios_base::app);
    for (auto element : sepsetCount) {
        Triple t = element.first;

	// logfile << t << ":   {" << element.second.first << "," << element.second.second << "}\n";
        
        if (isCollider(t)) {
            // colliderAllowed (knowledge)
            if (true) {
                graph.setEndpoint(t.x, t.y, ENDPOINT_ARROW);
                graph.setEndpoint(t.z, t.y, ENDPOINT_ARROW);
            }
        } else if (!isNonCollider(t)) {
            graph.addAmbiguousTriple(t.x, t.y, t.z);
        }

        allTriples.insert(t);
    }
    // logfile.close();
}

/**
 * Constructs a CPC algorithm that uses the given independence test as oracle. This does not make a copy of the
 * independence test, for fear of duplicating the data set!
 */
CpcStable::CpcStable(IndependenceTest *independenceTest) {
    if (independenceTest == NULL) {
        throw std::invalid_argument("independenceTest may not be NULL.");
    } 

    this->independenceTest = independenceTest;

    if (parallelism == 0) {
        parallelism = 4;
        Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
    }
}

/**
 * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
 * checked.
 *
 * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
 *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
 *              machines.
 */
void CpcStable::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

std::unordered_set<Edge> CpcStable::getAdjacencies() {
    return graph.getEdges();
}

std::unordered_set<Edge> CpcStable::getNonadjacencies() {
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
 * See PC for caveats. The number of possible cycles and bidirected edges is far less with CPC than with PC.
 */
EdgeListGraph CpcStable::search() {
    return search(independenceTest->getVariables());
}

EdgeListGraph CpcStable::search(const std::vector<Variable*>& nodes) {
    FasStableProducerConsumer fas(initialGraph, independenceTest);
    return search(fas, nodes);
}

EdgeListGraph CpcStable::search(FasStableProducerConsumer& fas, const std::vector<Variable*>& nodes) {
    if (verbose) Rcpp::Rcout << "Starting CPC algorithm" << std::endl;

    allTriples = {};

    if (independenceTest == NULL)
        throw std::invalid_argument("independenceTest of CpcStable may not be NULL.");

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<Variable*> allNodes = independenceTest->getVariables();

    for (Variable* node : nodes) {
        if (std::find(allNodes.begin(), allNodes.end(), node) == allNodes.end())
            throw std::invalid_argument("All of the given nodes must be in the domain of the independence test provided.");
    }

    fas.setDepth(depth);
    fas.setVerbose(verbose);

    // Note that we are ignoring the sepset map returned by this method
    // on purpose; it is not used in this search.
    graph = fas.search();
    sepsets = fas.getSepsets();

    orientUnshieldedTriples();

    MeekRules meekRules;
    meekRules.orientImplied(graph);

    // Set algorithm and type
    std::ostringstream alg;
    alg << "CpcStable: alpha = " << independenceTest->getAlpha();
    graph.setAlgorithm(alg.str());
    graph.setGraphType("markov equivalence class");

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    if (verbose) {
        Rcpp::Rcout.precision(2);
        Rcpp::Rcout << "CpcStable Elapsed time =  " << (elapsedTime / 1000.0) << " s" << std::endl;
    }

    return graph;
}
