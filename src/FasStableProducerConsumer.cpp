#include "FasStableProducerConsumer.hpp"
#include "IndTestMultiScaledZ.hpp"

#include "ChoiceGenerator.hpp"
#include "GraphUtils.hpp"

// Node FasStableProducerConsumer::nullNode = Node();

FasStableProducerConsumer::FasStableProducerConsumer(EdgeListGraph *initialGraph, IndependenceTest *test, int threads) : FasStableProducerConsumer(test, threads) 
{
    this->initialGraph = initialGraph;
}

FasStableProducerConsumer::FasStableProducerConsumer(IndependenceTest *test, int threads) : taskQueue(MAX_QUEUE_SIZE) 
{
    this->test = test;
    this->nodes = test->getVariables();

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }
}

void FasStableProducerConsumer::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    this->depth = depth;
}

/**
 * Discovers all adjacencies in data.  The procedure is to remove edges in the graph which connect pairs of
 * variables which are independent conditional on some other set of variables in the graph (the "sepset"). These are
 * removed in tiers.  First, edges which are independent conditional on zero other variables are removed, then edges
 * which are independent conditional on one other variable are removed, then two, then three, and so on, until no
 * more edges can be removed from the graph.  The edges which remain in the graph after this procedure are the
 * adjacencies in the data.
 *
 * @return a SepSet, which indicates which variables are independent conditional on which other variables
 */
EdgeListGraph FasStableProducerConsumer::search() {
    auto startTime = std::chrono::high_resolution_clock::now();
    
    if (verbose) Rcpp::Rcout << "  Starting FAS Stable..." << std::endl;

    sepset = SepsetMap();
    sepset.setReturnEmptyIfNotSet(sepsetsReturnEmptyIfNotFixed);

    int _depth = depth;

    if (_depth == -1) _depth = 1000;

    adjacencies.clear();

    for (const Node& node : nodes) {
        adjacencies[node] = {};
    }

    int maxDepth;

    for (int d = 0; d <= _depth; d++) {
        bool more;

	if (verbose) Rcpp::Rcout << "    Searching at depth " << d << "...\r";

        if (d == 0) {
            more = searchAtDepth0();
        } else {

	    // for (Node n : nodes) {
	    // 	if (n.isCensored()) {
	    // 	    std::vector<Node> adjVars(adjacencies[n].begin(), adjacencies[n].end());
	    // 	    if (adjVars.size() == 0)
	    // 		continue;
	    // 	    // RcppThread::Rcout << "Calling scaledZ..." << std::endl;
	    // 	    ((IndTestMultiScaledZ*)test)->scaledZ(n, adjVars);
	    // 	}
	    // }
	    
            more = searchAtDepth(d);
        }

	maxDepth = d;
	
        if (!more) break;
    }

    graph = EdgeListGraph(nodes);

    // int numEdges = nodes.size() * (nodes.size() - 1) / 2;
    // if (initialGraph != NULL) {
    // 	numEdges = initialGraph->getNumEdges();
    // }

    // if (true) {
    // 	std::vector<NodePair> edgeVec;
    // 	std::unordered_set<NodePair,
    // 			   boost::hash<NodePair>> edgeSet;

    // 	for (int i = 0; i < nodes.size(); i++) {
    // 	    const Node& x = nodes[i];
	    
    // 	    for (int j = i+1; j < nodes.size(); j++) {
    // 		const Node& y = nodes[j];

    // 		if (initialGraph != NULL) {
    // 		    const Node& x2 = initialGraph->getNode(x.getName());
    // 		    const Node& y2 = initialGraph->getNode(y.getName());

    // 		    if (!initialGraph->isAdjacentTo(x2, y2))
    // 			continue;
    // 		}

    // 		edgeSet.insert(std::minmax(x,y));
    // 	    }
    // 	}

    // 	edgeVec = std::vector<NodePair>(edgeSet.begin(), edgeSet.end());

	
    // 	std::sort(edgeVec.begin(), edgeVec.end(),
    // 		  [&](const NodePair& e1,
    // 		      const NodePair& e2) {
    // 		      return edgePvals[e1] < edgePvals[e2];
    // 		  });

    // 	// Rcpp::Function qvalue("qvalue");

    // 	// arma::rowvec pvals(numEdges, arma::fill::ones);

    // 	// for (arma::uword i = 0; i < numEdges; i++) {
    // 	//     pvals(i) = edgePvals[edgeVec[i]];
    // 	// }

    // 	// pvals.transform( [](double val) { return (val < 1e-16) ? 1e-16 : val; } );
	
    // 	// Rcpp::Rcout << "\np-vals:\n" << pvals << std::endl;

    // 	// arma::vec lambda = arma::linspace(0.01, 0.99, 99);

    // 	// // lambda = nodes.size() * lambda / (2 *  numEdges);
    // 	// // lambda = 1 - lambda;

    // 	// // Rcpp::Rcout << "lambda: " << lambda.t() << std::endl;

    // 	// Rcpp::List qobj = qvalue(pvals,
    // 	// 			 Rcpp::_["pfdr"]=true,
    // 	// 			 Rcpp::_["lambda"]=lambda,
    // 	// 			 Rcpp::_["pi0.method"]="smoother",
    // 	// 			 Rcpp::_["smooth.log.pi0"]=true);

    // 	// arma::rowvec qvals = Rcpp::as<arma::rowvec>(qobj["qvalues"]);

    // 	// double pi0 = Rcpp::as<double>(qobj["pi0"]);

    // 	// Rcpp::Rcout << "pi0:\n" << pi0 << std::endl;

    // 	// Rcpp::Rcout << "q-vals:\n" << qvals << std::endl;

    // 	// numEdges = edgeVec.size();
	
    // 	// double harmonicSum = 0;
    // 	// if (initialGraph != NULL) {
    // 	//     for (int i = 1; i <= numEdges; i++) harmonicSum += (1 / (double) i);
    // 	// }

    // 	// // Rcpp::Rcout << "Number of edges: " << numEdges << std::endl;
    // 	// // if (initialGraph != NULL) {
    // 	// //     Rcpp::Rcout << "Multiplier = " << harmonicSum << std::endl;
    // 	// // }
	
    // 	// double maxFdrpval = 0;
    // 	// for (int i = 0; i < edgeVec.size(); i++) {
    // 	//     NodePair edgePair = edgeVec.at(i);
    // 	//     double pval = edgePvals[edgeVec.at(i)];
    // 	//     double fdrpval = numEdges / ((double) i+1) * pval;
	    
    // 	//     if (initialGraph != NULL) {
    // 	// 	fdrpval *= harmonicSum;
    // 	//     }

    // 	//     maxFdrpval = std::max(maxFdrpval, fdrpval);


    // 	//     if (pval <= test->getAlpha()) {
    // 	// 	Rcpp::Rcout << "  Edge " << i << ": " << edgePair.first << " --- "
    // 	// 		    << edgePair.second << "\n    p = " << pval
    // 	// 		    << "\n    FDR p = " << maxFdrpval << std::endl;
	    
    // 	// 	if (maxFdrpval > test->getAlpha()) {
    // 	// 	    sepset.set(edgePair.first, edgePair.second, edgeMaxPSet[edgePair], pval);
    // 	// 	    adjacencies[edgePair.first].erase(edgePair.second);
    // 	// 	    adjacencies[edgePair.second].erase(edgePair.first);
    // 	// 	}
    // 	//     }
    // 	// }
    // }
    
    // Should be more efficient than the above code
    for (const Node& x : nodes) {
	for (const Node& y : adjacencies[x]) {
	    graph.addUndirectedEdge(x, y);
	}
    }
    

    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    if (verbose) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << std::endl << "  FAS Stable Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

    // if (verbose) Rcpp::Rcout << "Fas graph: \n" << graph << std::endl;

    if (initialGraph != NULL) {
	graph.setHyperParam("lambda", Rcpp::clone(initialGraph->getHyperParam("lambda")));
	// graph.setHyperParam("penalty", Rcpp::clone(initialGraph.getHyperParam("penalty")));
    }
    graph.setHyperParam("alpha", Rcpp::NumericVector::create(test->getAlpha()));

    return graph;
}

std::unordered_map<Node, std::unordered_set<Node>> FasStableProducerConsumer::searchMapOnly() {
    if (verbose) Rcpp::Rcout << "Starting FasStableProducerConsumer Adjacency Search." << std::endl;

    graph.removeEdges(graph.getEdgeList());

    sepset = SepsetMap();

    int _depth = depth;

    if (_depth == -1) _depth = 1000;

    std::unordered_map<Node, std::unordered_set<Node>> adjacencies;
    std::vector<Node> nodes = graph.getNodes();

    for (const Node& node : nodes) {
        adjacencies[node] = {};
    }

    for (int d = 0; d <= _depth; d++) {
        bool more;

	if (verbose) Rcpp::Rcout << "Searching at depth " << d << "...\r";

        if (d == 0) {
            more = searchAtDepth0();
        } else {
	    // for (Node n : nodes) {
	    // 	if (n.isCensored()) {
	    // 	    std::vector<Node> adjVars(adjacencies[n].begin(), adjacencies[n].end());
	    // 	    if (adjVars.size() == 0)
	    // 		continue;
	    // 	    // RcppThread::Rcout << "Calling scaledZ..." << std::endl;
	    // 	    ((IndTestMultiScaledZ*)test)->scaledZ(n, adjVars);
	    // 	}
	    // }
            more = searchAtDepth(d);
        }

	int edgeCount = 0;

	for (const Node& node : nodes) {
	    edgeCount += adjacencies[node].size();
	}

	if (verbose) Rcpp::Rcout << "\t" << edgeCount/2 << " edges remaining..." << std::endl;

        if (!more) break;
    }

    if (verbose) Rcpp::Rcout << std::endl << "Finishing FasStableProducerConsumer Adjacency Search." << std::endl;
    return adjacencies;
}

bool FasStableProducerConsumer::searchAtDepth0() {

    std::vector<RcppThread::Thread> threads;

    int numEdges = nodes.size() * (nodes.size() - 1) / 2;
    if (initialGraph != NULL) {
    	numEdges = initialGraph->getNumEdges();
    }

    threads.push_back(RcppThread::Thread( [this] { producerDepth0(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( [this] { consumerDepth0(); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        if (threads[i].joinable()) {
	    threads[i].join();
	} else {
	    Rcpp::Rcout << "#### THREAD " << i << " NOT JOINABLE ####\n";
	}
    }

    if (fdr) {
	std::vector<NodePair> edgeVec;
	std::unordered_set<NodePair,
			   boost::hash<NodePair>> edgeSet;

	for (int i = 0; i < nodes.size(); i++) {
	    const Node& x = nodes[i];
	    
	    for (int j = i+1; j < nodes.size(); j++) {
		const Node& y = nodes[j];

		if (initialGraph != NULL) {
		    const Node& x2 = initialGraph->getNode(x.getName());
		    const Node& y2 = initialGraph->getNode(y.getName());

		    if (!initialGraph->isAdjacentTo(x2, y2))
			continue;
		}

		edgeSet.insert(std::minmax(x,y));
	    }
	}

	edgeVec = std::vector<NodePair>(edgeSet.begin(), edgeSet.end());

	
	std::sort(edgeVec.begin(), edgeVec.end(),
		  [&](const NodePair& e1,
		      const NodePair& e2) {
		      return edgePvals[e1] < edgePvals[e2];
		  });

	// numEdges = edgeVec.size();
	
	double harmonicSum = 0;
	if (initialGraph != NULL) {
	    for (int i = 1; i <= numEdges; i++) harmonicSum += (1 / (double) i);
	}

	// Rcpp::Rcout << "Number of  edges: " << numEdges << std::endl;
	// if (initialGraph != NULL) {
	//     Rcpp::Rcout << "Multiplier = " << harmonicSum << std::endl;
	// }

	// harmonicSum = 1;
	
	double maxFdrpval = 0;
	for (int i = 0; i < edgeVec.size(); i++) {
	    NodePair edgePair = edgeVec.at(i);
	    double pval = edgePvals[edgeVec.at(i)];
	    double fdrpval = numEdges / ((double) i+1) * pval;
	    
	    if (initialGraph != NULL) {
		fdrpval *= harmonicSum;
	    }

	    maxFdrpval = std::max(maxFdrpval, fdrpval);


	    if (pval <= test->getAlpha()) {
		// Rcpp::Rcout << "  Edge: " << edgePair.first << " --- "
		// 	    << edgePair.second << "\n    p = " << pval
		// 	    << "\n    FDR p = " << maxFdrpval << std::endl;
	    
		if (maxFdrpval > test->getAlpha()) {
		    sepset.set(edgePair.first, edgePair.second, edgeMaxPSet[edgePair], pval);
		    adjacencies[edgePair.first].erase(edgePair.second);
		    adjacencies[edgePair.second].erase(edgePair.first);
		}
	    }
	}
    }

    return freeDegree() > 0;
}

void FasStableProducerConsumer::producerDepth0() {
    std::vector<Node> empty = {};

    for (int i = 0; i < nodes.size(); i++) {
        const Node& x = nodes[i];

        for (int j = i+1; j < nodes.size(); j++) {
            const Node& y = nodes[j];

            if (initialGraph != NULL) {
                const Node& x2 = initialGraph->getNode(x.getName());
                const Node& y2 = initialGraph->getNode(y.getName());

                if (!initialGraph->isAdjacentTo(x2, y2))
                    continue;
            }

            taskQueue.push(IndependenceTask(x, y, empty));
        }

	// RcppThread::checkUserInterrupt();
	if (RcppThread::isInterrupted()) {
	  break;
	}
    }

    // RcppThread::Rcout << "\nCreating poison pills\n";

    // add poison pill to stop consumers
    IndependenceTask poisonPill; // (Node(), Node(), empty);

    // RcppThread::Rcout << "Poison pills created\n";
    
    // RcppThread::Rcout << "X is NULL:  " << ((poisonPill.x.isNull()) ? "true" : "false")
    // 		      << "\nY is NULL:  " << ((poisonPill.y.isNull()) ? "true" : "false")
    // 		      << std::endl;

    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }

}

void FasStableProducerConsumer::consumerDepth0() {
    while(true) {
        IndependenceTask task = taskQueue.pop();

	// RcppThread::Rcout << "X = " << task.x << "\nY = " << task.y << "\n";
        //If poison, return
        if (task.x.isNull() && task.y.isNull())  return;

        numIndependenceTests++;
	double pval = 0.0;
        bool independent = test->isIndependent(task.x, task.y, task.z, &pval);

	if (independent) {
            numIndependenceJudgements++;
        } else {
            numDependenceJudgement++;
        }

        // Knowledge
        bool noEdgeRequired = true;
        bool forbiddenEdge = false;
	
	if (true) {
	    std::lock_guard<std::mutex> pvalLock(pvalMutex);
	    NodePair edgePair = std::minmax(task.x, task.y);
	    if (edgePvals.find(edgePair) == edgePvals.end() ||
		(edgePvals[edgePair] < pval && edgePvals[edgePair] < test->getAlpha())) {
		edgePvals[edgePair] = pval;
		edgeMaxPSet[edgePair] = task.z;
	    }
	}

	if (independent && noEdgeRequired) {
	    std::lock_guard<std::mutex> adjacencyLock(adjacencyMutex);
	    
	    if (!sepset.isReturnEmptyIfNotSet()) {
		sepset.set(task.x, task.y, task.z, pval);
	    }
	} else if (!forbiddenEdge) {
	    std::lock_guard<std::mutex> adjacencyLock(adjacencyMutex);
	    
	    adjacencies[task.x].insert(task.y);
	    adjacencies[task.y].insert(task.x);
	}
    }
}

void FasStableProducerConsumer::consumerDepth(int depth) {
    while(true) {
        IndependenceTask task = taskQueue.pop();

        //If poison, return
        if (task.x.isNull() && task.y.isNull()) return;

	bool edgeExists;
	{
	    std::lock_guard<std::mutex> adjacencyLock(adjacencyMutex);
	    edgeExists = adjacencies[task.x].count(task.y) && adjacencies[task.y].count(task.x);
	}
	
        if (!edgeExists) continue; // Skip if the edge no longer exists

        numIndependenceTests++;
	double pval = 0.0;
        bool independent;
        independent = test->isIndependent(task.x, task.y, task.z, &pval);

	if (independent) {
            numIndependenceJudgements++;
        } else {
            numDependenceJudgement++;
        }

        // Knowledge
        bool noEdgeRequired = true;

	if (true) {
	    std::lock_guard<std::mutex> pvalLock(pvalMutex);
	    NodePair edgePair = std::minmax(task.x, task.y);
	    if (edgePvals.find(edgePair) == edgePvals.end() ||
		(edgePvals[edgePair] < pval && edgePvals[edgePair] < test->getAlpha())) {
		edgePvals[edgePair] = pval;
		edgeMaxPSet[edgePair] = task.z;
	    }
	}

	if (independent && noEdgeRequired) {
	    std::lock_guard<std::mutex> adjacencyLock(adjacencyMutex);
	    
	    adjacencies[task.x].erase(task.y);
	    adjacencies[task.y].erase(task.x);
	    sepset.set(task.x, task.y, task.z, pval);

	}
    }
}

void FasStableProducerConsumer::producerDepth(int depth, std::unordered_map<Node, std::unordered_set<Node>>& adjacenciesCopy) {
    for (const Node& x : nodes) {

        std::unordered_set<Node> adjx = adjacenciesCopy[x];

        for (const Node& y : adjx) {

            std::vector<Node> _adjx(adjx.begin(), adjx.end());
            _adjx.erase(std::remove(_adjx.begin(), _adjx.end(), y), _adjx.end());

            // Knowledge: possible parents
            std::vector<Node> ppx = _adjx;

	    std::sort(ppx.begin(), ppx.end());

            if (ppx.size() >= depth) {
                ChoiceGenerator cg(ppx.size(), depth);
                std::vector<int> *choice;

                for (choice = cg.next(); choice != NULL; choice = cg.next()) {
                    std::vector<Node> condSet = GraphUtils::asList(*choice, ppx);

		    // std::sort(condSet.begin(), condSet.end(),
		    // 	      []( const Node& lhs, const Node& rhs )
		    // 	      {
		    // 		  return lhs->getName() < rhs->getName();
		    // 	      }
		    // 	);
                    
                    taskQueue.push(IndependenceTask(x, y, condSet));
                }
            }
	    
	    // RcppThread::checkUserInterrupt();
	    if (RcppThread::isInterrupted()) {
	      break;
	    }
        }
	
	if (RcppThread::isInterrupted()) {
	  break;
	}
    }

    // add poison pill to stop consumers
    // RcppThread::Rcout << "Creating poison pills\n";
    // std::vector<Node> empty = {};
    IndependenceTask poisonPill; // (Node(), Node(), empty);
    
    // RcppThread::Rcout << "X is NULL:  " << ((poisonPill.x.isNull()) ? "true" : "false")
    // 		      << "\nY is NULL:  " << ((poisonPill.y.isNull()) ? "true" : "false")
    // 		      << std::endl;

    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }
}


int FasStableProducerConsumer::freeDegree() {
    int max = 0;

    for (const Node& x : nodes) {
        std::unordered_set<Node> opposites = adjacencies[x];

        for (const Node& y : opposites) {
            std::unordered_set<Node> adjx(opposites);
            adjx.erase(y);

            if (adjx.size() > max) {
                max = adjx.size();
            }
        }
    }

    return max;
}

bool FasStableProducerConsumer::searchAtDepth(int depth) {

    std::unordered_map<Node, std::unordered_set<Node>> adjacenciesCopy = adjacencies;

    int numEdges = 0;
    
    if (fdr) {
	for (auto it = adjacencies.begin(); it != adjacencies.end(); it++)
	    numEdges += it->second.size();
	numEdges /= 2;
    }

    std::vector<RcppThread::Thread> threads;

    threads.push_back(RcppThread::Thread( [&] { producerDepth(depth, adjacenciesCopy); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( [&] { consumerDepth(depth); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
	if (threads[i].joinable()) {
	    threads[i].join();
	} else {
	    Rcpp::Rcout << "#### THREAD " << i << " NOT JOINABLE ####\n";
	}
    }

    if (fdr) {
        std::vector<NodePair> edgeVec;
	std::unordered_set<NodePair,
			   boost::hash<NodePair>> edgeSet;

        for (const Node& x : nodes) {

	    std::unordered_set<Node> adjx = adjacenciesCopy[x];
	    
	    for (const Node& y : adjx) {
		edgeSet.insert(std::minmax(x,y));
	    }
	}

	edgeVec = std::vector<NodePair>(edgeSet.begin(), edgeSet.end());
	
	std::sort(edgeVec.begin(), edgeVec.end(),
		  [&](const NodePair& e1,
		      const NodePair& e2) {
		      return edgePvals[e1] < edgePvals[e2];
		  });

	
	double harmonicSum = 0;
	for (int i = 1; i <= numEdges; i++) harmonicSum += (1 / (double) i);

	// Rcpp::Rcout << "\nNumber of Edges = " << numEdges << std::endl;
	// Rcpp::Rcout << "Multiplier = " << harmonicSum << std::endl;

	if (numEdges * harmonicSum > nodes.size() * (nodes.size() - 1) / 2) {
	    numEdges = nodes.size() * (nodes.size() - 1) / 2;
	    harmonicSum = 1;
	}

	// harmonicSum = 1;
	
	double maxFdrpval = 0;
	for (int i = 0; i < edgeVec.size(); i++) {
	    NodePair edgePair = edgeVec.at(i);
	    double pval = edgePvals[edgeVec.at(i)];
	    double fdrpval = numEdges * harmonicSum / ((double) i+1) * pval;

	    maxFdrpval = std::max(maxFdrpval, fdrpval);
	    
	    if (pval <= test->getAlpha()) {
		// Rcpp::Rcout << "  Edge " << i << ": " << edgePair.first << " --- "
		// 	    << edgePair.second << "\n    p = " << pval
		// 	    << "\n    FDR p = " << maxFdrpval << std::endl;

		if (maxFdrpval > test->getAlpha()) {
		    sepset.set(edgePair.first, edgePair.second, edgeMaxPSet[edgePair], pval);
		    adjacencies[edgePair.first].erase(edgePair.second);
		    adjacencies[edgePair.second].erase(edgePair.first);
		}
	    }
	}
    }

    return freeDegree() > depth;
}
