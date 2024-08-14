#include "FasStableProducerConsumer.hpp"
#include "IndTestMultiCox.hpp"

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

    for (int d = 0; d <= _depth; d++) {
        bool more;

	if (verbose) Rcpp::Rcout << "    Searching at depth " << d << "...\n";

        if (d == 0) {
            more = searchAtDepth0();
        } else {

	    if (test->getData().isCensored()) {

		((IndTestMultiCox*)test)->clearWZmap();
		
		for (Node n : nodes) {
		    if (n.isCensored()) {
			// std::vector<Node> adjVars(adjacencies[n].begin(), adjacencies[n].end());
			// if (adjVars.size() == 0)
			//     continue;
			// RcppThread::Rcout << "Calling resetWZ for " + n.getName() + "..." << std::endl;
			std::vector<Node> emptySet = {};
			((IndTestMultiCox*)test)->resetWZ(n, emptySet);
		    }
		}
	    }
	    
            more = searchAtDepth(d);
        }

	// int edgeCount = 0;

	// for (const Node& node : nodes) {
	//     edgeCount += adjacencies[node].size();
	// }

	// if (verbose) Rcpp::Rcout << "        " << edgeCount/2 << " edges remaining..." << std::endl;

        if (!more) break;
    }

    graph = EdgeListGraph(nodes);

    
    // Should be more efficient than the above code
    for (const Node& x : nodes) {
	for (const Node& y : adjacencies[x]) {
	    graph.addUndirectedEdge(x, y);
	}
    }

    // if (fdr) {
    // 	int numEdges = graph.getNumEdges();
	
    // 	double harmonicSum = 0;
    // 	for (int i = 1; i <= numEdges; i++) harmonicSum += (1 / (double) i);

    // 	Rcpp::Rcout << "\nNumber of Edges = " << numEdges << std::endl;
    // 	Rcpp::Rcout << "Multiplier = " << harmonicSum << std::endl;

    // 	double alphaStar = test->getAlpha() / harmonicSum;

    // 	Rcpp::Rcout << "Final alpha* = " << alphaStar << std::endl;

    // 	test->setAlpha(alphaStar);

    // }
    

    // auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    double elapsedTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();

    if (verbose) {
        double factor = (elapsedTime < 10) ? std::pow(10, 2 - std::ceil(std::log10(std::abs(elapsedTime)))) : 1.0;
	elapsedTime = std::round(elapsedTime * factor) / factor;
        Rcpp::Rcout << "  FAS Stable Elapsed Time =  " << elapsedTime << " s" << std::endl;
    }

    // if (verbose) Rcpp::Rcout << "Fas graph: \n" << graph << std::endl;

    if (initialGraph != NULL) {
	graph.setHyperParam("lambda", initialGraph->getHyperParam("lambda"));
	// graph.setHyperParam("penalty", Rcpp::clone(initialGraph.getHyperParam("penalty")));
    }
    graph.setHyperParam("alpha", { test->getAlpha() });

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

	    if (test->getData().isCensored()) {

		((IndTestMultiCox*)test)->clearWZmap();
		
		for (Node n : nodes) {
		    if (n.isCensored()) {
			std::vector<Node> adjVars(adjacencies[n].begin(), adjacencies[n].end());
			if (adjVars.size() == 0)
			    continue;
			// RcppThread::Rcout << "\nCalling resetWZ for " + n.getName() + "..." << std::endl;
			((IndTestMultiCox*)test)->resetWZ(n, adjVars);
		    }
		}
	    }
	    
            more = searchAtDepth(d);
        }

	// int edgeCount = 0;

	// for (const Node& node : nodes) {
	//     edgeCount += adjacencies[node].size();
	// }

	// if (verbose) Rcpp::Rcout << "\t" << edgeCount/2 << " edges remaining..." << std::endl;

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
		      return edgePvals[e1] > edgePvals[e2];
		  });

	// numEdges = edgeVec.size();
	
	double harmonicSum = 1;
	if (initialGraph != NULL) {
	    harmonicSum = 0;
	    for (int i = 1; i <= numEdges; i++) harmonicSum += (1 / (double) i);
	}

	// Rcpp::Rcout << "\nNumber of  edges: " << numEdges << std::endl;
	// if (initialGraph != NULL) {
	//     Rcpp::Rcout << "Multiplier = " << harmonicSum << std::endl;
	// }

	if (numEdges * harmonicSum > nodes.size() * (nodes.size() - 1) / 2) {
	    numEdges = nodes.size() * (nodes.size() - 1) / 2;
	    harmonicSum = 1;
	}
	
	double maxFdrpval = 0;
	int firstRejectIdx = -1;
	for (int i = 0; i < edgeVec.size(); i++) {
	    NodePair edgePair = edgeVec.at(i);
	    double pval = edgePvals[edgeVec.at(i)];
	    double fdrpval = harmonicSum * numEdges / ((double) edgeVec.size()-i) * pval;
	    double alphaStar = ((double) edgeVec.size()-i) /
		(((double) numEdges) * harmonicSum) * test->getAlpha();
	    
	    // if (initialGraph != NULL) {
	    // 	fdrpval *= harmonicSum;
	    // 	alphaStar /= harmonicSum;
	    // }

	    // Rcpp::Rcout << "  Edge: " << edgePair.first << " --- "
	    // 		<< edgePair.second << "\n    p = " << pval
	    // 		<< "\n    FDR p = " << fdrpval
	    // 		<< "\n    alpha* = " << alphaStar << std::endl;

	    if (fdrpval < test->getAlpha() && firstRejectIdx < 0) {
		firstRejectIdx = i;
		break;
	    }

	    // maxFdrpval = std::max(maxFdrpval, fdrpval);

	    // if (pval <= test->getAlpha()) {
	    // 	// Rcpp::Rcout << "  Edge: " << edgePair.first << " --- "
	    // 	// 	    << edgePair.second << "\n    p = " << pval
	    // 	// 	    << "\n    FDR p = " << maxFdrpval << std::endl;
	    
	    // 	if (maxFdrpval > test->getAlpha()) {
	    // 	    sepset.set(edgePair.first, edgePair.second, edgeMaxPSet[edgePair], pval);
	    // 	    adjacencies[edgePair.first].erase(edgePair.second);
	    // 	    adjacencies[edgePair.second].erase(edgePair.first);
	    // 	}
	    // }
	}

	if (firstRejectIdx < 0)
	    firstRejectIdx = edgeVec.size();

	double minFdrpval = 1.0;
	for (int i = 0; i < firstRejectIdx; i++) {
	    NodePair edgePair = edgeVec.at(i);
	    double pval = edgePvals[edgeVec.at(i)];
	    // double fdrpval = numEdges * harmonicSum / ((double) edgeVec.size()-i) * pval;
	    // minFdrpval = std::min(minFdrpval, fdrpval);
	    sepset.set(edgePair.first, edgePair.second, edgeMaxPSet[edgePair], pval);
	    adjacencies[edgePair.first].erase(edgePair.second);
	    adjacencies[edgePair.second].erase(edgePair.first);
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

                if (!initialGraph->isAdjacentTo(x2, y2)) {
		    // std::vector<Node> nodesCopy(nodes);
		    // std::vector<Node> nodesXY = { x, y };
		    
		    // auto it = std::set_difference(nodesCopy.begin(), nodesCopy.end(),
		    // 				  nodesXY.begin(), nodesXY.end(),
		    // 				  nodesCopy.begin());

		    // nodesCopy.resize(it-nodesCopy.begin());
		    // // nodesCopy.erase(x);
		    // // nodesCopy.erase(y);
		  
		    // std::lock_guard<std::mutex> adjacencyLock(adjacencyMutex);
		    // // if (!sepset.isReturnEmptyIfNotSet()) {
		    // sepset.set(x, y, nodesCopy, 1.0);
		    // // }
                    continue;
		}
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
        bool noEdgeRequired = knowledge.noEdgeRequired(task.x, task.y);
        bool forbiddenEdge = (knowledge.isForbidden(task.x, task.y) && knowledge.isForbidden(task.y, task.x));
	
	if (fdr) {
	    std::lock_guard<std::mutex> pvalLock(pvalMutex);
	    NodePair edgePair = std::minmax(task.x, task.y);
	    if (edgePvals.find(edgePair) == edgePvals.end()) {
		edgePvals[edgePair] = pval;
		edgeMaxPSet[edgePair] = task.z;
	    } else if (edgePvals[edgePair] < pval) {
		edgePvals[edgePair] = pval;
		edgeMaxPSet[edgePair] = task.z;
	    }
	}

	if (independent && noEdgeRequired) {
	    std::lock_guard<std::mutex> adjacencyLock(adjacencyMutex);
	    
	    // if (!sepset.isReturnEmptyIfNotSet()) {
	    sepset.set(task.x, task.y, task.z, pval);
	    // }
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

	// bool edgeExists;
	// {
	//     std::lock_guard<std::mutex> adjacencyLock(adjacencyMutex);
	//     edgeExists = adjacencies[task.x].count(task.y) && adjacencies[task.y].count(task.x);
	// }
	
        // if (!edgeExists) continue; // Skip if the edge no longer exists

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
        bool noEdgeRequired = knowledge.noEdgeRequired(task.x, task.y);

	if (fdr) {
	    std::lock_guard<std::mutex> pvalLock(pvalMutex);
	    NodePair edgePair = std::minmax(task.x, task.y);
	    if (edgePvals.find(edgePair) == edgePvals.end()) {
		edgePvals[edgePair] = pval;
		edgeMaxPSet[edgePair] = task.z;
	    } else if (edgePvals[edgePair] < pval) {
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
            std::vector<Node> ppx = possibleParents(x, _adjx, y);

	    // std::sort(ppx.begin(), ppx.end());

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

std::vector<Node> FasStableProducerConsumer::possibleParents(const Node& x,
							     std::vector<Node>& adjx,
							     const Node& y) {
    std::vector<Node> possibleParents;

    for (const Node& z : adjx) {
	if (z==x) continue;
	if (z==y) continue;

	if (possibleParentOf(x, z)) {
	    possibleParents.push_back(z);
	}
    }

    return possibleParents;
}

bool FasStableProducerConsumer::possibleParentOf(const Node& x, const Node& z) {
    return !knowledge.isForbidden(z, x) && !knowledge.isRequired(x, z);
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
		      return edgePvals[e1] > edgePvals[e2];
		  });

	
	double harmonicSum = 0;
	for (int i = 1; i <= numEdges; i++) harmonicSum += (1 / (double) i);

	// Rcpp::Rcout << "\nNumber of Edges = " << numEdges << std::endl;
	// Rcpp::Rcout << "Multiplier = " << harmonicSum << std::endl;

	if (numEdges * harmonicSum > nodes.size() * (nodes.size() - 1) / 2) {
	    numEdges = nodes.size() * (nodes.size() - 1) / 2;
	    harmonicSum = 1;
	}
	
	double maxFdrpval = 0;
	int firstRejectIdx = -1;
	for (int i = 0; i < edgeVec.size(); i++) {
	    NodePair edgePair = edgeVec.at(i);
	    double pval = edgePvals[edgeVec.at(i)];
	    double fdrpval = numEdges * harmonicSum / ((double) edgeVec.size()-i) * pval;
	    double alphaStar = ((double) edgeVec.size()-i) /
		((double) numEdges * harmonicSum) * test->getAlpha();
	    
	    // Rcpp::Rcout << "  Edge: " << edgePair.first << " --- "
	    // 		<< edgePair.second << "\n    p = " << pval
	    // 		<< "\n    FDR p = " << fdrpval
	    // 		<< "\n    alpha* = " << alphaStar << std::endl;

	    if (fdrpval < test->getAlpha() && firstRejectIdx < 0) {
		firstRejectIdx = i;
		break;
	    }

	    // maxFdrpval = std::max(maxFdrpval, fdrpval);
	    
	    // if (pval <= test->getAlpha()) {
	    // 	// Rcpp::Rcout << "  Edge " << i << ": " << edgePair.first << " --- "
	    // 	// 	    << edgePair.second << "\n    p = " << pval
	    // 	// 	    << "\n    FDR p = " << maxFdrpval << std::endl;

	    // 	if (maxFdrpval > test->getAlpha()) {
	    // 	    sepset.set(edgePair.first, edgePair.second, edgeMaxPSet[edgePair], pval);
	    // 	    adjacencies[edgePair.first].erase(edgePair.second);
	    // 	    adjacencies[edgePair.second].erase(edgePair.first);
	    // 	}
	    // }
	}
	
	if (firstRejectIdx < 0)
	    firstRejectIdx = edgeVec.size();
	
	double minFdrpval = 1.0;
	for (int i = 0; i < firstRejectIdx; i++) {
	    NodePair edgePair = edgeVec.at(i);
	    double pval = edgePvals[edgeVec.at(i)];
	    // double fdrpval = numEdges * harmonicSum / ((double) edgeVec.size()-i) * pval;
	    // minFdrpval = std::min(minFdrpval, fdrpval);
	    sepset.set(edgePair.first, edgePair.second, edgeMaxPSet[edgePair], pval);
	    adjacencies[edgePair.first].erase(edgePair.second);
	    adjacencies[edgePair.second].erase(edgePair.first);
	}
    }

    return freeDegree() > depth;
}


