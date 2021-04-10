#include "FasStableProducerConsumer.hpp"

#include "ChoiceGenerator.hpp"
#include "GraphUtils.hpp"

FasStableProducerConsumer::FasStableProducerConsumer(EdgeListGraph *initialGraph, IndependenceTest *test) : FasStableProducerConsumer(test) 
{
    this->initialGraph = initialGraph;
    this->test = test;
    this->nodes = test->getVariables();

    // Rcpp::Rcout << "FAS paralellism = " << parallelism << std::endl;
    if (parallelism == 0) {
        parallelism = 4;
        Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
    }
}

FasStableProducerConsumer::FasStableProducerConsumer(IndependenceTest *test) : taskQueue(MAX_QUEUE_SIZE) 
{
    this->test = test;
    this->nodes = test->getVariables();

    // Rcpp::Rcout << "FAS paralellism = " << parallelism << std::endl;
    if (parallelism == 0) {
        parallelism = 4;
        Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
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
    if (verbose) Rcpp::Rcout << "Starting FasStableProducerConsumer Adjacency Search." << std::endl;

    sepset = SepsetMap();
    sepset.setReturnEmptyIfNotSet(sepsetsReturnEmptyIfNotFixed);

    int _depth = depth;

    if (_depth == -1) _depth = 1000;

    adjacencies.clear();

    for (Variable* node : nodes) {
        adjacencies[node] = {};
    }

    for (int d = 0; d <= _depth; d++) {
        bool more;

	if (verbose) Rcpp::Rcout << "Searching at depth " << d << "..." << std::endl;

        if (d == 0) {
            more = searchAtDepth0();
        } else {
            more = searchAtDepth(d);
        }

	int edgeCount = 0;

	for (Variable* node : nodes) {
	    edgeCount += adjacencies[node].size();
	}

	if (verbose) Rcpp::Rcout << "\t" << edgeCount/2 << " edges remaining..." << std::endl;

        if (!more) break;
    }

    graph = EdgeListGraph(nodes);

    // for (int i = 0; i < nodes.size(); i++) {
    //     for (int j = i+1; j < nodes.size(); j++) {
    //         Variable* x = nodes[i];
    //         Variable* y = nodes[j];

    //         if (adjacencies[x].count(y)) {
    //             graph.addUndirectedEdge(x, y);
    //         }
    //     }
    // }

    // Should be more efficient than the above code
    for (Variable* x : nodes) {
        for (Variable* y : adjacencies[x]) {
            graph.addUndirectedEdge(x, y);
        }
    }

    if (verbose) Rcpp::Rcout << "Finishing FasStableProducerConsumer Adjacency Search." << std::endl;

    // if (verbose) Rcpp::Rcout << "Fas graph: \n" << graph << std::endl;

    return graph;
}

std::unordered_map<Variable*, std::unordered_set<Variable*>> FasStableProducerConsumer::searchMapOnly() {
    if (verbose) Rcpp::Rcout << "Starting FasStableProducerConsumer Adjacency Search." << std::endl;

    graph.removeEdges(graph.getEdgeList());

    sepset = SepsetMap();

    int _depth = depth;

    if (_depth == -1) _depth = 1000;

    std::unordered_map<Variable*, std::unordered_set<Variable*>> adjacencies;
    std::vector<Variable*> nodes = graph.getNodes();

    for (Variable* node : nodes) {
        adjacencies[node] = {};
    }

    for (int d = 0; d <= _depth; d++) {
        bool more;

	if (verbose) Rcpp::Rcout << "Searching at depth " << d << "..." << std::endl;

        if (d == 0) {
            more = searchAtDepth0();
        } else {
            more = searchAtDepth(d);
        }

	int edgeCount = 0;

	for (Variable* node : nodes) {
	    edgeCount += adjacencies[node].size();
	}

	if (verbose) Rcpp::Rcout << "\t" << edgeCount/2 << " edges remaining..." << std::endl;

        if (!more) break;
    }

    if (verbose) Rcpp::Rcout << "Finishing FasStableProducerConsumer Adjacency Search." << std::endl;
    return adjacencies;
}

bool FasStableProducerConsumer::searchAtDepth0() {

    std::vector<std::thread> threads;

    threads.push_back(std::thread( [this] { producerDepth0(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( [this] { consumerDepth0(); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    return freeDegree() > 0;
}

void FasStableProducerConsumer::producerDepth0() {
    std::vector<Variable*> empty = {};

    for (int i = 0; i < nodes.size(); i++) {
        Variable* x = nodes[i];

        for (int j = i+1; j < nodes.size(); j++) {
            Variable* y = nodes[j];

            if (initialGraph != NULL) {
                Variable* x2 = initialGraph->getNode(x->getName());
                Variable* y2 = initialGraph->getNode(y->getName());

                if (!initialGraph->isAdjacentTo(x2, y2))
                    continue;
            }

            taskQueue.push(IndependenceTask(x, y, empty));
        }
    }

    // add poison pill to stop consumers
    IndependenceTask poisonPill(NULL, NULL, empty);
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }

}

void FasStableProducerConsumer::consumerDepth0() {
    std::ofstream logfile;

    while(true) {
        IndependenceTask task = taskQueue.pop();

        //If poison, return
        if (task.x == NULL && task.y == NULL) return;

        numIndependenceTests++;
        bool independent = test->isIndependent(task.x, task.y, task.z);

        if (independent) {
            numIndependenceJudgements++;
        } else {
            numDependenceJudgement++;
        }

        // Knowledge
        bool noEdgeRequired = true;
        bool forbiddenEdge = false;

	{
	    std::unique_lock<std::mutex> adjacencyLock(adjacencyMutex);
	    adjacencyCondition.wait(adjacencyLock, [this] { return !adjacencyModifying; });
	    adjacencyModifying = true;
	    logfile.open("debug2.log", std::ios_base::app);
	    logfile << adjacencyLock.owns_lock() << "\t" << task.x->getName() << " ? "
		    << task.y->getName() << " | []\t" << independent << "\n";
	    if (!adjacencyLock.owns_lock()) {
		Rcpp::Rcout << "lock not owned\n";
		adjacencyLock.lock();
	    }
//	    if (adjacencyModifying) Rcpp::Rcout << "adjacency being modified elsewhere\n";
	    // adjacencyModifying = true;
	    if (independent && noEdgeRequired) {
		if (!sepset.isReturnEmptyIfNotSet()) {
		    sepset.set(task.x, task.y, task.z);
		}
	    } else if (!forbiddenEdge) {
		adjacencies[task.x].insert(task.y);
		adjacencies[task.y].insert(task.x);
	    }
	    logfile.close();
	    adjacencyModifying = false;
	    // adjacencyCondition.notify_one();
	    // adjacencyLock.unlock();
	    // adjacencyCondition.notify_one();
	}
	adjacencyCondition.notify_one();
    }
    // adjacencyModifying = false;
}

void FasStableProducerConsumer::consumerDepth(int depth) {
    std::ofstream logfile;

    while(true) {
        IndependenceTask task = taskQueue.pop();

        //If poison, return
        if (task.x == NULL && task.y == NULL) return;

	// bool edgeExists;
	// {
	//     std::unique_lock<std::mutex> adjacencyLock(adjacencyMutex);
	//     adjacencyCondition.wait(adjacencyLock, [this] { return !adjacencyModifying; });
	//     adjacencyModifying = true;
	//     edgeExists = adjacencies[task.x].count(task.y) && adjacencies[task.y].count(task.x);
	//     adjacencyModifying = false;
	// }
	// adjacencyCondition.notify_one();
        // adjacencyLock.unlock();

        // if (!edgeExists) continue; // Skip if the edge no longer exists

        numIndependenceTests++;
        bool independent;
        independent = test->isIndependent(task.x, task.y, task.z);

        if (independent) {
            numIndependenceJudgements++;
        } else {
            numDependenceJudgement++;
        }

        // Knowledge
        bool noEdgeRequired = true;

	{
	     std::unique_lock<std::mutex> adjacencyLock(adjacencyMutex);
	     adjacencyCondition.wait(adjacencyLock, [this] { return !adjacencyModifying; });
	     adjacencyModifying = true;
	     if (!adjacencyLock.owns_lock()) {
		 Rcpp::Rcout << "lock not owned\n";
		 adjacencyLock.lock();
	     }
	     // if (adjacencyModifying) Rcpp::Rcout << "adjacency being modified elsewhere\n";
	     // adjacencyModifying = true;
	     logfile.open("debug2.log", std::ios_base::app);
	     logfile << adjacencyLock.owns_lock() << "\t" << task.x->getName() << " ? " << task.y->getName() << " | [";
	     for (Variable* n : task.z) logfile << n->getName() << ",";
	     logfile << "]\t" << independent << "\n";
	     // adjacencyLock.lock();
	     if (independent && noEdgeRequired) {
		 // if (!adjacencyLock.owns_lock()) adjacencyLock.lock();
		 // if (adjacencyModifying) Rcpp::Rcout << "adjacency being modified elsewhere\n";

		 adjacencies[task.x].erase(task.y);
		 adjacencies[task.y].erase(task.x);
		 sepset.set(task.x, task.y, task.z);

		 // adjacencyLock.unlock();
		 
	     }
	     logfile.close();
	     adjacencyModifying = false;
	}
	adjacencyCondition.notify_one();
        // adjacencyLock.unlock();
    }
    // adjacencyModifying = false;
}

void FasStableProducerConsumer::producerDepth(int depth, std::unordered_map<Variable*, std::unordered_set<Variable*>>& adjacenciesCopy) {
    for (Variable* x : nodes) {

        std::unordered_set<Variable*> adjx = adjacenciesCopy[x];

        for (Variable* y : adjx) {

            std::vector<Variable*> _adjx(adjx.begin(), adjx.end());
            _adjx.erase(std::remove(_adjx.begin(), _adjx.end(), y), _adjx.end());

            // Knowledge: possible parents
            std::vector<Variable*> ppx = _adjx;

            if (ppx.size() >= depth) {
                ChoiceGenerator cg(ppx.size(), depth);
                std::vector<int> *choice;

                for (choice = cg.next(); choice != NULL; choice = cg.next()) {
                    std::vector<Variable*> condSet = GraphUtils::asList(*choice, ppx);

		    std::sort(condSet.begin(), condSet.end());
                    
                    taskQueue.push(IndependenceTask(x, y, condSet));
                }
            }
        }
    }

    // add poison pill to stop consumers
    std::vector<Variable*> empty = {};
    IndependenceTask poisonPill(NULL, NULL, empty);
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }
}


int FasStableProducerConsumer::freeDegree() {
    int max = 0;

    for (Variable* x : nodes) {
        std::unordered_set<Variable*> opposites = adjacencies[x];

        for (Variable* y : opposites) {
            std::unordered_set<Variable*> adjx(opposites);
            adjx.erase(y);

            if (adjx.size() > max) {
                max = adjx.size();
            }
        }
    }

    return max;
}

bool FasStableProducerConsumer::searchAtDepth(int depth) {

    std::unordered_map<Variable*, std::unordered_set<Variable*>> adjacenciesCopy = adjacencies;

    std::vector<std::thread> threads;

    threads.push_back(std::thread( [&] { producerDepth(depth, adjacenciesCopy); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( [&] { consumerDepth(depth); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    return freeDegree() > depth;
}
