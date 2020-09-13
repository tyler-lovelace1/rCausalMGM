#include "FasStableProducerConsumer.hpp"

FasStableProducerConsumer::FasStableProducerConsumer(EdgeListGraph *initialGraph, IndependenceTest *test) {
    this->initialGraph = initialGraph;
    this->test = test;
    this->nodes = test->getVariables();
}

FasStableProducerConsumer::FasStableProducerConsumer(IndependenceTest *test) {
    this->test = test;
    this->nodes = test->getVariables();
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
    Rcpp::Rcout << "Starting FasStableProducerConsumer Adjacency Search." << std::endl;

    sepset = SepsetMap();
    sepset.setReturnEmptyIfNotSet(sepsetsReturnEmptyIfNotFixed);

    int _depth = depth;

    if (_depth == -1) _depth = 1000;

    adjacencies.clear();

    for (Variable* node : nodes) {
        adjacencies[node] = {};
    }

    taskQueue = TaskQueue(boost::circular_buffer<IndependenceTask>(MAX_QUEUE_SIZE));

    for (int d = 0; d <= _depth; d++) {
        bool more;

        if (d == 0) {
            more = searchAtDepth0();
        } else {
            more = searchAtDepth(d);
        }

        if (!more) break;
    }

    graph = EdgeListGraph(nodes);

    for (int i = 0; i < nodes.size(); i++) {
        for (int j = i+1; j < nodes.size(); j++) {
            Variable* x = nodes[i];
            Variable* y = nodes[j];

            if (adjacencies[x].count(y)) {
                graph.addUndirectedEdge(x, y);
            }
        }
    }

    Rcpp::Rcout << "Finishing FasStableProducerConsumer Adjacency Search." << std::endl;

    Rcpp::Rcout << "Fas graph: \n" << graph << std::endl;

    return graph;
}

std::unordered_map<Variable*, std::unordered_set<Variable*>> FasStableProducerConsumer::searchMapOnly() {
    Rcpp::Rcout << "Starting FasStableProducerConsumer Adjacency Search." << std::endl;

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

        if (d == 0) {
            more = searchAtDepth0();
        } else {
            more = searchAtDepth(d);
        }

        if (!more) break;
    }

    Rcpp::Rcout << "Finishing FasStableProducerConsumer Adjacency Search." << std::endl;
    return adjacencies;
}

bool FasStableProducerConsumer::searchAtDepth0() {
    stillProducing = true;

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

            std::unique_lock<std::mutex> queueLock(queueMutex);
            producerGo.wait(queueLock, [&]() { return taskQueue.size() < MAX_QUEUE_SIZE; }); // Wait for space in buffer
            taskQueue.push(IndependenceTask(x, y, empty));
            queueLock.unlock();
            consumerGo.notify_one();
        }
    }

    stillProducing = false;
    consumerGo.notify_all(); // Tell all consumers to end

}

void FasStableProducerConsumer::consumerDepth0() {
    while(true) {
        std::unique_lock<std::mutex> queueLock(queueMutex);
        consumerGo.wait(queueLock, [&]() { return !taskQueue.empty() || !stillProducing; });
        if (taskQueue.empty()) return; // If the queue is empty at this point, then we must be done procuding
        IndependenceTask task = taskQueue.front();
        taskQueue.pop();
        queueLock.unlock();
        producerGo.notify_one();

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

        std::unique_lock<std::mutex> adjacencyLock(adjacencyMutex);
        if (independent && noEdgeRequired) {
            if (!sepset.isReturnEmptyIfNotSet()) {
                sepset.set(task.x, task.y, task.z);
            }
        } else if (!forbiddenEdge) {
            adjacencies[task.x].insert(task.y);
            adjacencies[task.y].insert(task.x);
        } 
        adjacencyLock.unlock();
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
    int count = 0;

    std::unordered_map<Variable*, std::unordered_set<Variable*>> adjacenciesCopy = adjacencies;

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
                    
                    numIndependenceTests++;
                    bool independent;
                    independent = test->isIndependent(x, y, condSet);

                    if (independent) {
                        numIndependenceJudgements++;
                    } else {
                        numDependenceJudgement++;
                    }

                    // Knowledge
                    bool noEdgeRequired = true;

                    if (independent && noEdgeRequired) {
                        // Rcpp::Rcout << "x = " << x->getName() << " y = " << y->getName() << std::endl;
                        // Rcpp::Rcout << "adjacencies[x] BEFORE = { ";
                        // for (Variable* v : adjacencies[x]) { Rcpp::Rcout << v->getName() << " "; }
                        // Rcpp::Rcout << "}" << std::endl;
                        adjacencies[x].erase(y);
                        adjacencies[y].erase(x);
                        // Rcpp::Rcout << "adjacencies[x] AFTER = { ";
                        // for (Variable* v : adjacencies[x]) { Rcpp::Rcout << v->getName() << " "; }
                        // Rcpp::Rcout << "}" << std::endl;

                        sepset.set(x, y, condSet);

                        goto EDGE_CONTINUE; // No need to test other combinations
                    }
                }
            }
            EDGE_CONTINUE:;
        }
    }

    return freeDegree() > depth;
}
