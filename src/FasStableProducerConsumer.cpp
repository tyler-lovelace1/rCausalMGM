#include "FasStableProducerConsumer.hpp"

FasStableProducerConsumer::FasStableProducerConsumer(EdgeListGraph *initialGraph, IndependenceTest *test) : FasStableProducerConsumer(test) 
{
    this->initialGraph = initialGraph;
}

FasStableProducerConsumer::FasStableProducerConsumer(IndependenceTest *test) : taskQueue(MAX_QUEUE_SIZE) 
{
    this->test = test;
    this->nodes = test->getVariables();

    Rcpp::Rcout << "paralellism = " << parallelism << std::endl;
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
    Rcpp::Rcout << "Starting FasStableProducerConsumer Adjacency Search." << std::endl;

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

void FasStableProducerConsumer::consumerDepth(int depth) {
    while(true) {
        IndependenceTask task = taskQueue.pop();

        //If poison, return
        if (task.x == NULL && task.y == NULL) return;

        std::unique_lock<std::mutex> adjacencyLock(adjacencyMutex);
        bool edgeExists = adjacencies[task.x].count(task.y) && adjacencies[task.y].count(task.x);
        adjacencyLock.unlock();

        if (!edgeExists) continue; // Skip if the edge no longer exists

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

        adjacencyLock.lock();
        if (independent && noEdgeRequired) {
            adjacencies[task.x].erase(task.y);
            adjacencies[task.y].erase(task.x);
            sepset.set(task.x, task.y, task.z);
        }
        adjacencyLock.unlock();
    }
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
