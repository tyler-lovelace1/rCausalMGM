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
    this->poisonEdge = Edge(new ContinuousVariable("EJWMX3RCpPi0qbp"),
			    new ContinuousVariable("nLtWU7DmeZyYPZs"),
			    ENDPOINT_NULL, ENDPOINT_NULL);
}

//========================PRIVATE METHODS==========================//


//This method is based on background knowledge aka IKnowledge which is not yet implemented
/**
* Removes from the list of nodes any that cannot be parents of x given the background knowledge.
*/
std::vector<Variable*> PossibleDsepFciConsumerProducer::possibleParents(Variable* x, std::vector<Variable*> nodes /*, IKnowledge knowledge*/) {

    std::vector<Variable*> possibleParents;
    std::string x_ = x->getName();

    for (Variable* z : nodes) {
        std::string z_ = z->getName();

        if (possibleParentOf(z_, x_ /* , knowledge*/)) {
            possibleParents.push_back(z);
        }
    }
    return possibleParents;
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
std::unordered_set<Variable*> PossibleDsepFciConsumerProducer::getPossibleDsep(Variable* node1, Variable* node2, int maxPathLength) {

    std::unordered_set<Variable*>  dsep = GraphUtils::possibleDsep(node1, node2, graph, maxPathLength);

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
    std::unordered_map<Edge, std::vector<Variable*>> edgeCondsetMap;

    concurrentSearch(graph, edgeCondsetMap);

    for (std::pair<Edge, std::vector<Variable*>> entry : edgeCondsetMap) {
        Edge edge = entry.first;
        std::vector<Variable*> condSet = entry.second;
        Variable* x = edge.getNode1();
        Variable* y = edge.getNode2();

        graph.removeEdge(x, y);
        if (verbose) Rcpp::Rcout << "Removed " << x->getName() << " --- " << y->getName() << "\n";
        sepset.set(x, y, condSet);
    }
    return sepset;
}

void PossibleDsepFciConsumerProducer::concurrentSearch(EdgeListGraph& graph, std::unordered_map<Edge, std::vector<Variable*>>& edgeCondsetMap) {
    const std::unordered_set<Edge> edges(graph.getEdges());
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

void PossibleDsepFciConsumerProducer::PossibleDsepProducer(std::unordered_set<Edge> edges) {
    PossibleDsepTask poisonPill(poisonEdge, std::vector<Variable*>());

    for (Edge edge : edges) {
        Variable* x = edge.getNode1();
        Variable* y = edge.getNode2();

        std::unordered_set<Variable*> possibleDsepSet = getPossibleDsep(x, y, maxReachablePathLength);
        std::vector<Variable*> possibleDsep;

        possibleDsep.insert(possibleDsep.end(), possibleDsepSet.begin(), possibleDsepSet.end());

        // bool noEdgeRequired = getKnowledge().noEdgeRequired(x.getName(), y.getName());
        bool noEdgeRequired = true;

        if (!noEdgeRequired)
            continue;

        //possible parents is not fully implemented without background knowledge
        std::vector<Variable*> possParents = possibleParents(x, possibleDsep /*, getKnowledge() */ );
        int depth_ = getDepth() == -1 ? 1000 : getDepth();

        for (int d = 0; d <= std::min((std::size_t) depth_, possParents.size()); d++) {
            ChoiceGenerator cg(possParents.size(), d);
            std::vector<int> *combination;

            for (combination = cg.next(); combination != NULL; combination = cg.next()) {
                std::vector<Variable*> condSet = GraphUtils::asList(*combination, possParents);
                PossibleDsepTask newTask(edge, condSet);
                taskQueue.push(newTask);
            }
        }

        possibleDsepSet = getPossibleDsep(y, x, maxReachablePathLength);
        possibleDsep.clear();
        possibleDsep.insert(possibleDsep.end(), possibleDsepSet.begin(), possibleDsepSet.end());

        //possible parents is not fully implemented without background knowledge
        possParents = possibleParents(y, possibleDsep /*, getKnowledge() */);

        for (int d = 0; d <= std::min((std::size_t) depth_, possParents.size()); d++) {
            ChoiceGenerator cg (possParents.size(), d);
            std::vector<int> *combination;

            for (combination = cg.next(); combination != NULL; combination = cg.next()) {
                std::vector<Variable*> condSet = GraphUtils::asList(*combination, possParents);
                PossibleDsepTask newTask(edge, condSet);
                taskQueue.push(newTask);
            }
        }

	if (RcppThread::isInterrupted()) {
	  break;
	}
    }

    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }
    return;
}

void PossibleDsepFciConsumerProducer::PossibleDsepConsumer(std::unordered_map<Edge, std::vector<Variable*>>& edgeCondsetMap) {
    PossibleDsepTask task = taskQueue.pop();
    while (task.edge.getNode1()->getName() != "EJWMX3RCpPi0qbp"
	   && task.edge.getNode2()->getName() != "nLtWU7DmeZyYPZs") {
	
        if (edgeCondsetMap.count(task.edge) == 0) {
	    Variable* x = task.edge.getNode1();
	    Variable* y = task.edge.getNode2();
	    if (x->getName() > y->getName()) {
		x = task.edge.getNode2();
		y = task.edge.getNode1();
	    }
	    
            if (test->isIndependent(x, y, task.condSet)) {
                std::lock_guard<std::mutex> edgeLock(edgeMutex);
		edgeCondsetMap[task.edge] = task.condSet;
                // edgeCondsetMap.insert(std::pair<Edge,
		// 		      std::vector<Variable*>>(task.edge, task.condSet));
	    }
        }
        task = taskQueue.pop();
    }
    return;
}
