#ifndef POSSIBLEDSEPFCICONSUMERPRODUCER_HPP_
#define POSSIBLEDSEPFCICONSUMERPRODUCER_HPP_

// [[Rcpp::depends(RcppThread)]]

/**
 * This class implements the Possible-D-Sep search step of Spirtes, et al's (1993) FCI algorithm (pp 144-145).
 * Specifically, the methods in this class perform step D. of the algorithm. </p> The algorithm implemented by this
 * class is a bit broader, however, because it allows for the possibility that some pairs of variables have already been
 * compared by a different algorithm. Specifically, if the <code>prevCheck</code> variable is provided in the
 * constructor, then the algorithm pairwise checks every variable in the graph with every variable in v \
 * <code>prevCheck</code> (that is, the unchecked variables). This feature is used by the CIVI algorithm of Danks's
 * "Efficient Inclusion of Novel Variables."
 *
 * @author David Danks
 */

#include "EdgeListGraph.hpp"
#include "Knowledge.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"
#include "BlockingQueue.hpp"
#include "GraphUtils.hpp"
#include "ChoiceGenerator.hpp"
#include "RcppThread.h"
#include <limits>
#include <mutex>
#include <thread>


class PossibleDsepFciConsumerProducer {

private:

    EdgeListGraph graph;
    IndependenceTest *test;

    SepsetMap sepset;
    int depth = -1;

    /**
    * The background knowledge.
    */
    Knowledge knowledge;


    int maxReachablePathLength = -1;

    bool verbose = false;

    /**
    * Concurrency variables
    */
    int parallelism;

    // Edge poisonEdge;

    struct PossibleDsepTask {
        Edge edge;
        std::vector<Node> condSet;

        PossibleDsepTask() :
	    edge(),
	    condSet(std::vector<Node>()) {}

        PossibleDsepTask(const Edge& edge_, std::vector<Node> condSet_) :
	    edge(edge_),
	    condSet(condSet_) {}

        PossibleDsepTask(const PossibleDsepTask& it) = default;
	PossibleDsepTask& operator=(const PossibleDsepTask& other) = default;
	PossibleDsepTask(PossibleDsepTask&& it) = default;
	PossibleDsepTask& operator=(PossibleDsepTask&& other) = default;
	~PossibleDsepTask() = default;

    };

    const int MAX_QUEUE_SIZE = 100;
    BlockingQueue<PossibleDsepTask> taskQueue;
    std::mutex edgeMutex;

    void PossibleDsepProducer(std::set<Edge> edges);

    void PossibleDsepConsumer(std::unordered_map<Edge, std::vector<Node>>& edgeCondsetMap);

          //========================PRIVATE METHODS==========================//

    /**
     * Removes from the list of nodes any that cannot be parents of x given the background knowledge.
     */
    std::vector<Node> possibleParents(const Node& x, std::vector<Node>& nodes);

    bool possibleParentOf(const Node& x, const Node& z);

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
    std::unordered_set<Node> getPossibleDsep(const Node& node1,
					     const Node& node2, int maxPathLength);

public:

      //============================CONSTRUCTORS============================//

    /**
     * Creates a new SepSet and assumes that none of the variables have yet been checked.
     *
     * @param graph The GaSearchGraph on which to work
     * @param test  The IndependenceChecker to use as an oracle
     */
    PossibleDsepFciConsumerProducer(EdgeListGraph& graph, IndependenceTest *test, int threads = -1);

    PossibleDsepFciConsumerProducer(IndependenceTest *test, int threads = -1);

    ~PossibleDsepFciConsumerProducer() = default;

     //========================PUBLIC METHODS==========================//

    /**
     * Performs pairwise comparisons of each variable in the graph with the variables that have not already been
     * checked. We get the Possible-D-Sep sets for the pair of variables, and we check to see if they are independent
     * conditional on some subset of the union of Possible-D-Sep sets. This method returns the SepSet passed in the
     * constructor (if any), possibly augmented by some edge removals in this step. The GaSearchGraph passed in the
     * constructor is directly changed.
     */
    SepsetMap& search();

    void concurrentSearch(EdgeListGraph& graph, std::unordered_map<Edge, std::vector<Node>>& edgeCondsetMap);

    int getDepth() { return depth; }

    EdgeListGraph getGraph() { return graph; }

    void setDepth(int depth);

    void setVerbose(bool verbose) { this->verbose = verbose; }

    Knowledge getKnowledge() { return knowledge; }

    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }

    int getMaxReachablePathLength() { return (maxReachablePathLength == std::numeric_limits<int>::max()) ? -1 : maxReachablePathLength; }

    void setMaxPathLength(int maxReachablePathLength);

 };

 #endif /* POSSIBLEDSEPFCICONSUMERPRODUCER_HPP_ */
