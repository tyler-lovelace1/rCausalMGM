#ifndef FASSTABLEPRODUCERCONSUMER_HPP_
#define FASSTABLEPRODUCERCONSUMER_HPP_

// [[Rcpp::depends(BH, RcppThread)]]

#include "EdgeListGraph.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"
#include "BlockingQueue.hpp"
#include "RcppThread.h"
#include "Knowledge.hpp"
#include <boost/functional/hash.hpp>
#include <mutex>
#include <thread>
#include <atomic>

class FasStableProducerConsumer {

private:

    /**
     * The search graph. It is assumed going in that all of the true adjacencies of x are in this graph for every node
     * x. It is hoped (i.e. true in the large sample limit) that true adjacencies are never removed.
     */
    EdgeListGraph graph;

    /**
     * The search nodes.
     */
    std::vector<Node> nodes;

    // /**
    //  * Default constructor node will serve as a null node
    //  */
    // static Node nullNode;

    /**
     * The independence test. This should be appropriate to the types
     */
    IndependenceTest *test;

    //Knowledge
    /**
     * Specification of which edges are forbidden or required.
     */
    Knowledge knowledge;

    /**
     * The maximum number of variables conditioned on in any conditional independence test. If the depth is -1, it will
     * be taken to be the maximum value, which is 1000. Otherwise, it should be set to a non-negative integer.
     */
    int depth = 1000;

    /**
     * The number of independence tests.
     */
    std::atomic<int> numIndependenceTests;

    //TODO - Logger?

    /**
     * The true graph, for purposes of comparison. Temporary.
     */
    // EdgeListGraph trueGraph;

    /**
     * The number of dependence judgements. Temporary.
     */
    std::atomic<int> numDependenceJudgement;

    std::atomic<int> numIndependenceJudgements;

    /**
     * The sepsets found during the search.
     */
    SepsetMap sepset;

    /**
     * True if this is being run by FCI--need to skip the knowledge forbid step.
     */
    // bool fci = false;

    std::unordered_map<Node, std::unordered_set<Node>> adjacencies;

    std::unordered_map<std::pair<Node, Node>, double,
		       boost::hash<std::pair<Node, Node>>> edgePvals;

    std::unordered_map<std::pair<Node, Node>, std::vector<Node>,
		       boost::hash<std::pair<Node, Node>>> edgeMaxPSet;

    /**
     * The depth 0 graph, specified initially.
     */
    EdgeListGraph *initialGraph = NULL;

    bool verbose = false;

    bool sepsetsReturnEmptyIfNotFixed = false;

    bool fdr = false;

    bool searchAtDepth0();

    int freeDegree();

    bool searchAtDepth(int depth);

    // Knowlegde
    // std::vector<Node> possibleParents(Node x, std::vector<Node>& adjx, IKnowledge knowledge);
    // bool possibleParentOf(std::string z, std::string x, IKnowledge knowledge);
    // bool forbiddenEdge(Node x, Node y);

    // Concurrency
    struct IndependenceTask {
        Node x;
        Node y;
        std::vector<Node> z;
	IndependenceTask() : x(),
			     y(),
			     z() {}
        IndependenceTask(const Node& _x, const Node& _y, std::vector<Node>& _z) : x(_x),
										  y(_y),
										  z(_z) {}
	IndependenceTask(const IndependenceTask& it) = default;
	IndependenceTask& operator=(const IndependenceTask& other) = default;
	IndependenceTask(IndependenceTask&& it) = default;
	IndependenceTask& operator=(IndependenceTask&& other) = default;
	~IndependenceTask() = default;
        // IndependenceTask(const IndependenceTask& it) { x = it.x; y = it.y; z = it.z; }
    };

    const int MAX_QUEUE_SIZE = 100;
    BlockingQueue<IndependenceTask> taskQueue;
    
    int parallelism;
    
    std::mutex adjacencyMutex, pvalMutex;

    void consumerDepth0();
    void producerDepth0();

    void consumerDepth(int depth);
    void producerDepth(int depth, std::unordered_map<Node, std::unordered_set<Node>>& adjacenciesCopy);

public:

    FasStableProducerConsumer(EdgeListGraph *initialGraph, IndependenceTest *test, int threads = -1);
    FasStableProducerConsumer(IndependenceTest *test, int threads = -1);

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
    EdgeListGraph search();

    std::unordered_map<Node, std::unordered_set<Node>> searchMapOnly();

    std::vector<Node> possibleParents(const Node& x, std::vector<Node>& adjx, const Node& y);
    bool possibleParentOf(const Node& x, const Node& z);

    int getDepth() { return depth; }
    void setDepth(int depth);

    void setVerbose(bool v) { verbose = v; }

    void setFDR(bool fdr) { this->fdr = fdr; }

    // TODO
    Knowledge getKnowledge() { return knowledge; }
    void setKnowledge(Knowledge knowledge) { this->knowledge = knowledge; }

    int getNumIndependenceTests() { return numIndependenceTests; }

    // void setTrueGraph(EdgeListGraph& trueGraph) { this->trueGraph = trueGraph; }

    int getNumDependenceJudgments() { return numDependenceJudgement; }

    SepsetMap getSepsets() { return sepset; }

    // TODO - Override?
    // bool isAggressivelyPreventCycles() { return false; }
    // void setAggressivelyPreventCycles(bool aggressivelyPreventCycles) {}

    std::vector<Node> getNodes() { return test->getVariables(); }

    int getNumIndependenceJudgements() { return numIndependenceJudgements; }

    bool isSepsetsReturnEmptyIfNotFixed() { return sepsetsReturnEmptyIfNotFixed; }

    void setSepsetsReturnEmptyIfNotFixed(bool sepsetsReturnEmptyIfNotFixed) { this->sepsetsReturnEmptyIfNotFixed = sepsetsReturnEmptyIfNotFixed; }

};

#endif /* FASSTABLEPRODUCERCONSUMER_HPP_ */
