#ifndef FASSTABLEPRODUCERCONSUMER_HPP_
#define FASSTABLEPRODUCERCONSUMER_HPP_

#include "EdgeListGraph.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"
#include "ChoiceGenerator.hpp"
#include "IndependenceTask.hpp"
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <boost/circular_buffer.hpp>

using TaskQueue = std::queue<IndependenceTask, boost::circular_buffer<IndependenceTask>>;

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
    std::vector<Variable*> nodes;

    /**
     * The independence test. This should be appropriate to the types
     */
    IndependenceTest *test;

    //Knowledge
    /**
     * Specification of which edges are forbidden or required.
     */
    // private IKnowledge knowledge = new Knowledge2();

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

    std::unordered_map<Variable*, std::unordered_set<Variable*>> adjacencies;

    /**
     * The depth 0 graph, specified initially.
     */
    EdgeListGraph *initialGraph = NULL;

    bool sepsetsReturnEmptyIfNotFixed = true;

    bool searchAtDepth0();

    int freeDegree();

    bool searchAtDepth(int depth);

    // Knowlegde
    // std::vector<Variable*> possibleParents(Variable* x, std::vector<Variable*>& adjx, IKnowledge knowledge);
    // bool possibleParentOf(std::string z, std::string x, IKnowledge knowledge);
    // bool forbiddenEdge(Variable* x, Variable* y);

    // Concurrency
    const int MAX_QUEUE_SIZE = 10;
    TaskQueue taskQueue;
    
    const int parallelism = 1; //TODO - change to number of processors
    
    std::mutex queueMutex;
    std::mutex adjacencyMutex;

    std::condition_variable producerGo;
    std::condition_variable consumerGo;

    std::atomic<bool> stillProducing;

    void consumerDepth0();
    void producerDepth0();

public:

    FasStableProducerConsumer(EdgeListGraph *initialGraph, IndependenceTest *test);
    FasStableProducerConsumer(IndependenceTest *test);

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

    std::unordered_map<Variable*, std::unordered_set<Variable*>> searchMapOnly();

    int getDepth() { return depth; }
    void setDepth(int depth);

    // TODO
    // IKnowledge getKnowledge() { return knowledge; }
    // void setKnowledge(IKnowledge knowledge);

    int getNumIndependenceTests() { return numIndependenceTests; }

    // void setTrueGraph(EdgeListGraph& trueGraph) { this->trueGraph = trueGraph; }

    int getNumDependenceJudgments() { return numDependenceJudgement; }

    SepsetMap getSepsets() { return sepset; }

    // TODO - Override?
    // bool isAggressivelyPreventCycles() { return false; }
    // void setAggressivelyPreventCycles(bool aggressivelyPreventCycles) {}

    std::vector<Variable*> getNodes() { return test->getVariables(); }

    int getNumIndependenceJudgements() { return numIndependenceJudgements; }

    bool isSepsetsReturnEmptyIfNotFixed() { return sepsetsReturnEmptyIfNotFixed; }

    void setSepsetsReturnEmptyIfNotFixed(bool sepsetsReturnEmptyIfNotFixed) { this->sepsetsReturnEmptyIfNotFixed = sepsetsReturnEmptyIfNotFixed; }

};

#endif /* FASSTABLEPRODUCERCONSUMER_HPP_ */