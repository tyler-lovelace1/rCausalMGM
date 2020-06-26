#ifndef FAS_HPP_
#define FAS_HPP_

#include "Graph.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"

class Fas {

private:

    // TODO - Graph
    /**
     * The search graph. It is assumed going in that all of the true adjacencies of x are in this graph for every node
     * x. It is hoped (i.e. true in the large sample limit) that true adjacencies are never removed.
     */
    Graph graph;

    /**
     * The search nodes.
     */
    std::vector<Variable*> nodes;

    /**
     * The independence test. This should be appropriate to the types
     */
    IndependenceTest *test;

    //TODO
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
    int numIndependenceTests;

    //TODO - Logger?

    //TODO - Graph
    /**
     * The true graph, for purposes of comparison. Temporary.
     */
    Graph trueGraph;

    /**
     * The number of false dependence judgements, judged from the true graph using d-separation. Temporary.
     */
    int numFalseDependenceJudgments;

    /**
     * The number of dependence judgements. Temporary.
     */
    int numDependenceJudgement;

    int numIndependenceJudgements;

    /**
     * The sepsets found during the search.
     */
    SepsetMap sepset;

    /**
     * True if this is being run by FCI--need to skip the knowledge forbid step.
     */
    bool fci = false;

    //TODO - Graph
    /**
     * The depth 0 graph, specified initially.
     */
    Graph initialGraph;

    bool sepsetsReturnEmptyIfNotFixed;

    bool searchAtDepth0(std::vector<Variable*>& nodes, IndependenceTest *test, std::unordered_map<Variable*, std::unordered_set<Variable*>>& adjacencies);

    int freeDegree(std::vector<Variable*>& nodes, std::unordered_map<Variable*, std::unordered_set<Variable*>>& adjacencies);

    bool forbiddenEdge(Variable* x, Variable* y);

    bool searchAtDepth(std::vector<Variable*>& nodes, const IndependenceTest *test, std::unordered_map<Variable*, std::unordered_set<Variable*>>& adjacencies, int depth);

    // TODO
    // std::vector<Variable*> possibleParents(Variable* x, std::vector<Variable*>& adjx, IKnowledge knowledge);

    // TODO
    // bool possibleParentOf(std::string z, std::string x, IKnowledge knowledge);

public:

    Fas(Graph& initialGraph, IndependenceTest *test);
    Fas(IndependenceTest *test);

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
    Graph search();

    std::unordered_map<Variable*, std::unordered_set<Variable*>> searchMapOnly();

    int getDepth() { return depth; }
    void setDepth(int depth);

    // TODO
    // IKnowledge getKnowledge() { return knowledge; }
    // void setKnowledge(IKnowledge knowledge);

    int getNumIndependenceTests() { return numIndependenceTests; }

    void setTrueGraph(Graph& trueGraph) { this->trueGraph = trueGraph; }

    int getNumFalseDependenceJudgments() { return numFalseDependenceJudgments; }

    int getNumDependenceJudgments() { return numDependenceJudgement; }

    SepsetMap getSepsets() { return sepset; }

    // TODO?
    // bool isVerbose() { return verbose; }
    // void setVerbose(bool verbose) { this->verbose = verbose; }

    // TODO - Override?
    // bool isAggressivelyPreventCycles() { return false; }
    // void setAggressivelyPreventCycles(bool aggressivelyPreventCycles) {}

    // TODO - Jack has been using std::list
    std::vector<Variable*> getNodes() { return test->getVariables(); }

    int getNumIndependenceJudgements() { return numIndependenceJudgements; }

    bool isSepsetsReturnEmptyIfNotFixed() { return sepsetsReturnEmptyIfNotFixed; }

    void setSepsetsReturnEmptyIfNotFixed(bool sepsetsReturnEmptyIfNotFixed) { this->sepsetsReturnEmptyIfNotFixed = sepsetsReturnEmptyIfNotFixed; }

};

#endif /* FAS_HPP_ */