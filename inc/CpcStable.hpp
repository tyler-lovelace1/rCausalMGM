#ifndef CPCSTABLE_HPP_
#define CPCSTABLE_HPP_

#include "IndependenceTest.hpp"
#include "EdgeListGraph.hpp"
#include "Variable.hpp"
#include "SepsetMap.hpp"
#include "FasStableProducerConsumer.hpp"
#include "SearchGraphUtils.hpp"
#include "MeekRules.hpp"
#include <chrono>
#include <thread>

/**
 * Implements a convervative version of PC, in which the Markov condition is assumed but faithfulness is tested
 * locally. Uses the PC-Stable adjacency search.
 *
 * @author Joseph Ramsey (this version).
 * @author Max Dudek (conversion to C++)
 */
class CpcStable {

private:

    /**
     * The independence test used for the PC search.
     */
    IndependenceTest *independenceTest;

    /**
     * Forbidden and required edges for the search.
     */
    // IKnowledge knowledge();

    /**
     * The sepsets
     */
    SepsetMap sepsets;

    /**
     * The maximum number of nodes conditioned on in the search. The default it 1000.
     */
    int depth = 1000;

    /**
     * The graph that's constructed during the search.
     */
    EdgeListGraph graph;

    /**
     * The initial graph for the Fast Adjacency Search, or null if there is none.
     */
    EdgeListGraph *initialGraph = NULL;

    /**
     * Elapsed time of the most recent search.
     */
    long elapsedTime;

    /**
     * The list of all unshielded triples.
     */
    std::unordered_set<Triple> allTriples;

    /**
     * True if cycles are to be aggressively prevented. May be expensive for large graphs (but also useful for large
     * graphs).
     */
    bool aggressivelyPreventCycles = false;

    /**
     * The logger for this class. The config needs to be set.
     */ 
    //TetradLogger logger = TetradLogger.getInstance();

    bool verbose = false;

    // First int in pair = # sepsets containing y
    // Second int in pair = # sepsets without y
    std::unordered_map<Triple, std::pair<int, int>> sepsetCount;

    struct ColliderTask {
        Triple t;
        std::vector<Variable*> sepset;
        ColliderTask(Triple _t, const std::vector<Variable*>& _sepset): t(_t), sepset(_sepset) {}
        ColliderTask(const ColliderTask& ct): t(ct.t), sepset(ct.sepset) {}
        ColliderTask(): t(NULL, NULL, NULL), sepset({}) {}
    };

    int parallelism = std::thread::hardware_concurrency();
    std::mutex mapMutex;
    std::condition_variable mapCondition;
    bool mapModifying = false;

    void orientUnshieldedTriples();

    bool isCollider(Triple t);

    bool isNonCollider(Triple t);

public:

    /**
     * Constructs a CPC algorithm that uses the given independence test as oracle. This does not make a copy of the
     * independence test, for fear of duplicating the data set!
     */
    CpcStable(IndependenceTest *independenceTest);

    /**
     * @return true iff edges will not be added if they would create cycles.
     */
    bool isAggressivelyPreventCycles() const { return aggressivelyPreventCycles; }

    /**
     * @param aggressivelyPreventCycles Set to true just in case edges will not be addeds if they would create cycles.
     */
    void setAggressivelyPreventCycles(bool aggressivelyPreventCycles) { this->aggressivelyPreventCycles = aggressivelyPreventCycles; }

    /**
     * @return the independence test being used in the search.
     */
    IndependenceTest *getIndependenceTest() const { return independenceTest; }

    /**
     * @return the elapsed time of the search, in milliseconds.
     */
    long getElapsedTime() const { return elapsedTime; }

    /**
     * @return the knowledge specification used in the search. Non-null.
     */
    // IKnowledge getKnowledge() { return knowledge; }

    /**
     * @return the sepset map from the most recent search. Non-null after the first call to <code>search()</code>.
     */
    SepsetMap getSepsets() { return sepsets; }

    /**
     * @return the set of all triples found during the most recent run of the algorithm. Non-null after a call to
     * <code>search()</code>.
     */
    std::unordered_set<Triple> getAllTriples() const { return allTriples; }

    std::unordered_set<Edge> getAdjacencies();

    std::unordered_set<Edge> getNonadjacencies();

    /**
     * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
     * checked.
     *
     * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
     *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
     *              machines.
     */
    void setDepth(int depth);

    std::vector<Variable*> getNodes() { return graph.getNodes(); }

    void setInitialGraph(EdgeListGraph *initialGraph) { this->initialGraph = initialGraph; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    /**
     * Runs PC starting with a fully connected graph over all of the variables in the domain of the independence test.
     * See PC for caveats. The number of possible cycles and bidirected edges is far less with CPC than with PC.
     */
    EdgeListGraph search();

    EdgeListGraph search(const std::vector<Variable*>& nodes);

    EdgeListGraph search(FasStableProducerConsumer& fas, const std::vector<Variable*>& nodes);

};

#endif /* CPCSTABLE_HPP_ */
