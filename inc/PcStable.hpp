#ifndef PCSTABLE_HPP_
#define PCSTABLE_HPP_

#include "IndependenceTest.hpp"
#include "EdgeListGraph.hpp"
#include "Variable.hpp"
#include "SepsetMap.hpp"
#include "FasStableProducerConsumer.hpp"
#include "SearchGraphUtils.hpp"
#include "MeekRules.hpp"
#include <chrono>

class PcStable {

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
     * Sepset information accumulated in the search.
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
     * True if cycles are to be aggressively prevented. May be expensive for large graphs (but also useful for large
     * graphs).
     */
    bool aggressivelyPreventCycles = false;

    // TODO?
    /**
     * The logger for this class. The config needs to be set.
     */ 
    //TetradLogger logger = TetradLogger.getInstance();

    bool verbose = false;

public:

    /**
     * Constructs a new PC search using the given independence test as oracle.
     *
     * @param independenceTest The oracle for conditional independence facts. This does not make a copy of the
     *                         independence test, for fear of duplicating the data set!
     */
    PcStable(IndependenceTest *independenceTest);

    /**
     * @return true iff edges will not be added if they would create cycles.
     */
    bool isAggressivelyPreventCycles() { return aggressivelyPreventCycles; }

    /**
     * @param aggressivelyPreventCycles Set to true just in case edges will not be addeds if they would create cycles.
     */
    void setAggressivelyPreventCycles(bool aggressivelyPreventCycles) { this->aggressivelyPreventCycles = aggressivelyPreventCycles; }

    /**
     * @return the independence test being used in the search.
     */
    IndependenceTest *getIndependenceTest() { return independenceTest; }

    /**
     * @return the elapsed time of the search, in milliseconds.
     */
    long getElapsedTime() { return elapsedTime; }

    /**
     * @return the knowledge specification used in the search. Non-null.
     */
    // IKnowledge getKnowledge() { return knowledge; }

    /**
     * @return the sepset map from the most recent search. Non-null after the first call to <code>search()</code>.
     */
    SepsetMap getSepsets() { return sepsets; }

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
     * Runs PC starting with a complete graph over all nodes of the given conditional independence test, using the given
     * independence test and knowledge and returns the resultant graph. The returned graph will be a pattern if the
     * independence information is consistent with the hypothesis that there are no latent common causes. It may,
     * however, contain cycles or bidirected edges if this assumption is not born out, either due to the actual presence
     * of latent common causes, or due to statistical errors in conditional independence judgments.
     */
    EdgeListGraph search();

    /**
     * Runs PC starting with a commplete graph over the given list of nodes, using the given independence test and
     * knowledge and returns the resultant graph. The returned graph will be a pattern if the independence information
     * is consistent with the hypothesis that there are no latent common causes. It may, however, contain cycles or
     * bidirected edges if this assumption is not born out, either due to the actual presence of latent common causes,
     * or due to statistical errors in conditional independence judgments.
     * <p>
     * All of the given nodes must be in the domain of the given conditional independence test.
     */
    EdgeListGraph search(const std::vector<Variable*>& nodes);

};

#endif /* PCSTABLE_HPP_ */