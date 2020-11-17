#ifndef PCMAX_HPP_
#define PCMAX_HPP_

#include "IndependenceTest.hpp"
#include "EdgeListGraph.hpp"
#include "Variable.hpp"
#include "FasStableProducerConsumer.hpp"
#include "OrientCollidersMaxP.hpp"
#include "MeekRules.hpp"
#include <chrono>

/**
 * Implements a modification of the the PC ("Peter/Clark") algorithm, as specified in Chapter 6 of
 * Spirtes, Glymour, and Scheines, Causation, Prediction, and Search, 2nd edition, using the rule set
 * in step D due to Chris Meek. For the modified rule set, see Chris Meek (1995), "Causal inference and
 * causal explanation with background knowledge." The modification is to replace the collider orientation
 * step by a scoring step.
 *
 * @author Joseph Ramsey.
 * @author Max Dudek - Conversion to C++
 */
class PcMax {

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
     * The maximum number of nodes conditioned on in the search. The default it 1000.
     */
    int depth = -1;

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

    bool useHeuristic = false;

    int maxPathLength = 3;

    bool verbose = false;

public:

    /**
     * Constructs a new PC search using the given independence test as oracle.
     *
     * @param independenceTest The oracle for conditional independence facts. This does not make a copy of the
     *                         independence test, for fear of duplicating the data set!
     */
    PcMax(IndependenceTest* independenceTest);

    /**
     * @return the independence test being used in the search.
     */
    IndependenceTest* getIndependenceTest() { return independenceTest; }

    /**
     * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
     * checked.
     *
     * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
     *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
     *              machines.
     */
    void setDepth(int depth);

    /**
     * @return the current depth of search--that is, the maximum number of conditioning nodes for any conditional
     * independence checked.
     */
    int getDepth() const { return depth; }

    /**
     * @return the elapsed time of the search, in milliseconds.
     */
    long getElapsedTime() const { return elapsedTime; }

    std::unordered_set<Edge> getAdjacencies();

    std::unordered_set<Edge> getNonadjacencies();

    std::vector<Variable*> getNodes() { return graph.getNodes(); }

    void setInitialGraph(EdgeListGraph *initialGraph) { this->initialGraph = initialGraph; }

    void setVerbose(bool verbose) { this->verbose = verbose; }
    bool isVerbose() { return verbose; }

    void setUseHeuristic(bool useHeuristic) { this->useHeuristic = useHeuristic; }
    bool isUseHeuristic() { return useHeuristic; }

    void setMaxPathLength(bool maxPathLength) { this->maxPathLength = maxPathLength; }
    bool getMaxPathLength() { return maxPathLength; }

    /**
     * Runs PC starting with a fully connected graph over all of the variables in the domain of the independence test.
     * See PC for caveats. The number of possible cycles and bidirected edges is far less with CPC than with PC.
     */
    EdgeListGraph search();

    EdgeListGraph search(const std::vector<Variable*>& nodes);

};

#endif /* PCMAX_HPP_ */