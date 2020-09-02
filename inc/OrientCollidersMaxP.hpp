#ifndef ORIENTCOLLIDERSMAXP_HPP_
#define ORIENTCOLLIDERSMAXP_HPP_

#include "IndependenceTest.hpp"
#include "EdgeListGraph.hpp"

/**
 * This is an optimization of the CCD (Cyclic Causal Discovery) algorithm by Thomas Richardson.
 *
 * @author Joseph Ramsey
 * @author Max Dudek - conversion to C++
 */
class OrientCollidersMaxP {

private:
    /**
     * The independence test used for the PC search.
     */
    IndependenceTest *independenceTest;

    /**
     * The maximum number of nodes conditioned on in the search. The default it 1000.
     */
    int depth = -1;

    long elapsedTime = 0;

    bool useHeuristic = true;

    int maxPathLength = 3;

    void doNode(EdgeListGraph& graph, std::unordered_map<Triple, double>& scores, Variable* b);

    void testColliderMaxP(EdgeListGraph& graph, std::unordered_map<Triple, double>& scores, Variable* a, Variable* b, Variable* c);

    void testColliderHeuristic(EdgeListGraph& graph, std::unordered_map<Triple, double>& colliders, Variable* a, Variable* b, Variable* c);

    void orientCollider(EdgeListGraph& graph, Variable* a, Variable* b, Variable* c);

    bool wouldCreateBadCollider(EdgeListGraph& graph, Variable* x, Variable* y);

    // Returns a sepset containing the nodes in 'containing' but not the nodes in 'notContaining', or
    // an empty list if there is no such sepset.
    std::vector<Variable*> sepset(EdgeListGraph& graph, Variable* a, Variable* c, std::unordered_set<Variable*>& containing, std::unordered_set<Variable*>& notContaining);

    // Returns true if there is an undirected path from x to either y or z within the given number of steps.
    bool existsShortPath(Variable* x, Variable* z, int bound, EdgeListGraph& graph);

public:
    OrientCollidersMaxP(IndependenceTest *test);

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
     * @param depth The depth of search for the Fast Adjacency Search.
     */
    int getDepth() const { return depth; }

    /**
     * @return the elapsed time of the search, in milliseconds.
     */
    long getElapsedTime() const { return elapsedTime; }

    void setUseHeuristic(bool useHeuristic) { this->useHeuristic = useHeuristic; }
    bool isUseHeuristic() { return useHeuristic; }

    void setMaxPathLength(bool maxPathLength) { this->maxPathLength = maxPathLength; }
    bool getMaxPathLength() { return maxPathLength; }

    void orient(EdgeListGraph& graph) { addColliders(graph); }

    void addColliders(EdgeListGraph& graph);


};

#endif /* ORIENTCOLLIDERSMAXP_HPP_ */