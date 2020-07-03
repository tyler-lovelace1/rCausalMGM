#ifndef GRAPHUTILS_HPP_
#define GRAPHUTILS_HPP_

#include "EdgeListGraph.hpp"

class GraphUtils {

private:

    static std::vector<Variable*> getSepsetVisit(Variable* x, Variable* y, EdgeListGraph& graph, int bound);

    static bool sepsetPathFound(Variable* a, Variable* b, Variable* y, std::unordered_set<Variable*>& path, std::vector<Variable*>& z, 
                                EdgeListGraph& graph, std::unordered_set<Triple>& colliders, int bound);

public:
    /**
     * @return true just in case there is a nonempty path from one node to
     * another. Because the path needs to be non-empty, this can distinguish
     * cycles. The case where from = to but there is no cycle from from to to
     * needs to be checked separately.
     */
    static bool existsDirectedPathFromToBreathFirst(Variable* from, Variable* to, EdgeListGraph& G);

    // Breadth first.
    static bool isDConnectedTo(Variable* x, Variable* y, std::vector<Variable*>& z, EdgeListGraph& graph);

    static std::vector<Variable*> getSepset(Variable* x, Variable* y, EdgeListGraph& graph);

    static std::vector<Variable*> getCausalOrdering(EdgeListGraph& graph);

    /**
     * @return A list of triples of the form <X, Y, Z>, where <X, Y, Z> is a
     * definite noncollider in the given graph.
     */
    static std::vector<Triple> getAmbiguousTriplesFromGraph(Variable* node, EdgeListGraph& graph);

    /**
     * @return A list of triples of the form <X, Y, Z>, where <X, Y, Z> is a
     * definite noncollider in the given graph.
     */
    static std::vector<Triple> getUnderlinedTriplesFromGraph(Variable* node, EdgeListGraph& graph);

    /**
     * @return A list of triples of the form <X, Y, Z>, where <X, Y, Z> is a
     * definite noncollider in the given graph.
     */
    static std::vector<Triple> getDottedUnderlinedTriplesFromGraph(Variable* node, EdgeListGraph& graph);

};

#endif /* GRAPHUTILS_HPP_ */