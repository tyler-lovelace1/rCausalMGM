#ifndef GRAPHUTILS_HPP_
#define GRAPHUTILS_HPP_

// [[Rcpp::depends(BH, RcppThread)]]

#include <thread>

class EdgeListGraph; // Forward declaration

#include "EdgeListGraph.hpp"
#include "RcppThread.h"
#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>
#include <queue>

class GraphUtils {

private:

    static std::vector<Variable*> getSepsetVisit(Variable* x, Variable* y, EdgeListGraph& graph, int bound);

    static bool sepsetPathFound(Variable* a, Variable* b, Variable* y, std::unordered_set<Variable*>& path, std::vector<Variable*>& z,
                                EdgeListGraph& graph, std::unordered_set<Triple>& colliders, int bound);

    static void addToList(std::unordered_map<Variable*, std::vector<Variable*>> previous, Variable* b, Variable* c);

    static bool existOnePathWithPossibleParents(std::unordered_map<Variable*, std::vector<Variable*>> previous, Variable* w,
                                Variable* x, Variable* b, EdgeListGraph& graph);

    static bool existsSemidirectedPath(Variable* from, Variable* to, EdgeListGraph& G);
    
public:

    static std::vector<std::string> splitString(std::string s, const std::string& delim);

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

    /**
     * Constructs a list of nodes from the given <code>nodes</code> list at the
     * given indices in that list.
     *
     * @param indices The indices of the desired nodes in <code>nodes</code>.
     * @param nodes The list of nodes from which we select a sublist.
     * @return the The sublist selected.
     */
    static std::vector<Variable*> asList(std::vector<int>& indices, std::vector<Variable*>& nodes);

    static std::unordered_set<Variable*> asSet(std::vector<int>& indices, std::vector<Variable*>& nodes);

    static EdgeListGraph completeGraph(EdgeListGraph& graph);

    static EdgeListGraph undirectedGraph(EdgeListGraph& graph);

    static std::unordered_set<Variable*> possibleDsep(Variable* x, Variable* y, EdgeListGraph& graph, int maxPathLength);

};

#endif /* GRAPHUTILS_HPP_ */
