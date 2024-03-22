#ifndef GRAPHUTILS_HPP_
#define GRAPHUTILS_HPP_

// [[Rcpp::depends(BH, RcppThread)]]

// #include <thread>

class EdgeListGraph; // Forward declaration

#include "EdgeListGraph.hpp"
#include "SepsetMap.hpp"
#include "RcppThread.h"
#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>
#include <queue>
#include <stack>
#include <list>
#include <set>

class GraphUtils {

private:

  // static std::vector<Node> getSepsetVisit(Node x, Node y, EdgeListGraph& graph, int bound);

    static bool sepsetPathFound(Node a, Node b, Node y,
				std::unordered_set<Node>& path,
				std::vector<Node>& z, EdgeListGraph& graph,
				std::unordered_set<Triple>& colliders, int bound);

    static void addToList(std::unordered_map<Node, std::vector<Node>> previous,
			  Node b, Node c);

    static bool existOnePathWithPossibleParents(std::unordered_map<Node, std::vector<Node>> previous, Node w, Node x, Node b, EdgeListGraph& graph);

    static bool existOnePathWithPossibleParents(std::map<Node,std::set<Node>>& previous,
						EdgeListGraph& graph, const Node& w,
						const Node& x, const Node& b);

    static bool existsSemidirectedPath(Node from, Node to, EdgeListGraph& G);
    
public:

    static std::vector<std::string> splitString(std::string s, const std::string& delim);

    /**
     * @return true just in case there is a nonempty path from one node to
     * another. Because the path needs to be non-empty, this can distinguish
     * cycles. The case where from = to but there is no cycle from from to to
     * needs to be checked separately.
     */
    static bool existsDirectedPathFromToBreathFirst(Node from, Node to, EdgeListGraph& G);

    // Breadth first.
    static bool isDConnectedTo(Node x, Node y, std::vector<Node>& z, EdgeListGraph& graph);

    static std::vector<Node> getSepset(const Node& x, const Node& y, EdgeListGraph& graph);

    static bool sepsetPathFound(const Node& a, const Node& b, const Node& y, std::set<Node>& path, std::vector<Node>& z, std::set<Triple>& colliders, EdgeListGraph& graph);

    static std::vector<Node> getPassNodes(const Node& a, const Node& b,
					  std::vector<Node>& z,
					  EdgeListGraph& graph);

    static bool reachable(const Node& a, const Node& b, const Node& c,
			  std::vector<Node>& z, EdgeListGraph& graph);

    static bool isAncestor(const Node& b, std::vector<Node>& z, EdgeListGraph& graph);

    static std::vector<Node> getInducingPath(const Node& x, const Node& y,
					     EdgeListGraph& graph);

    static std::vector<Node> getCausalOrdering(EdgeListGraph& graph);

    /**
     * @return A list of triples of the form <X, Y, Z>, where <X, Y, Z> is a
     * definite noncollider in the given graph.
     */
    static std::vector<Triple> getAmbiguousTriplesFromGraph(Node node, EdgeListGraph& graph);

    /**
     * @return A list of triples of the form <X, Y, Z>, where <X, Y, Z> is a
     * definite noncollider in the given graph.
     */
    static std::vector<Triple> getUnderlinedTriplesFromGraph(Node node, EdgeListGraph& graph);

    /**
     * @return A list of triples of the form <X, Y, Z>, where <X, Y, Z> is a
     * definite noncollider in the given graph.
     */
    static std::vector<Triple> getDottedUnderlinedTriplesFromGraph(Node node, EdgeListGraph& graph);

    /**
     * Constructs a list of nodes from the given <code>nodes</code> list at the
     * given indices in that list.
     *
     * @param indices The indices of the desired nodes in <code>nodes</code>.
     * @param nodes The list of nodes from which we select a sublist.
     * @return the The sublist selected.
     */
    static std::vector<Node> asList(std::vector<int>& indices, std::vector<Node>& nodes);

    static std::unordered_set<Node> asSet(std::vector<int>& indices, std::vector<Node>& nodes);

    static EdgeListGraph completeGraph(EdgeListGraph& graph);

    static EdgeListGraph undirectedGraph(EdgeListGraph& graph);

    static std::unordered_set<Node> possibleDsep(Node x, Node y, EdgeListGraph& graph, int maxPathLength);

    static std::set<Node> possibleDsep2(const Node& x, const Node& y, EdgeListGraph& graph, int maxPathLength);

    static bool existsPossibleColliderPath(Node from, Node to, std::unordered_set<Node> Z, EdgeListGraph& G);

};

#endif /* GRAPHUTILS_HPP_ */
