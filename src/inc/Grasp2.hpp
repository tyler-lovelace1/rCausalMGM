#ifndef GRASP_HPP_
#define GRASP_HPP_

// [[Rcpp::depends(RcppThread)]]

#include "GrowShrink.hpp"
#include "EdgeListGraph.hpp"
#include "Knowledge.hpp"
#include "RcppThread.h"
#include <chrono>

typedef std::pair<Node, Node> NodePair;

class Grasp {
private:

    /**
     * The maximum depth for the DFS. The default it 2.
     */
    int depth = 2;

    int numStarts = 3;

    double penalty = 1;

    DataSet ds;

    int parallelism;

    std::mutex scoreMutex;

    /**
     * The graph that's constructed during the search.
     */
    EdgeListGraph graph;

    /**
     * The initial graph for the for constraining the set of possible
     * parents, or null if there is none.
     */
    EdgeListGraph *initialGraph = NULL;

    Knowledge knowledge;

    /**
     * Elapsed time of the most recent search.
     */
    long elapsedTime;

    /**
     * The number of consumer threads to create for multi-threaded
     * steps. -1 to set automatically
     */ 
    int threads = -1;

    bool verbose = false;

    GrowShrink growShrink;

    struct OrderGraph {
	double score, edgeCount;
	int p;
	std::vector<Node> nodes;
	std::vector<int> order; // indices of nodes to be sorted in their topological order
	std::vector<int> updateIndices;
	std::map<Node, int> orderMap; // nodes to order idx
	std::map<Node, std::set<Node>> parentMap;
	std::map<Node, double> scoreMap;

	OrderGraph() {}

	OrderGraph(std::vector<Node> nodes) : nodes(nodes) {
	    score = 0;
	    edgeCount = 0;
	    p = nodes.size();
	    for (int i = 0; i < p; i++) {
		scoreMap[nodes.at(i)] = 0;
		orderMap[nodes.at(i)] = i;
		order.push_back(i);
	    }
	}

	bool isCovered(const Node& x, const Node& y) {
	    std::set<Node> xPar(parentMap[x]);
	    std::set<Node> yPar(parentMap[y]);

	    auto it1 = xPar.find(y);
	    auto it2 = yPar.find(x);
	    if (it1 != xPar.end()) xPar.erase(it1);
	    if (it2 != yPar.end()) yPar.erase(it2);

	    if (xPar.size() != yPar.size()) return false;

	    for (const Node& n : xPar) {
		if (yPar.count(n)==0) {
		    return false;
		}
	    }
	    
	    return true;
	}

	// DFS to determine if the edge x --> y is singular (no other
	// path from x to y)
	bool isSingular(const Node& x, const Node& y) {
	    int start = orderMap.at(x);
	    // int stop = orderMap.at(y);

	    // bool singular = true;

	    Rcpp::Rcout << "Checking if " << x << " --> " << y << " is singular\n  ";

	    std::set<Node> visited;
	    std::stack<Node> S;
	    
	    for (const Node& pa : parentMap[n]) {
		if (orderMap.at(pa) > start) {
		    S.push(pa);
		}
	    }

	    while (!S.empty()) {
		Node n = S.top();
		S.pop();
		if (orderMap.at(n) <= start) continue;

		Rcpp::Rcout << n << " ";

		visited.insert(n);

		for (Node pa : parentMap[n]) {
		    if (pa == y) {
			Rcpp::Rcout << "\n";
			return false;
		    }
		    if (visited.count(pa)==0) {
			Q.push(pa);
		    }
		}
	    }

	    Rcpp::Rcout << "\n";

	    return true;
	}

	bool isValid(Knowledge& knowledge) {   
	    if (knowledge.isEmpty()) return true;	    
	    for (int i = 1; i < order.size(); i++) {
		for (int j = 0; j < i; j++) {
		    if (knowledge.isForbidden(nodes.at(order.at(j)), nodes.at(order.at(i)))) {
			return false;
		    }
		}
	    }
	    return true;
	}

	// covered tuck operation : move y right in front of x
        void tuckCovered(const Node& x, const Node& y) {
	    int j = orderMap.at(x);
	    int k = orderMap.at(y);

	    
	    int temp = order.at(k);

	    for (int i = k; i > j; i--) {
		orderMap.at(nodes.at(order.at(i-1))) = i;
	    	order.at(i) = order.at(i-1);
	    }

	    orderMap.at(y) = j;
	    order.at(j) = temp;

	    // std::vector<int> updateIndices;
	    updateIndices.push_back(j);
	    updateIndices.push_back(k);
	    return updateIndices;
	}

	// uncovered tuck operation : move y right in front of x. Move
	// all variables in gamma in front of y's new position, in
	// order.
        void tuckUncovered(const Node& x, const Node& y) {
	    int j = orderMap.at(x);
	    int k = orderMap.at(y);
	    
	    int temp = order.at(k);

	    for (int i = k; i > j; i--) {
		orderMap.at(nodes.at(order.at(i-1))) = i;
	    	order.at(i) = order.at(i-1);
	    }

	    orderMap.at(y) = j;
	    order.at(j) = temp;

	    std::vector<int> updateIndices;
	    updateIndices.push_back(j);
	    updateIndices.push_back(k);
	    return updateIndices;
	}
	
	// BFS to identify set of ancestors of y between x and y in
	// the order graph: gamma in the GrASP paper
	std::set<Node> ancestorSetGamma(const Node& x, const Node& y) {
	    int start = orderMap.at(x);
	    // int stop = orderMap.at(y);

	    Rcpp::Rcout << "Finding set gamma : ancestors of " << y << " before " << x << " in the order graph\n  ";

	    std::set<Node> ancestors;
	    std::queue<Node> Q;

	    for (const Node& pa : parentMap[n]) {
		if (orderMap.at(pa) > start) {
		    Q.push(pa);
		}
	    }

	    while (!Q.empty()) {
		Node n = Q.front();
		Q.pop();
		if (orderMap.at(n) <= start) continue;

		Rcpp::Rcout << n << " ";

		ancestors.insert(n);

		for (const Node& pa : parentMap[n]) {
		    if (ancestors.count(pa)==0) {
			Q.push(pa);
		    }
		}
	    }

	    Rcpp::Rcout << "\n";
	    
	    return ancestors;
	}

	EdgeListGraph toGraph() const {
	    // std::vector<Node> sortedNodes;

	    // for (int idx : order) {
	    // 	sortedNodes.push_back(nodes.at(idx));
	    // }
	    
	    EdgeListGraph graph(nodes);

	    for (const Node& node : nodes) {
		for (const Node& pa : parentMap.at(node)) {
		    graph.addDirectedEdge(pa, node);
		}
	    }
	    
	    return graph.getCPDAG();
	}

	std::vector<Node> getShuffledNodes() {
	    std::vector<Node> shuffledNodes(nodes.begin(), nodes.end());
	    std::random_shuffle(shuffledNodes.begin(), shuffledNodes.end(), randWrapper);
	    return shuffledNodes;
	}
	
	bool operator==(const OrderGraph& rhs) const;
	bool operator<(const OrderGraph& rhs) const;
	bool operator!=(const OrderGraph& rhs) const;
	bool operator> (const OrderGraph& rhs) const;
	bool operator<=(const OrderGraph& rhs) const;
	bool operator>=(const OrderGraph& rhs) const;
    };

    OrderGraph pi;

    std::set<OrderGraph> bestGraphs;

    double score;

    double edgeCount;

    static int randWrapper(const int n) { return std::floor(R::unif_rand()*n); }

    std::list<Node> initializeRandom();

public:

};

#endif /* GRASP_HPP_ */
