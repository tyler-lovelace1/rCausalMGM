#ifndef GRASP_HPP_
#define GRASP_HPP_

// [[Rcpp::depends(BH)]]

#include <boost/functional/hash.hpp>
#include <unordered_set>
#include <unordered_map>
#include <stack>
// #include "GrowShrink.hpp"
#include "GrowShrinkTree.hpp"
#include "EdgeListGraph.hpp"
#include "Knowledge.hpp"
#include <chrono>
#include <mutex>

typedef std::pair<Node, Node> NodePair;

class Grasp {

private:

    /**
     * The maximum depth for the DFS. The default it 1.
     */
    int depth = 1;

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
     * The initial graph for the for constraining the set of possible parents, or null if there is none.
     */
    EdgeListGraph *initialGraph = NULL;

    Knowledge knowledge;

    /**
     * Elapsed time of the most recent search.
     */
    long elapsedTime;

    /**
     * The number of consumer threads to create for multi-threaded steps. -1 to set automatically
     */ 
    int threads = -1;

    // TODO?
    /**
     * The logger for this class. The config needs to be set.
     */ 
    //TetradLogger logger = TetradLogger.getInstance();

    bool verbose = false;

    bool bossInit = false;

    std::list<Node> nodes;

    // GrowShrink growShrink;

    Score* scorer;

    std::map<Node, std::unique_ptr<GrowShrinkTree>> gstMap;

    std::map<Node, std::mutex> gstMutexMap;

    struct OrderGraph {
        double score=1e20, bic=1e20;
	std::list<Node> order;
	std::list<Node>::iterator start, stop;
	std::unordered_map<Node, std::unordered_set<Node>> parentMap;
	std::unordered_map<Node, double> scoreMap, bicMap;
	OrderGraph() {}
	OrderGraph(std::list<Node> order) : order(order) {
	    start = this->order.begin();
	    stop = this->order.end();
	    score = 0;
	    bic = 0;
	    for (const Node& n : this->order) scoreMap[n] = 0;
	    for (const Node& n : this->order) bicMap[n] = 0;
	}
	
	OrderGraph(const OrderGraph& other) : score(other.score),
					      bic(other.bic), 
					      parentMap(other.parentMap),
					      scoreMap(other.scoreMap),
					      bicMap(other.bicMap) {

	    // Rcpp::Rcout << "Executing copy constructor\n";

	    order = other.order;
	    
	    if (other.start == other.order.begin()) {
		start = order.begin();
	    } else {
		for (auto it = order.begin(); it != order.end(); it++) {
		    if (*it == *other.start) {
			start = it;
			break;
		    }
		}
	    }

	    if (other.stop == other.order.end()) {
		stop = order.end();
	    } else {
		for (auto it = start; it != order.end(); it++) {
		    if (*it == *other.stop) {
			stop = it;
			break;
		    }
		}
	    }
	}
	
	OrderGraph& operator=(const OrderGraph& other) {

	    // Rcpp::Rcout << "Executing copy assignment\n";
	    
	    score = other.score;
	    bic = other.bic;
	    order = other.order;
	    
	    if (other.start == other.order.begin()) {
		start = order.begin();
	    } else {
		for (auto it = order.begin(); it != order.end(); it++) {
		    if (*it == *other.start) {
			start = it;
			break;
		    }
		}
	    }

	    if (other.stop == other.order.end()) {
		stop = order.end();
	    } else {
		for (auto it = start; it != order.end(); it++) {
		    if (*it == *other.stop) {
			stop = it;
			break;
		    }
		}
	    }
	    
	    parentMap = other.parentMap;
	    scoreMap = other.scoreMap;
	    bicMap = other.bicMap;
	    return *this;
	}

	OrderGraph(OrderGraph&& other) : score(std::move(other.score)),
					 bic(std::move(other.bic)),
					 parentMap(std::move(other.parentMap)),
					 scoreMap(std::move(other.scoreMap)),
					 bicMap(std::move(other.bicMap)) {
	    // Rcpp::Rcout << "Executing move constructor\n";
	    bool flag = (other.stop == other.order.end());
	    std::swap(order, other.order);
	    start = other.start;
	    if (flag)  stop = this->order.end();
	    else       stop = other.stop;
	}
	
	OrderGraph& operator=(OrderGraph&& other) {
	    // Rcpp::Rcout << "Executing move assignment operator\n";
	    std::swap(score, other.score);
	    std::swap(bic, other.bic);
	    bool flag = (other.stop == other.order.end());
	    std::swap(order, other.order);
	    start = other.start;
	    if (flag)  stop = this->order.end();
	    else       stop = other.stop;
	    std::swap(parentMap, other.parentMap);
	    std::swap(scoreMap, other.scoreMap);
	    std::swap(bicMap, other.bicMap);
	    return *this;
	}
	
	// OrderGraph(std::vector<Node>&& order) : order(std::list<Node>(order.begin(), order.end())) {}
	bool isCovered(const Node& x, const Node& y) {
	    std::unordered_set<Node> xPar(parentMap[x]);
	    std::unordered_set<Node> yPar(parentMap[y]);

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

	std::vector<Node> getShuffledNodes() {
	    std::vector<Node> nodes(order.begin(), order.end());
	    std::random_shuffle(nodes.begin(), nodes.end(), randWrapper);
	    return nodes;
	}

	bool isValid(Knowledge& knowledge) {
	    
	    if (knowledge.isEmpty()) return true;
	    
	    for (auto it = order.begin(); it!=order.end(); it++) {
		for (auto jt = order.begin(); jt!=it; jt++) {
		    Node x = *jt;
		    Node y = *it;
		    if (knowledge.isForbidden(x, y)) {
			return false;
		    }
		}
	    }
	    return true;
	}

	void getAncestors(const Node& node, std::unordered_set<Node>& ancestors) {
	    if (ancestors.count(node)) return;

	    ancestors.insert(node);

	    for (const Node& pa : parentMap[node]) {
		getAncestors(pa, ancestors);
	    }
	}

	EdgeListGraph toDAG() const {
	    std::vector<Node> nodes(order.begin(), order.end());

	    EdgeListGraph graph(nodes);

	    for (const Node& node : nodes) {
		for (const Node& pa : parentMap.at(node)) {
		    graph.addDirectedEdge(pa, node);
		}
	    }

	    return graph;
	}

	EdgeListGraph toGraph() const {
	    std::vector<Node> nodes(order.begin(), order.end());

	    EdgeListGraph graph(nodes);

	    for (const Node& node : nodes) {
		for (const Node& pa : parentMap.at(node)) {
		    graph.addDirectedEdge(pa, node);
		}
	    }

	    return graph.getCPDAG();
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

    std::list<Node> initializeRandom();

    std::list<Node> initializeMinDegree();

    std::list<Node> initializeBOSS();

    // OrderGraph internalBOSS(Node n, OrderGraph tau);

    OrderGraph bestMove(Node n, OrderGraph tau);

    std::list<Node> minDegreeAlgorithm(EdgeListGraph graph);

public:

    /**
     * Constructs a new GRaSP search using the given independence test as oracle.
     *
     * @param independenceTest The oracle for conditional independence facts. This does not make a copy of the
     *                         independence test, for fear of duplicating the data set!
     */
    Grasp() {}

    Grasp(Score* scorer, int threads = -1);

    // ~Grasp() {
    // 	for (Node n : nodes) {
    // 	    delete gstMap[n];
    // 	}
    // }


    /**
     * @return the elapsed time of the search, in milliseconds.
     */
    long getElapsedTime() { return elapsedTime; }

    /**
     * @return the knowledge specification used in the search. Non-null.
     */
    // IKnowledge getKnowledge() { return knowledge; }

    /**
     * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
     * checked.
     *
     * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
     *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
     *              machines.
     */
    void setDepth(int depth);

    void setNumStarts(int numStarts);

    std::vector<Node> getNodes() { return graph.getNodes(); }

    void setInitialGraph(EdgeListGraph *initialGraph) { this->initialGraph = initialGraph; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setPenalty(double penalty) {
	this->penalty = penalty;
	this->scorer->setPenalty(penalty);
    }

    void setBossInit(bool bossInit) { this->bossInit = bossInit; }

    void setThreads(int threads) { this->threads = threads; }

    Knowledge getKnowledge() { return this->knowledge; }

    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }

    double update(OrderGraph& tau);

    double updateParallel(OrderGraph& tau);

    void resetTrees() {
	for (Node n : nodes) {
	    gstMap[n]->reset();
	}
    }

    // double update(OrderGraph& tau,
    // 		  std::list<Node>::iterator start,
    // 		  std::list<Node>::iterator stop);

    bool graspDFS(int tier);

    // bool graspDFSr(OrderGraph& tau, int curDepth, int tier,
    // 		   std::set<NodePair>& tucks,
    // 		   std::set<std::set<NodePair>>& dfsHistory);

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
    EdgeListGraph search(const std::vector<Node>& nodes);

};

#endif /* GRASP_HPP_ */
