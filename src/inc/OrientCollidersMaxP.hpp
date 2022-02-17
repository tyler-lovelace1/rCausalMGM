#ifndef ORIENTCOLLIDERSMAXP_HPP_
#define ORIENTCOLLIDERSMAXP_HPP_

// [[Rcpp::depends(RcppThread)]]

#include "IndependenceTest.hpp"
#include "EdgeListGraph.hpp"
#include "BlockingQueue.hpp"
#include "ChoiceGenerator.hpp"
#include "DepthChoiceGenerator.hpp"
#include "RcppThread.h"
#include <queue>
#include <thread>
#include <mutex>

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

    EdgeListGraph *graph;

    /**
     * The maximum number of nodes conditioned on in the search. The default it 1000.
     */
    int depth = -1;

    long elapsedTime = 0;

    bool useHeuristic = false;

    int maxPathLength = 3;

    std::unordered_map<Triple, double> scores;

    std::unordered_map<Triple, bool> colliders; // True if Triple is a collider, false otherwise

    void testColliderMaxP(const Node& a, const Node& b, const Node& c);

    void testColliderHeuristic(const Node& a, const Node& b, const Node& c);

    void orientCollider(const Node& a, const Node& b, const Node& c);

    bool wouldCreateBadCollider(const Node& x, const Node& y);

    // Finds a sepset containing the nodes in 'containing' but not the nodes in 'notContaining'
    // Returns true of that set exists, false if not
    // If the set exists, it will be placed into output
    bool sepset(const Node& a, const Node& c, std::unordered_set<Node>& containing, std::unordered_set<Node>& notContaining, std::vector<Node>* output = NULL);

    // Returns true if there is an undirected path from x to either y or z within the given number of steps.
    bool existsShortPath(const Node& x, const Node& z, int bound);

    // Concurrency
    /**
     * Represents testing if a and c are independent
     * b is used to determine the corresponding Triple
     */ 
    struct IndependenceTask {
        Node a;
        Node b;
        Node c;
        std::vector<Node> s;
	IndependenceTask() {}
        IndependenceTask(const Node& _a, const Node& _b, const Node& _c,
			 const std::vector<Node>& _s) : a(_a),
							b(_b),
							c(_c),
							s(_s) {}
	
	IndependenceTask(const IndependenceTask& it) = default;
	IndependenceTask& operator=(const IndependenceTask& other) = default;
	IndependenceTask(IndependenceTask&& it) = default;
	IndependenceTask& operator=(IndependenceTask&& other) = default;
	~IndependenceTask() = default;
        
	// IndependenceTask(const IndependenceTask& it) { a = it.a; b = it.b; c = it.c; s = it.s; }
    };

    void producer();
    void consumer();

    const int MAX_QUEUE_SIZE = 100;
    BlockingQueue<IndependenceTask> taskQueue;

    int parallelism;
    
    std::mutex mapMutex;

public:
    OrientCollidersMaxP(IndependenceTest *test, EdgeListGraph *graph, int threads = -1);

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

    void setMaxPathLength(int maxPathLength) { this->maxPathLength = maxPathLength; }
    int getMaxPathLength() { return maxPathLength; }

    void orient() { addColliders(); }

    void addColliders();


};

#endif /* ORIENTCOLLIDERSMAXP_HPP_ */
