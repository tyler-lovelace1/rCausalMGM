#ifndef ORIENTCOLLIDERSMAXP_HPP_
#define ORIENTCOLLIDERSMAXP_HPP_

#include "IndependenceTest.hpp"
#include "EdgeListGraph.hpp"
#include "BlockingQueue.hpp"
#include "ChoiceGenerator.hpp"
#include "DepthChoiceGenerator.hpp"
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

    void testColliderMaxP(Variable* a, Variable* b, Variable* c);

    void testColliderHeuristic(Variable* a, Variable* b, Variable* c);

    void orientCollider(Variable* a, Variable* b, Variable* c);

    bool wouldCreateBadCollider(Variable* x, Variable* y);

    // Finds a sepset containing the nodes in 'containing' but not the nodes in 'notContaining'
    // Returns true of that set exists, false if not
    // If the set exists, it will be placed into output
    bool sepset(Variable* a, Variable* c, std::unordered_set<Variable*>& containing, std::unordered_set<Variable*>& notContaining, std::vector<Variable*>* output = NULL);

    // Returns true if there is an undirected path from x to either y or z within the given number of steps.
    bool existsShortPath(Variable* x, Variable* z, int bound);

    // Concurrency
    /**
     * Represents testing if a and c are independent
     * b is used to determine the corresponding Triple
     */ 
    struct IndependenceTask {
        Variable* a;
        Variable* b;
        Variable* c;
        std::vector<Variable*> s;
	    IndependenceTask() {}
        IndependenceTask(Variable* _a, Variable* _b, Variable* _c, const std::vector<Variable*>& _s) : a(_a),
										   b(_b),
										   c(_c),
                                           s(_s) {} 
        IndependenceTask(const IndependenceTask& it) { a = it.a; b = it.b; c = it.c; s = it.s; }
    };

    void producer();
    void consumer();

    const int MAX_QUEUE_SIZE = 10000;
    BlockingQueue<IndependenceTask> taskQueue;

    int parallelism = std::thread::hardware_concurrency();
    
    std::mutex mapMutex;

public:
    OrientCollidersMaxP(IndependenceTest *test, EdgeListGraph *graph);

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