#ifndef SEPSETPRODUCERMAJORITY_HPP_
#define SEPSETPRODUCERMAJORITY_HPP_

// [[Rcpp::depends(RcppThread)]]

#include "Triple.hpp"
#include "SepsetProducer.hpp"
#include "SepsetMap.hpp"
#include "IndependenceTest.hpp"
#include "BlockingQueue.hpp"
#include "ChoiceGenerator.hpp"
#include "DepthChoiceGenerator.hpp"
#include "RcppThread.h"
#include <thread>
#include <mutex>

/**
 * @author Tyler Lovelace
 */

class SepsetProducerMajority : public SepsetProducer {

    std::unordered_map<Triple, std::pair<int, int>> sepsetCount;

    std::unordered_set<Triple> ambiguous;

    SepsetMap extraSepsets;

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
        IndependenceTask(const Node& _a,
			 const Node& _b,
			 const Node& _c,
			 const std::vector<Node>& _s) :
	    a(_a),
	    b(_b),
	    c(_c),
	    s(_s) {} 
	
	IndependenceTask(const IndependenceTask& it) = default;
	IndependenceTask& operator=(const IndependenceTask& other) = default;
	IndependenceTask(IndependenceTask&& it) = default;
	IndependenceTask& operator=(IndependenceTask&& other) = default;
	~IndependenceTask() = default;
    };

    void producer();
    void consumer();

    const int MAX_QUEUE_SIZE = 100;
    BlockingQueue<IndependenceTask> taskQueue;

    int parallelism;
    
    std::mutex mapMutex;

public:

    SepsetProducerMajority(EdgeListGraph& _graph,
			       IndependenceTest* _test,
			       SepsetMap& _extraSepsets,
                   int threads = -1) :
	SepsetProducer(_graph, _test),
	extraSepsets(_extraSepsets),
	taskQueue(MAX_QUEUE_SIZE) {
        if (threads > 0) parallelism = threads;
        else {
            parallelism = std::thread::hardware_concurrency();
            if (parallelism == 0) {
                parallelism = 4;
                Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
            }
        }
    }

    ~SepsetProducerMajority() {}

    std::vector<Node>  getSepset(const Node& a, const Node& b);

    void fillMap();

    bool isCollider(const Node& i, const Node& j, const Node& k);

    bool isNoncollider(const Node& i, const Node& j, const Node& k);

    bool isIndependent(const Node& a, const Node& b, std::vector<Node>& c) { return test->isIndependent(a,b,c); }

    double getPValue() { return test->getPValue(); }

    double getScore() { return -(test->getPValue() - test->getAlpha()); }

    std::vector<Node> getVariables() { return test->getVariables(); }

    std::unordered_set<Triple> getAmbiguousTriples() { return ambiguous; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setDepth(int depth);
};

#endif /* SEPSETPRODUCERMAJORITY_HPP_ */
