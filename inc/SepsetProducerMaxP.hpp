#ifndef SEPSETPRODUCERMAXP_HPP_
#define SEPSETPRODUCERMAXP_HPP_

#include "SepsetProducer.hpp"
#include "IndependenceTest.hpp"
#include "BlockingQueue.hpp"
#include "ChoiceGenerator.hpp"
#include "DepthChoiceGenerator.hpp"
#include <queue>
#include <thread>
#include <mutex>

/**
 * @author Tyler Lovelace
 */

class SepsetProducerMaxP : public SepsetProducer {

    std::unordered_map<Triple, double> scores;

    std::unordered_map<Triple, bool> colliders; // True if Triple is a collider, false otherwise

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
        IndependenceTask(Variable* _a, Variable* _b, Variable* _c, const std::vector<Variable*>& _s) :
	    a(_a),
	    b(_b),
	    c(_c),
	    s(_s) {} 
        IndependenceTask(const IndependenceTask& it) { a = it.a; b = it.b; c = it.c; s = it.s; }
    };

    void producer();
    void consumer();

    const int MAX_QUEUE_SIZE = 10000;
    BlockingQueue<IndependenceTask> taskQueue;

    int parallelism = 4;
    
    std::mutex mapMutex;

public:

    SepsetProducerMaxP(EdgeListGraph& _graph, IndependenceTest* _test) :
	SepsetProducer(_graph, _test),
	taskQueue(MAX_QUEUE_SIZE),
	parallelism(std::thread::hardware_concurrency()) {}

    std::vector<Variable*>  getSepset(Variable* a, Variable* b);

    void fillMap();

    bool isCollider(Variable* i, Variable* j, Variable* k);

    bool isNoncollider(Variable* i, Variable* j, Variable* k);

    bool isIndependent(Variable* a, Variable* b, std::vector<Variable*>  c) { return test->isIndependent(a,b,c); }

    double getPValue() { return test->getPValue(); }

    double getScore() { return -(test->getPValue() - test->getAlpha()); }

    std::vector<Variable*> getVariables() { return test->getVariables(); }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setDepth(int depth);
};

#endif /* SEPSETPRODUCERMAXP_HPP_ */
