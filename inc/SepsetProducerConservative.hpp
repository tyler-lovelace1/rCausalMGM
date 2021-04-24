#ifndef SEPSETPRODUCERCONSERVATIVE_HPP_
#define SEPSETPRODUCERCONSERVATIVE_HPP_

#include "Triple.hpp"
#include "SepsetProducer.hpp"
#include "SepsetMap.hpp"
#include "IndependenceTest.hpp"
#include "BlockingQueue.hpp"
#include "ChoiceGenerator.hpp"
#include "DepthChoiceGenerator.hpp"
#include <thread>
#include <mutex>

/**
 * @author Tyler Lovelace
 */

class SepsetProducerConservative : public SepsetProducer {

    std::unordered_map<Triple, std::pair<int, int>> sepsetCount;

    std::unordered_set<Triple> ambiguous;

    SepsetMap extraSepsets;

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

    const int MAX_QUEUE_SIZE = 100;
    BlockingQueue<IndependenceTask> taskQueue;

    int parallelism = 4;
    
    std::mutex mapMutex;
    std::condition_variable mapCondition;
    bool mapModifying = false;

public:

    SepsetProducerConservative(EdgeListGraph& _graph,
			       IndependenceTest* _test,
			       SepsetMap& _extraSepsets) :
	SepsetProducer(_graph, _test),
	extraSepsets(_extraSepsets),
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

    std::unordered_set<Triple> getAmbiguousTriples() { return ambiguous; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setDepth(int depth);
};

#endif /* SEPSETPRODUCERCONSERVATIVE_HPP_ */
