#ifndef SEPSETPRODUCER_HPP_
#define SEPSETPRODUCER_HPP_

// [[Rcpp::depends(RcppThread)]]

#include "Triple.hpp"
#include "SepsetProducer.hpp"
#include "Knowledge.hpp"
#include "SepsetMap.hpp"
#include "IndependenceTest.hpp"
#include "BlockingQueue.hpp"
#include "ChoiceGenerator.hpp"
#include "DepthChoiceGenerator.hpp"
#include "RcppThread.h"
#include "Node.hpp"
#include "EdgeListGraph.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include <mutex>

/**
 * @author Tyler Lovelace
 */

enum OrientRule {ORIENT_SEPSETS, ORIENT_MAJORITY, ORIENT_MAXP, ORIENT_CONSERVATIVE};

class SepsetProducer {
    
private:
    EdgeListGraph graph;
    IndependenceTest* test;

    std::unordered_map<Triple, std::pair<int, int>> sepsetCount;
    std::unordered_set<Triple> ambiguous;
    
    std::unordered_map<Triple, double> maxP;
    std::unordered_map<Triple, bool> maxPCollider;

    SepsetMap sepsets;

    Knowledge knowledge;

    OrientRule rule;

    int depth = -1;
    bool verbose = false;

    bool mapFilled = false;

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

    void producerSepsetMap();
    void consumerSepsetMap();

    static const int MAX_QUEUE_SIZE = 100;
    BlockingQueue<IndependenceTask> taskQueue;

    int parallelism = std::thread::hardware_concurrency();
    
    std::mutex mapMutex;

    std::vector<Node> possibleParents(const Node& x,
				      const std::vector<Node>& adjx,
				      const Node& y);
  
    bool possibleParentOf(const Node& x, const Node& z);
    

public:

    SepsetProducer() : taskQueue(MAX_QUEUE_SIZE) {}

    SepsetProducer(const EdgeListGraph& _graph, IndependenceTest* _test) : graph(_graph), test(_test), rule(ORIENT_MAJORITY), taskQueue(MAX_QUEUE_SIZE) {}

    SepsetProducer(const EdgeListGraph& _graph,
		   IndependenceTest* _test,
		   SepsetMap& _sepsets,
		   int threads = -1) :
	SepsetProducer(_graph, _test) {
	this->sepsets = _sepsets;
        if (threads > 0) parallelism = threads;
        else {
            parallelism = std::thread::hardware_concurrency();
            if (parallelism == 0) {
                parallelism = 4;
                Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
            }
        }
	mapFilled = false;
    }

    SepsetProducer(const EdgeListGraph& _graph, SepsetMap& sepsets, IndependenceTest *test) : taskQueue(MAX_QUEUE_SIZE){
	this->graph = _graph;
	this->sepsets = sepsets;
	this->test = test;
	this->rule = ORIENT_SEPSETS;
    }

    // SepsetProducer(EdgeListGraph& _graph, SepsetMap& sepsets,
    // 		   IndependenceTest *test) : taskQueue(MAX_QUEUE_SIZE) {
    // 	this->graph = _graph;
    // 	this->sepsets = sepsets;
    // 	this->test = test;
    // 	this->rule = ORIENT_SEPSETS;
    // }

    // SepsetProducer(const SepsetProducer& other) = default;

    SepsetProducer(const SepsetProducer& other) : graph(other.graph),
    						  test(other.test),
    						  sepsetCount(other.sepsetCount),
    						  ambiguous(other.ambiguous),
    						  maxP(other.maxP),
    						  maxPCollider(other.maxPCollider),
    						  sepsets(other.sepsets),
    						  knowledge(other.knowledge),
    						  rule(other.rule),
    						  depth(other.depth),
    						  mapFilled(other.mapFilled),
    						  verbose(other.verbose),
    						  parallelism(other.parallelism),
    						  taskQueue(other.taskQueue) {}
    
    SepsetProducer(SepsetProducer&& other) = default;

    // SepsetProducer& operator=(const SepsetProducer& other) = default;

    SepsetProducer& operator=(const SepsetProducer& other) {
    	graph = other.graph;
    	test = other.test;
    	sepsetCount = other.sepsetCount;
    	ambiguous = other.ambiguous;
    	maxP = other.maxP;
    	maxPCollider = other.maxPCollider;
    	sepsets = other.sepsets;
    	knowledge = other.knowledge;
    	rule = other.rule;
    	depth = other.depth;
    	mapFilled = other.mapFilled;
    	verbose = other.verbose;
    	parallelism = other.parallelism;
    	taskQueue = other.taskQueue;
    	return *this;
    }
    
    SepsetProducer& operator=(SepsetProducer&& other) = default;

    ~SepsetProducer() = default;

    std::vector<Node> getSepset(const Node& a, const Node& b);

    // std::vector<Node> getMinSepset(const Node& a, const Node& b, std::vector<Node>& sepset);

    void fillMap();

    bool isCollider(const Node& i, const Node& j, const Node& k);

    bool isNoncollider(const Node& i, const Node& j, const Node& k);

    bool isIndependent(const Node& a, const Node& b, std::vector<Node>& c) { return test->isIndependent(a,b,c); }

    double getPValue() { return test->getPValue(); }

    double getScore() { return -(test->getPValue() - test->getAlpha()); }

    std::vector<Node> getVariables() { return test->getVariables(); }

    std::vector<Triple> getOrderedColliders();

    std::unordered_set<Triple> getAmbiguousTriples() { return ambiguous; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setDepth(int depth);

    OrientRule getOrientRule() { return rule; }
    void setOrientRule(OrientRule rule) {this->rule = rule; }

    Knowledge getKnowledge() { return knowledge; }
    void setKnowledge(Knowledge knowledge) { this->knowledge = knowledge; }

    static OrientRule str2rule(std::string rule) {
	std::transform(rule.begin(), rule.end(), rule.begin(),
		       [](unsigned char c){ return std::tolower(c); });
	if (rule == "majority") {
	    return ORIENT_MAJORITY;
	} else if (rule == "conservative") {
	    return ORIENT_CONSERVATIVE;
	} else if (rule == "maxp") {
	    return ORIENT_MAXP;
	} else if (rule == "sepsets") {
	    return ORIENT_SEPSETS;
	} else {
	    throw std::invalid_argument("Orientation rule must be one of {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
	}
    }

    static std::string rule2str(OrientRule rule) {
	if (rule == ORIENT_MAJORITY) {
	    return "majority";
	} else if (rule == ORIENT_CONSERVATIVE) {
	    return "conservative";
	} else if (rule == ORIENT_MAXP) {
	    return "maxp";
	} else if (rule == ORIENT_SEPSETS) {
	    return "sepsets";
	} else {
	    throw std::invalid_argument("Invalid OrientRule value");
	}
    }


    // virtual ~SepsetProducer() {}

    // virtual void fillMap() = 0;

    // virtual std::vector<Node>  getSepset(const Node& a, const Node& b) = 0;

    // virtual bool isCollider(const Node& i, const Node& j, const Node& k) = 0;

    // virtual bool isNoncollider(const Node& i, const Node& j, const Node& k) = 0;

    // virtual bool isIndependent(const Node& a, const Node& b, std::vector<Node>& c) = 0;

    // virtual double getPValue() = 0;

    // virtual double getScore() = 0;

    // virtual std::vector<Node> getVariables() = 0;

    // virtual void setVerbose(bool verbose) = 0;

    // virtual void setDepth(int depth) = 0;
};

#endif /* SEPSETPRODUCER_HPP_ */
