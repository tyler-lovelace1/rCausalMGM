#ifndef SEPSETPRODUCER_HPP_
#define SEPSETPRODUCER_HPP_


#include "Node.hpp"
#include "EdgeListGraph.hpp"
#include "IndependenceTest.hpp"
#include <vector>

/**
 * @author Joseph Ramsey
 */

class SepsetProducer {
    
protected:
    EdgeListGraph graph;
    IndependenceTest* test;
    int depth = -1;
    bool verbose = false;

public:

    SepsetProducer() {}

    SepsetProducer(EdgeListGraph& _graph, IndependenceTest* _test) : graph(_graph), test(_test) {}

    // virtual ~SepsetProducer() {}

    virtual void fillMap() = 0;

    virtual std::vector<Node>  getSepset(const Node& a, const Node& b) = 0;

    virtual bool isCollider(const Node& i, const Node& j, const Node& k) = 0;

    virtual bool isNoncollider(const Node& i, const Node& j, const Node& k) = 0;

    virtual bool isIndependent(const Node& a, const Node& b, std::vector<Node>& c) = 0;

    virtual double getPValue() = 0;

    virtual double getScore() = 0;

    virtual std::vector<Node> getVariables() = 0;

    virtual void setVerbose(bool verbose) = 0;

    virtual void setDepth(int depth) = 0;
};

#endif /* SEPSETPRODUCER_HPP_ */
