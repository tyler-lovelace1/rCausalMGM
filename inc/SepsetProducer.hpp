#ifndef SEPSETPRODUCER_HPP_
#define SEPSETPRODUCER_HPP_


#include "Variable.hpp"
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

    virtual void fillMap() = 0;

    virtual std::vector<Variable*>  getSepset(Variable* a, Variable* b) = 0;

    virtual bool isCollider(Variable* i, Variable* j, Variable* k) = 0;

    virtual bool isNoncollider(Variable* i, Variable* j, Variable* k) = 0;

    virtual bool isIndependent(Variable* a, Variable* b, std::vector<Variable*>  c) = 0;

    virtual double getPValue() = 0;

    virtual double getScore() = 0;

    virtual std::vector<Variable*> getVariables() = 0;

    virtual void setVerbose(bool verbose) = 0;

    virtual void setDepth(int depth) = 0;
};

#endif /* SEPSETPRODUCER_HPP_ */
