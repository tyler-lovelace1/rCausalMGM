#ifndef SEPSETSSET_HPP_
#define SEPSETSSET_HPP_

#include "SepsetProducer.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"

class SepsetsSet : public SepsetProducer  {

private:
    SepsetMap sepsets;
    IndependenceTest *test;
    double p;
    bool verbose = false;

public:
    SepsetsSet(SepsetMap sepsets, IndependenceTest *test);

    std::vector<Variable*> getSepset(Variable* a, Variable* b) { return *sepsets.get(a, b); }

    bool isCollider(Variable* i, Variable* j, Variable* k);

    bool isNoncollider(Variable* i, Variable* j, Variable* k);

    bool isIndependent(Variable* a, Variable* b, std::vector<Variable*> c) { return test->isIndependent(a, b, c); }

    double getPValue() { return test->getPValue(); }

    double getScore() { return -(test->getPValue() - test->getAlpha()); }

    std::vector<Variable*> getVariables() { return test->getVariables(); }

    bool isVerbose() { return verbose; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

};

#endif /* SEPSETSSET_HPP_ */
