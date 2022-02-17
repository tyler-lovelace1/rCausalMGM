#ifndef SEPSETSSET_HPP_
#define SEPSETSSET_HPP_

#include "SepsetProducer.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"

class SepsetsSet : public SepsetProducer  {

private:
    SepsetMap sepsets;
    IndependenceTest* test;
    double p;
    bool verbose = false;

public:
    SepsetsSet(SepsetMap sepsets, IndependenceTest *test);

    ~SepsetsSet() {}

    void fillMap() {}

    std::vector<Node> getSepset(const Node& a, const Node& b) { return *sepsets.get(a, b); }

    bool isCollider(const Node& i, const Node& j, const Node& k);

    bool isNoncollider(const Node& i, const Node& j, const Node& k);

    bool isIndependent(const Node& a, const Node& b, std::vector<Node>& c) { return test->isIndependent(a, b, c); }

    double getPValue() { return test->getPValue(); }

    double getScore() { return -(test->getPValue() - test->getAlpha()); }

    std::vector<Node> getVariables() { return test->getVariables(); }

    bool isVerbose() { return verbose; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setDepth(int depth) { this->depth = depth; }

};

#endif /* SEPSETSSET_HPP_ */
