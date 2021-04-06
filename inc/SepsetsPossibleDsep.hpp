#ifndef SEPSETSPOSSIBLEDSEP_HPP_
#define SEPSETSPOSSIBLEDSEP_HPP_

#include "SepsetProducer.hpp"
#include "EdgeListGraph.hpp"
#include "IndependenceTest.hpp"
#include "Variable.hpp"
#include "ChoiceGenerator.hpp"
#include "GraphUtils.hpp"
#include <vector>
#include <algorithm>

class SepsetsPossibleDsep : public SepsetProducer  {

private:
    EdgeListGraph graph;
    IndependenceTest *test;
    int maxPathLength = 5;
    // IKnowledge knowledge = new Knowledge2();
    int depth = -1;
    bool verbose = false;

    std::vector<Variable*> getCondSet(Variable* node1, Variable* node2, int maxPathLength);

    std::unordered_set<Variable*> getPossibleDsep(Variable* x, Variable* y, int maxPathLength);

    /**
     * Removes from the list of nodes any that cannot be parents of x given the background knowledge.
     */
    std::vector<Variable*> possibleParents(Variable* x, std::vector<Variable*> nodes
                                            /*IKnowledge knowledge*/);

    // Is dependent on the implementation of background knowledge
    bool possibleParentOf(std::string _z, std::string _x /*, IKnowledge bk*/) { 
        return true /* !(bk.isForbidden(_z, _x) || bk.isRequired(_x, _z)) */;
    }

public:
    SepsetsPossibleDsep(EdgeListGraph graph, IndependenceTest *test, /*IKnowledge knowledge,*/
                                 int depth, int maxPathLength);

    /**
     * Pick out the sepset from among adj(i) or adj(k) with the highest p value.
     */
    std::vector<Variable*> getSepset(Variable* i, Variable* k);

    bool isCollider(Variable* i, Variable* j, Variable* k);

    bool isNoncollider(Variable* i, Variable* j, Variable* k);

    bool isIndependent(Variable* a, Variable* b, std::vector<Variable*> c) { return test->isIndependent(a, b, c); }

    double getPValue() { return test->getPValue(); }

    double getScore() { return -(test->getPValue() - test->getAlpha()); }

    std::vector<Variable*> getVariables() { return test->getVariables(); }

    void setVerbose(bool verbose) { this->verbose = verbose; }

};

#endif /* SEPSETSPOSSIBLEDSEP_HPP_ */
