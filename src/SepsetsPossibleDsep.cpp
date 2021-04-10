#include "SepsetsPossibleDsep.hpp"

SepsetsPossibleDsep::SepsetsPossibleDsep(EdgeListGraph& graph, IndependenceTest *test, /*IKnowledge knowledge,*/
					 int depth, int maxPathLength) : SepsetProducer(graph, test) {
    // this->graph = graph;
    // this->test = test;
    this->maxPathLength = maxPathLength;
    // this.knowledge = knowledge;
    this->depth = depth;
}

std::vector<Variable*> SepsetsPossibleDsep::getSepset(Variable* i, Variable* k) {
    std::vector<Variable*> condSet = getCondSet(i, k, maxPathLength);

    if (condSet.empty()) {
        condSet = getCondSet(k, i, maxPathLength);
    }

    return condSet;
}

bool SepsetsPossibleDsep::isCollider(Variable* i, Variable* j, Variable* k) {
    std::vector<Variable*> sepset = getSepset(i, k);
    return (!sepset.empty())  && !(std::count(sepset.begin(), sepset.end(), j) == 1);
}

bool SepsetsPossibleDsep::isNoncollider(Variable* i, Variable* j, Variable* k) {
    std::vector<Variable*> sepset = getSepset(i, k);
    return (!sepset.empty()) && (std::count(sepset.begin(), sepset.end(), j) == 1);
}

std::vector<Variable*> SepsetsPossibleDsep::getCondSet(Variable* node1, Variable* node2, int maxPathLength) {
    std::unordered_set<Variable*> possibleDsepSet = getPossibleDsep(node1, node2, maxPathLength);
    std::vector<Variable*> possibleDsep (possibleDsepSet.begin(), possibleDsepSet.end());

    // bool noEdgeRequired = knowledge.noEdgeRequired(node1.getName(), node2.getName());

    std::vector<Variable*> possParents = possibleParents(node1, possibleDsep /*, knowledge*/);

    size_t depth_ = depth == -1 ? 1000 : depth;

    for (int d = 0; d <= std::min(depth_, possParents.size()); d++) {


        ChoiceGenerator cg (possParents.size(), d);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            std::vector<Variable*> condSet = GraphUtils::asList(*choice, possParents);
            bool independent = test->isIndependent(node1, node2, condSet);

            if (independent /*&& noEdgeRequired*/) {
                return condSet;
            }
        }
    }
    std::vector<Variable*> nullSet;
    return nullSet;
}

std::unordered_set<Variable*> SepsetsPossibleDsep::getPossibleDsep(Variable* x, Variable* y, int maxPathLength) {
    std::unordered_set<Variable*> dsep = GraphUtils::possibleDsep(x, y, graph, maxPathLength);
    return dsep;
}

std::vector<Variable*> SepsetsPossibleDsep::possibleParents(Variable* x, std::vector<Variable*> nodes
                                                                /*IKnowledge knowledge*/) {
    std::vector<Variable*> possibleParents;
    std::string x_ = x->getName();

    for (Variable* z : nodes) {
        std::string z_ = z->getName();

        if (possibleParentOf(z_, x_ /*, knowledge*/)) {
            possibleParents.push_back(z);
        }
    }

    return possibleParents;
}
