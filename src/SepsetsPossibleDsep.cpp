#include "SepsetsPossibleDsep.hpp"

SepsetsPossibleDsep::SepsetsPossibleDsep(EdgeListGraph& graph, IndependenceTest *test, /*IKnowledge knowledge,*/
					 int depth, int maxPathLength) : SepsetProducer(graph, test) {
    // this->graph = graph;
    // this->test = test;
    this->maxPathLength = maxPathLength;
    // this.knowledge = knowledge;
    this->depth = depth;
}

std::vector<Node> SepsetsPossibleDsep::getSepset(const Node& i, const Node& k) {
    std::vector<Node> condSet = getCondSet(i, k, maxPathLength);

    if (condSet.empty()) {
        condSet = getCondSet(k, i, maxPathLength);
    }

    return condSet;
}

bool SepsetsPossibleDsep::isCollider(const Node& i, const Node& j, const Node& k) {
    std::vector<Node> sepset = getSepset(i, k);
    return (!sepset.empty())  && !(std::count(sepset.begin(), sepset.end(), j) == 1);
}

bool SepsetsPossibleDsep::isNoncollider(const Node& i, const Node& j, const Node& k) {
    std::vector<Node> sepset = getSepset(i, k);
    return (!sepset.empty()) && (std::count(sepset.begin(), sepset.end(), j) == 1);
}

std::vector<Node> SepsetsPossibleDsep::getCondSet(const Node& node1, const Node& node2, int maxPathLength) {
    std::unordered_set<Node> possibleDsepSet = getPossibleDsep(node1, node2, maxPathLength);
    std::vector<Node> possibleDsep (possibleDsepSet.begin(), possibleDsepSet.end());

    // bool noEdgeRequired = knowledge.noEdgeRequired(node1.getName(), node2.getName());

    std::vector<Node> possParents = possibleParents(node1, possibleDsep /*, knowledge*/);

    size_t depth_ = depth == -1 ? 1000 : depth;

    for (int d = 0; d <= std::min(depth_, possParents.size()); d++) {


        ChoiceGenerator cg (possParents.size(), d);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            std::vector<Node> condSet = GraphUtils::asList(*choice, possParents);
            bool independent = test->isIndependent(node1, node2, condSet);

            if (independent /*&& noEdgeRequired*/) {
                return condSet;
            }
        }
    }
    std::vector<Node> nullSet;
    return nullSet;
}

std::unordered_set<Node> SepsetsPossibleDsep::getPossibleDsep(const Node& x, const Node& y, int maxPathLength) {
    std::unordered_set<Node> dsep = GraphUtils::possibleDsep(x, y, graph, maxPathLength);
    return dsep;
}

std::vector<Node> SepsetsPossibleDsep::possibleParents(const Node& x, std::vector<Node> nodes
                                                                /*IKnowledge knowledge*/) {
    std::vector<Node> possibleParents;
    std::string x_ = x.getName();

    for (const Node& z : nodes) {
        std::string z_ = z.getName();

        if (possibleParentOf(z_, x_ /*, knowledge*/)) {
            possibleParents.push_back(z);
        }
    }

    return possibleParents;
}
