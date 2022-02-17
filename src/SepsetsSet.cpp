#include "SepsetsSet.hpp"

SepsetsSet::SepsetsSet(SepsetMap sepsets, IndependenceTest *test) {
    this->sepsets = sepsets;
    this->test = test;
}

bool SepsetsSet::isCollider(const Node& i, const Node& j, const Node& k) {
    std::vector<Node>* sepset = sepsets.get(i, k);
    // isIndependent(i, k, sepsets.get(i, k));
    return sepset != NULL && !(std::count(sepset->begin(), sepset->end(), j) >= 1);
}

bool SepsetsSet::isNoncollider(const Node& i, const Node& j, const Node& k) {
    std::vector<Node>* sepset = sepsets.get(i, k);
    // isIndependent(i, k, sepsets.get(i, k));
    return sepset != NULL && (std::count(sepset->begin(), sepset->end(), j) >= 1);
}
