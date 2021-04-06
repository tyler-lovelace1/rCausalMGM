#include "SepsetsSet.hpp"

SepsetsSet::SepsetsSet(SepsetMap sepsets, IndependenceTest *test) {
    this->sepsets = sepsets;
    this->test = test;
}

bool SepsetsSet::isCollider(Variable* i, Variable* j, Variable* k) {
    std::vector<Variable*>* sepset = sepsets.get(i, k);
    // isIndependent(i, k, sepsets.get(i, k));
    return sepset != NULL && !(std::count(sepset->begin(), sepset->end(), j) >= 1);
}

bool SepsetsSet::isNoncollider(Variable* i, Variable* j, Variable* k) {
    std::vector<Variable*>* sepset = sepsets.get(i, k);
    // isIndependent(i, k, sepsets.get(i, k));
    return sepset != NULL && (std::count(sepset->begin(), sepset->end(), j) >= 1);
}
