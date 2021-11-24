#ifndef SEARCHGRAPHUTILS_HPP_
#define SEARCHGRAPHUTILS_HPP_

#include "EdgeListGraph.hpp"
#include "SepsetMap.hpp"
#include "ChoiceGenerator.hpp"

class SearchGraphUtils {

public:
    /**
     * Step C of PC; orients colliders using specified sepset. That is, orients x *-* y *-* z as x *-> y <-* z just in
     * case y is in Sepset({x, z}).
     */
    static std::vector<Triple> orientCollidersUsingSepsets(SepsetMap& set, EdgeListGraph& graph, bool verbose = false);

    static bool orientCollider(SepsetMap& set, EdgeListGraph& graph, Variable* a, Variable* b, Variable* c);

    static bool wouldCreateBadCollider(SepsetMap& set, EdgeListGraph& graph, Variable* x, Variable* y);

};

#endif /* SEARCHGRAPHUTILS_HPP_ */
