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
    static std::vector<Triple> orientCollidersUsingSepsets(SepsetMap& set, EdgeListGraph& graph);

};

#endif /* SEARCHGRAPHUTILS_HPP_ */