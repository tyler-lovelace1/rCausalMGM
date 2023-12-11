#ifndef SEARCHGRAPHUTILS_HPP_
#define SEARCHGRAPHUTILS_HPP_

#include "EdgeListGraph.hpp"
#include "Knowledge.hpp"
#include "SepsetMap.hpp"
#include "SepsetProducer.hpp"
#include "ChoiceGenerator.hpp"

class SearchGraphUtils {

public:

    static void pcOrientbk(Knowledge& knowledge, EdgeListGraph& graph);

    static bool isArrowheadAllowed(const Node& from, const Node& to, Knowledge& knowledge);
    
    /**
     * Step C of PC; orients colliders using specified sepset. That is, orients x *-* y *-* z as x *-> y <-* z just in
     * case y is in Sepset({x, z}).
     */
    static std::vector<Triple> orientCollidersUsingSepsets(SepsetMap& set,
							   EdgeListGraph& graph,
							   Knowledge& knowledge,
							   bool verbose = false);

    static void orientCollidersUsingOrderedColliders(std::vector<Triple>& orderedColliders,
						     EdgeListGraph& graph,
						     Knowledge& knowledge,
						     bool verbose = false);

    static bool orientCollider(SepsetMap& set, EdgeListGraph& graph, const Node& a, const Node& b, const Node& c);

    static bool wouldCreateBadCollider(SepsetMap& set, EdgeListGraph& graph, const Node& x, const Node& y);

    static bool orientCollider(std::vector<Triple>& colliders, EdgeListGraph& graph, const Node& a, const Node& b, const Node& c);

    static bool wouldCreateBadCollider(std::vector<Triple>& colliders, EdgeListGraph& graph, const Node& x, const Node& y);

};

#endif /* SEARCHGRAPHUTILS_HPP_ */
