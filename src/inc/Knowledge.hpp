#ifndef KNOWLEDGE_HPP_
#define KNOWLEDGE_HPP_

#include "armaLapack.hpp"

#include "Node.hpp"
#include <vector>
#include <set>

class Knowledge {

private:

    // A list of the nodes in the graph, in the order in which they
    // were added.
    std::vector<Node> nodes;

    // Tier based causal ordering
    std::vector<std::set<Node>> tiers;

    // Whether or not edges are allowed within tiers
    std::vector<bool> forbiddenWithinTier;

    // Forbidden directed pairs: forbids Node1 *-> Node2, but does not
    // forbid Node2 *-> Node1. To forbid all adjacencies, variables
    // must be added forbidding both directions
    std::set<std::pair<Node, Node>> forbidden;

    // Required directed pairs: requires Node1 *-> Node2
    std::set<std::pair<Node, Node>> required;

public:

    Knowledge() {}

    Knowledge(std::vector<Node> nodes, Rcpp::List& knowledge);

    bool isEmpty() { return tiers.empty() && forbidden.empty() && required.empty(); }

    bool isForbidden(const Node& node1, const Node& node2);

    bool isForbiddenByTiers(const Node& node1, const Node& node2);

    bool isRequired(const Node& node1, const Node& node2);

    bool noEdgeRequired(const Node& node1, const Node& node2);

    std::vector<std::set<Node>> getTiers() { return tiers; }

    std::set<std::pair<Node, Node>> getForbidden() { return forbidden; }

    std::set<std::pair<Node, Node>> getRequired() { return required; }

};

#endif /* KNOWLEDGE_HPP_ */
