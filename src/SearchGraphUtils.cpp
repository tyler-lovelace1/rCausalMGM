#include "SearchGraphUtils.hpp"


/**
 * Step C of PC; orients colliders using specified sepset. That is, orients x *-* y *-* z as x *-> y <-* z just in
 * case y is in Sepset({x, z}).
 */
std::vector<Triple> SearchGraphUtils::orientCollidersUsingSepsets(SepsetMap& set, EdgeListGraph& graph, bool verbose) {

    if (verbose) Rcpp::Rcout << "Starting Collider Orientation" << std::endl;

    std::vector<Triple> colliders;

    std::vector<Variable*> nodes = graph.getNodes();

    for (Variable* b : nodes) {
        std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) continue;

        ChoiceGenerator cg(adjacentNodes.size(), 2);

        std::vector<int> *combination;

        for (combination = cg.next(); combination != NULL; combination = cg.next()) {
            Variable* a = adjacentNodes[(*combination)[0]];
            Variable* c = adjacentNodes[(*combination)[1]];

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c)) continue;

            std::vector<Variable*>* sepset = set.get(a, c);

            if (sepset != NULL) {
                if (std::find(sepset->begin(), sepset->end(), b) == sepset->end()) { // !sepset.contains(b)
                    graph.setEndpoint(a, b, ENDPOINT_ARROW);
                    graph.setEndpoint(c, b, ENDPOINT_ARROW);
                    colliders.push_back(Triple(a, b, c));
                }
            }
        }
    }

    if (verbose) Rcpp::Rcout << "Finishing Collider Orientation" << std::endl;

    return colliders;
}