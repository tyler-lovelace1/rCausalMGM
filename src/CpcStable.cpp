#include "CpcStable.hpp"

bool CpcStable::isColliderSepset(Variable* j, std::vector<std::vector<Variable*>>& sepsets) {
    if (sepsets.size() == 0) return false;

    for (std::vector<Variable*> sepset : sepsets) {
        // if sepset.contains(j)
        if (std::find(sepset.begin(), sepset.end(), j) != sepset.end()) return false;
    }

    return true;
}

std::vector<std::vector<Variable*>> CpcStable::getSepsets(Variable* i, Variable* k, EdgeListGraph& g) {
    std::vector<Variable*> adji = g.getAdjacentNodes(i);
    std::vector<Variable*> adjk = g.getAdjacentNodes(k);
    std::vector<std::vector<Variable*>> sepsets;

    for (int d = 0; d <= std::max(adji.size(), adjk.size()); d++) {
        if (adji.size() >= 2 && d <= adji.size()) {
            ChoiceGenerator gen(adji.size(), d);

            std::vector<int> *choice;
            for (choice = gen.next(); choice != NULL; choice = gen.next()) {
                std::vector<Variable*> v = GraphUtils::asList(*choice, adji);
                if (independenceTest->isIndependent(i, k, v)) sepsets.push_back(v);
            }
        }

        if (adjk.size() >= 2 && d <= adjk.size()) {
            ChoiceGenerator gen(adjk.size(), d);

            std::vector<int> *choice;
            for (choice = gen.next(); choice != NULL; choice = gen.next()) {
                std::vector<Variable*> v = GraphUtils::asList(*choice, adjk);
                if (independenceTest->isIndependent(i, k, v)) sepsets.push_back(v);
            }

        }
    }

    return sepsets;
}

void CpcStable::orientUnshieldedTriples() {
    std::vector<Variable*> nodes = graph.getNodes();

    for (Variable* y : nodes) {
        std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(y);

        if (adjacentNodes.size() < 2)
            continue;

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *combination;
        for (combination = cg.next(); combination != NULL; combination = cg.next()) {
            Variable* x = adjacentNodes[(*combination)[0]];
            Variable* z = adjacentNodes[(*combination)[1]];

            if (graph.isAdjacentTo(x, z))
                continue;

            std::vector<std::vector<Variable*>> sepsetsxz = getSepsets(x, z, graph);

            if (isColliderSepset(y, sepsetsxz)) {
                // colliderAllowed (knowledge)
                if (true) {
                    graph.setEndpoint(x, y, ENDPOINT_ARROW);
                    graph.setEndpoint(z, y, ENDPOINT_ARROW);
                }
            } else {
                graph.addAmbiguousTriple(x, y, z);
            }

            allTriples.insert(Triple(x, y, z));
        }
    }
}

/**
 * Constructs a CPC algorithm that uses the given independence test as oracle. This does not make a copy of the
 * independence test, for fear of duplicating the data set!
 */
CpcStable::CpcStable(IndependenceTest *independenceTest) {
    if (independenceTest == NULL) {
        throw std::invalid_argument("independenceTest may not be NULL.");
    } 

    this->independenceTest = independenceTest;
}

/**
 * Sets the depth of the search--that is, the maximum number of conditioning nodes for any conditional independence
 * checked.
 *
 * @param depth The depth of the search. The default is 1000. A value of -1 may be used to indicate that the depth
 *              should be high (1000). A value of Integer.MAX_VALUE may not be used, due to a bug on multi-core
 *              machines.
 */
void CpcStable::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

std::unordered_set<Edge> CpcStable::getAdjacencies() {
    return graph.getEdges();
}

std::unordered_set<Edge> CpcStable::getNonadjacencies() {
    EdgeListGraph complete = GraphUtils::completeGraph(graph);
    std::unordered_set<Edge> nonAdjacencies = complete.getEdges();
    EdgeListGraph undirected = GraphUtils::undirectedGraph(graph);
    for (Edge edge : undirected.getEdges()) {
        nonAdjacencies.erase(edge);
    }
    return nonAdjacencies;
}

/**
 * Runs PC starting with a fully connected graph over all of the variables in the domain of the independence test.
 * See PC for caveats. The number of possible cycles and bidirected edges is far less with CPC than with PC.
 */
EdgeListGraph CpcStable::search() {
    return search(independenceTest->getVariables());
}

EdgeListGraph CpcStable::search(const std::vector<Variable*>& nodes) {
    FasStable fas(initialGraph, independenceTest);
    return search(fas, nodes);
}

EdgeListGraph CpcStable::search(FasStable& fas, const std::vector<Variable*>& nodes) {
    Rcpp::Rcout << "Starting CPC algorithm" << std::endl;

    allTriples = {};

    if (independenceTest == NULL)
        throw std::invalid_argument("independenceTest of CpcStable may not be NULL.");

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<Variable*> allNodes = independenceTest->getVariables();

    for (Variable* node : nodes) {
        if (std::find(allNodes.begin(), allNodes.end(), node) == allNodes.end())
            throw std::invalid_argument("All of the given nodes must be in the domain of the independence test provided.");
    }

    fas.setDepth(depth);

    // Note that we are ignoring the sepset map returned by this method
    // on purpose; it is not used in this search.
    graph = fas.search();
    sepsets = fas.getSepsets();

    orientUnshieldedTriples();

    MeekRules meekRules;
    meekRules.orientImplied(graph);

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    Rcpp::Rcout << "Returning this graph: " << graph << std::endl;
    Rcpp::Rcout << "CpcStable Elapsed time =  " << elapsedTime << " ms" << std::endl;
    Rcpp::Rcout << "Finishing CPC Algorithm" << std::endl;

    return graph;
}