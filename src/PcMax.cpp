#include "PcMax.hpp"

#include "GraphUtils.hpp"

/**
 * Constructs a new PC search using the given independence test as oracle.
 *
 * @param independenceTest The oracle for conditional independence facts. This does not make a copy of the
 *                         independence test, for fear of duplicating the data set!
 */
PcMax::PcMax(IndependenceTest *independenceTest) {
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
void PcMax::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

std::unordered_set<Edge> PcMax::getAdjacencies() {
    return graph.getEdges();
}

std::unordered_set<Edge> PcMax::getNonadjacencies() {
    EdgeListGraph complete = GraphUtils::completeGraph(graph);
    std::unordered_set<Edge> nonAdjacencies = complete.getEdges();
    EdgeListGraph undirected = GraphUtils::undirectedGraph(graph);
    for (Edge edge : undirected.getEdges()) {
        nonAdjacencies.erase(edge);
    }
    return nonAdjacencies;
}

/**
 * Runs PC search, returning the output pattern.
 */
EdgeListGraph PcMax::search() {
    return search(independenceTest->getVariables());
}

/**
 * Runs PC search, returning the output pattern, over the given nodes.
 */
EdgeListGraph PcMax::search(const std::vector<Variable*>& nodes) {
    Rcpp::Rcout << "Starting PCMax algorithm" << std::endl;

    if (independenceTest == NULL)
        throw std::invalid_argument("independenceTest of PCMax may not be NULL.");

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<Variable*> allNodes = independenceTest->getVariables();

    for (Variable* node : nodes) {
        if (std::find(allNodes.begin(), allNodes.end(), node) == allNodes.end())
            throw std::invalid_argument("All of the given nodes must be in the domain of the independence test provided.");
    }

    FasStableProducerConsumer fas(initialGraph, independenceTest);
    fas.setDepth(depth);
    graph = fas.search();

    OrientCollidersMaxP orientCollidersMaxP(independenceTest, &graph);

    orientCollidersMaxP.setUseHeuristic(useHeuristic);
    orientCollidersMaxP.setMaxPathLength(maxPathLength);
    orientCollidersMaxP.orient();

    Rcpp::Rcout << "Graph before Meek Rules: " << graph << std::endl;

    MeekRules meekRules;
    meekRules.orientImplied(graph);

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    Rcpp::Rcout << "PCMax Elapsed time =  " << elapsedTime << " ms" << std::endl;
    Rcpp::Rcout << "Finishing PCM Algorithm" << std::endl;

    return graph;
}

