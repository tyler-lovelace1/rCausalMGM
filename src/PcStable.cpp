#include "PcStable.hpp"

/**
 * Constructs a new PC search using the given independence test as oracle.
 *
 * @param independenceTest The oracle for conditional independence facts. This does not make a copy of the
 *                         independence test, for fear of duplicating the data set!
 */
PcStable::PcStable(IndependenceTest *independenceTest) {
    if (independenceTest == NULL) 
        throw std::invalid_argument("independenceTest may not be NULL.");

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
void PcStable::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

/**
 * Runs PC starting with a complete graph over all nodes of the given conditional independence test, using the given
 * independence test and knowledge and returns the resultant graph. The returned graph will be a pattern if the
 * independence information is consistent with the hypothesis that there are no latent common causes. It may,
 * however, contain cycles or bidirected edges if this assumption is not born out, either due to the actual presence
 * of latent common causes, or due to statistical errors in conditional independence judgments.
 */
EdgeListGraph PcStable::search() {
    return search(independenceTest->getVariables());
}

/**
 * Runs PC starting with a commplete graph over the given list of nodes, using the given independence test and
 * knowledge and returns the resultant graph. The returned graph will be a pattern if the independence information
 * is consistent with the hypothesis that there are no latent common causes. It may, however, contain cycles or
 * bidirected edges if this assumption is not born out, either due to the actual presence of latent common causes,
 * or due to statistical errors in conditional independence judgments.
 * <p>
 * All of the given nodes must be in the domain of the given conditional independence test.
 */
EdgeListGraph PcStable::search(const std::vector<Variable*>& nodes) {
    if (verbose) Rcpp::Rcout << "Starting PC algorithm" << std::endl;

    if (independenceTest == NULL)
        throw std::invalid_argument("independenceTest of PcStable may not be NULL.");

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<Variable*> allNodes = independenceTest->getVariables();

    for (Variable* node : nodes) {
        if (std::find(allNodes.begin(), allNodes.end(), node) == allNodes.end())
            throw std::invalid_argument("All of the given nodes must be in the domain of the independence test provided.");
    }
    
    // graph = EdgeListGraph(nodes);

    FasStableProducerConsumer fas(initialGraph, independenceTest);
    fas.setDepth(depth);
    fas.setVerbose(verbose);

    graph = fas.search();
    sepsets = fas.getSepsets();

    SearchGraphUtils::orientCollidersUsingSepsets(sepsets, graph, verbose);

    MeekRules rules;
    rules.setAggressivelyPreventCycles(aggressivelyPreventCycles);
    rules.orientImplied(graph);

    // Set algorithm and type
    std::ostringstream alg;
    alg << "PcStable: alpha = " << independenceTest->getAlpha();
    graph.setAlgorithm(alg.str());
    graph.setGraphType("markov equivalence class");

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    // Rcpp::Rcout << "Returning this graph: " << graph << std::endl;
    if (verbose) {
        Rcpp::Rcout.precision(2);
        Rcpp::Rcout << "PcStable Elapsed time =  " << (elapsedTime / 1000.0) << " s" << std::endl;
    }

    return graph;
}
