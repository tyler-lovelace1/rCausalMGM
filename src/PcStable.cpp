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


EdgeListGraph PcStable::search(const std::vector<Node>& nodes) {
    FasStableProducerConsumer fas(initialGraph, independenceTest, threads);

    fas.setDepth(depth);
    fas.setVerbose(verbose);
    fas.setFDR(fdr);
    fas.setKnowledge(knowledge);

    return search(fas, nodes);
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
EdgeListGraph PcStable::search(FasStableProducerConsumer& fas, const std::vector<Node>& nodes) {
    if (verbose) Rcpp::Rcout << "Starting PC-Stable algorithm..." << std::endl;

    if (independenceTest == NULL)
        throw std::invalid_argument("independenceTest of PcStable may not be NULL.");

    auto startTime = std::chrono::high_resolution_clock::now();

    std::vector<Node> allNodes = independenceTest->getVariables();

    for (const Node& node : nodes) {
        if (std::find(allNodes.begin(), allNodes.end(), node) == allNodes.end())
            throw std::invalid_argument("All of the given nodes must be in the domain of the independence test provided.");
    }
    
    // graph = EdgeListGraph(nodes);

    // FasStableProducerConsumer fas(initialGraph, independenceTest, threads);

    graph = fas.search();
    sepsets = fas.getSepsets();

    SepsetProducer sp;
    SepsetMap nullSepsets;
    
    // if (orientRule == ORIENT_SEPSETS) {
    // 	Rcpp::Rcout << "Sepset Map Orientations\n";
    // 	// sp = SepsetProducer(sepsets, independenceTest);
    // 	sp = SepsetProducer(graph, independenceTest, sepsets, threads);
    // } else {
    sp = SepsetProducer(graph, independenceTest, sepsets, threads);
    // }

    sp.setOrientRule(orientRule);
    sp.setDepth(depth);
    sp.setVerbose(verbose);
    sp.setKnowledge(knowledge);
    
    if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;

    sp.fillMap();
    
    if (verbose) Rcpp::Rcout << "  Orienting edges with knowledge..." << std::endl;

    SearchGraphUtils::pcOrientbk(knowledge, graph);
    // SearchGraphUtils::orientCollidersUsingSepsets(sepsets, graph, knowledge, verbose);

    MeekRules rules;
    rules.setAggressivelyPreventCycles(aggressivelyPreventCycles);
    rules.setKnowledge(knowledge);

    if (verbose) Rcpp::Rcout << "  Orienting colliders..." << std::endl;
    // rules.meekR0(graph);
    
    std::vector<Triple> orderedColliders = sp.getOrderedColliders();

    // if (false) {
    // 	SearchGraphUtils::orientCollidersUsingSepsets(sepsets, graph, knowledge, verbose);
    // } else {
    SearchGraphUtils::orientCollidersUsingOrderedColliders(orderedColliders, graph,
							   knowledge, verbose);
    // }
    
    for (auto t : sp.getAmbiguousTriples())
	graph.addAmbiguousTriple(t.x, t.y, t.z);

    if (verbose) Rcpp::Rcout << "  Orienting implied edges..." << std::endl;
    
    rules.orientImplied(graph);

    std::set<Edge> edgeSet = graph.getEdges();

    for (Edge edge : edgeSet) {
        if (edge.isBidirected()) {
	    Node node1 = edge.getNode1();
	    Node node2 = edge.getNode2();

	    graph.removeEdge(node1, node2);
	    graph.addUndirectedEdge(node1, node2);
	}
    }

    score = 0.0;
    for (int i = 0; i < allNodes.size(); i++) {
	for (int j = 0; j < i; j++) {
	    if (sepsets.isInSepsetMap(allNodes[i], allNodes[j])) {
		score += std::log(sepsets.getPValue(allNodes[i], allNodes[j]));
	    }	    
	}
    }

    // Set algorithm and type
    std::string algString;
    if (orientRule == ORIENT_SEPSETS)
	algString = "PC-Stable";
    else if (orientRule == ORIENT_MAXP)
	algString = "PC-Max";
    else if (orientRule == ORIENT_MAJORITY)
	algString = "MPC-Stable";
    else if (orientRule == ORIENT_CONSERVATIVE)
	algString = "CPC-Stable";
    else
	algString = "PC-Stable";
    
    std::ostringstream alg;
    if (initialGraph==NULL) {
	alg << algString;
    } else {
	alg << initialGraph->getAlgorithm() << "-" << algString;
    }
    graph.setAlgorithm(alg.str());
    graph.setGraphType("completed partially directed acyclic graph");

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    if (verbose) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "PC-Stable Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

    return graph;
}
