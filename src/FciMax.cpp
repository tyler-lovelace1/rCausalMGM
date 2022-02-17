#include "FciMax.hpp"

//============================CONSTRUCTORS============================//

/**
 * Constructs a new FCI search for the given independence test and background knowledge.
 */
FciMax::FciMax(IndependenceTest* test) {
    if (test == NULL /*|| knowledge == null*/) {
        throw std::invalid_argument("independenceTest cannot be null");
    }

    this->test = test;
    std::vector<Node> vars = test->getVariables();
    this->variables.insert(variables.end(), vars.begin(), vars.end());
    buildIndexing(vars);
}

/**
 * Constructs a new FCI search for the given independence test and background knowledge and a list of variables to
 * search over.
 */
FciMax::FciMax(IndependenceTest* test, std::vector<Node> searchVars) {
    if (test == NULL /*|| knowledge == null*/) {
        throw std::invalid_argument("independenceTest cannot be null");
    }

    this->test = test;
    std::vector<Node> vars = test->getVariables();
    this->variables.insert(variables.end(), vars.begin(), vars.end());

    std::unordered_set<Node> remVars;
    for (const Node& node1 : this->variables) {
        bool search = false;
        for (const Node& node2 : searchVars) {
            if (node1 == node2) {
                search = true;
            }
        }
        if (!search) {
            remVars.insert(node1);
        }
    }
    for (const Node& var : remVars) {
	std::remove(this->variables.begin(), this->variables.end(), var);
    }
}

//========================PUBLIC METHODS==========================//

EdgeListGraph FciMax::search() {
    return search(test->getVariables());
}

EdgeListGraph FciMax::search(const std::vector<Node>& nodes) {
    FasStableProducerConsumer fas(initialGraph, test, threads);

    return search(fas, nodes);
}

EdgeListGraph FciMax::search(FasStableProducerConsumer& fas, const std::vector<Node>& nodes) {

    auto startTime = std::chrono::high_resolution_clock::now();
    
    if (verbose) Rcpp::Rcout << "Starting FCI-Max algorithm..." << std::endl;

    whyOrient = std::unordered_map<std::string, std::string>();

    if (test == NULL)
        throw std::invalid_argument("independenceTest of FCI-Max may not be NULL.");

    // fas->setKnowledge(getKnowledge());
    fas.setDepth(depth);
    fas.setVerbose(verbose);
    fas.setFDR(fdr);
    graph = fas.search();
    sepsets = fas.getSepsets();
    graph.reorientAllWith(ENDPOINT_CIRCLE);

    // SepsetsPossibleDsep sp(graph, test,/* knowledge,*/ depth, maxPathLength);

    // The original FCI, with or without JiJi Zhang's orientation rules
    // Optional step: Possible Dsep. (Needed for correctness but very time consuming.)
    if (isPossibleDsepSearchDone()) {
	if (verbose) Rcpp::Rcout << "  Starting Posssible DSep search..." << std::endl;
	// SepsetsSet ssset(sepsets, test);
	// FciMaxOrient orienter(&ssset);

	PossibleDsepFciConsumerProducer possibleDSep(graph, test, threads);
	// possibleDSep.setKnowledge(getKnowledge());
	possibleDSep.setDepth(getDepth());
	possibleDSep.setMaxPathLength(maxPathLength);
	possibleDSep.setVerbose(verbose);

	SepsetMap searchResult = possibleDSep.search();
	sepsets.addAll(searchResult);

	graph = possibleDSep.getGraph();

	// Step FCI D.

	// Reorient all edges as o-o.
	graph.reorientAllWith(ENDPOINT_CIRCLE);
    }

    if (verbose) Rcpp::Rcout << "  Starting Orientations..." << std::endl;

    SepsetProducerMaxP sepsetsMaxP(graph, test, sepsets, threads);
    sepsetsMaxP.setDepth(depth);
    sepsetsMaxP.setVerbose(verbose);
    sepsetsMaxP.fillMap();

    // Step CI C (Zhang's step F3.)
    //fciOrientbk(getKnowledge(), graph, independenceTest.getVariables());    - Robert Tillman 2008
    //        fciOrientbk(getKnowledge(), graph, variables);
    //        new FciMaxOrient(graph, new Sepsets(this.sepsets)).ruleR0(new Sepsets(this.sepsets));
    //SepsetsSet sepsetsset_(sepsets, test);
    FciOrient fciorient_(&sepsetsMaxP, whyOrient);
    fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
    fciorient_.setMaxPathLength(maxPathLength);
    // fciOrient.setKnowledge(knowledge);
    // if (verbose) Rcpp::Rcout << "FCIOrient Set Up" << std::endl;
    fciorient_.ruleR0(graph);
    // if (verbose) Rcpp::Rcout << "Rule 0 finished" << std::endl;
    fciorient_.doFinalOrientation(graph);

    // // Set algorithm and type
    // std::ostringstream alg;
    // alg << "FCI-Max: alpha = " << test->getAlpha();
    // graph.setAlgorithm(alg.str());
    // graph.setGraphType("partial ancestral graph");

    // Set algorithm and type
    std::ostringstream alg;
    if (initialGraph==NULL) {
	alg << "FCI-Max";
    } else {
	alg << initialGraph->getAlgorithm() << "-" << "FCI-Max";
    }
    graph.setAlgorithm(alg.str());
    graph.setGraphType("partial ancestral graph");
    
    
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    if (verbose) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "FCI-Max Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }
    
    return graph;
}

/**
 * @param maxPathLength the maximum length of any discriminating path, or -1 if unlimited.
 */
void FciMax::setMaxPathLength(int maxPathLength) {
    if (maxPathLength < -1) {
        throw std::invalid_argument("Max path length must be -1 (unlimited) or >= 0: ");
    }
    this->maxPathLength = maxPathLength;
}

//========================PRIVATE METHODS==========================//

void FciMax::buildIndexing(std::vector<Node> nodes) {
    this->hashIndices =  std::unordered_map<Node, int>();
    for (const Node& node : nodes) {
        auto itr = find(variables.begin(), variables.end(), node);
        int index = std::distance(variables.begin(), itr);
        this->hashIndices.insert(std::pair<Node, int>(node, index));
    }
}

/**
 * Orients according to background knowledge
 */
void FciMax::fciOrientbk(/*IKnowledge bk,*/ EdgeListGraph graph, std::vector<Node> variables) {

    // for (Iterator<KnowledgeEdge> it =
    //      bk.forbiddenEdgesIterator(); it.hasNext(); ) {
    //     KnowledgeEdge edge = it.next();
    //
    //     //match strings to variables in the graph.
    //     Node from = SearchGraphUtils.translate(edge.getFrom(), variables);
    //     Node to = SearchGraphUtils.translate(edge.getTo(), variables);
    //
    //
    //     if (from == null || to == null) {
    //         continue;
    //     }
    //
    //     if (graph.getEdge(from, to) == null) {
    //         continue;
    //     }
    //
    //     // Orient to*->from
    //     graph.setEndpoint(to, from, Endpoint.ARROW);
    //     graph.setEndpoint(from, to, Endpoint.CIRCLE);
    //     logger.log("knowledgeOrientation", SearchLogUtils.edgeOrientedMsg("Knowledge", graph.getEdge(from, to)));
    // }
    //
    // for (Iterator<KnowledgeEdge> it =
    //      bk.requiredEdgesIterator(); it.hasNext(); ) {
    //     KnowledgeEdge edge = it.next();
    //
    //     //match strings to variables in this graph
    //     Node from = SearchGraphUtils.translate(edge.getFrom(), variables);
    //     Node to = SearchGraphUtils.translate(edge.getTo(), variables);
    //
    //     if (from == null || to == null) {
    //         continue;
    //     }
    //
    //     if (graph.getEdge(from, to) == null) {
    //         continue;
    //     }
    //
    //     graph.setEndpoint(to, from, Endpoint.TAIL);
    //     graph.setEndpoint(from, to, Endpoint.ARROW);
    // }
    return;
}
