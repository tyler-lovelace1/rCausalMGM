#include "Fci.hpp"

//============================CONSTRUCTORS============================//

/**
 * Constructs a new FCI search for the given independence test and background knowledge.
 */
Fci::Fci(IndependenceTest* test) {
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
Fci::Fci(IndependenceTest* test, std::vector<Node> searchVars) {
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
    for (const Node& var: remVars) {
	std::remove(this->variables.begin(), this->variables.end(), var);
    }
}

//========================PUBLIC METHODS==========================//

EdgeListGraph Fci::search() {
    return search(test->getVariables());
}

EdgeListGraph Fci::search(const std::vector<Node>& nodes) {
    FasStableProducerConsumer fas(initialGraph, test, threads);

    fas.setKnowledge(knowledge);
    fas.setDepth(depth);
    fas.setVerbose(verbose);
    fas.setFDR(fdr);

    return search(fas, nodes);
}

EdgeListGraph Fci::search(FasStableProducerConsumer& fas, const std::vector<Node>& nodes) {
    if (verbose) Rcpp::Rcout << "Starting FCI-Stable algorithm..." << std::endl;

    whyOrient = std::unordered_map<std::string, std::string>();

    if (test == NULL)
        throw std::invalid_argument("independenceTest of FCI Stable may not be NULL.");

    auto startTime = std::chrono::high_resolution_clock::now();
    
    graph = fas.search();
    sepsets = fas.getSepsets();
    graph.reorientAllWith(ENDPOINT_CIRCLE);
    fasGraph = graph;
    fasSepsets = sepsets;
    
    // SepsetProducer spTemp(graph, sepsets, test);
    // SepsetMap possDsepSepsets;

    // The original FCI, with or without JiJi Zhang's orientation rules
    // Optional step: Possible Dsep. (Needed for correctness but very time consuming.)
    if (isPossibleDsepSearchDone()) {
	if (verbose) Rcpp::Rcout << "  Starting Posssible DSep search" << std::endl;

	FciOrient fciorient_;
	posDmapSp = SepsetProducer(graph, sepsets, test);
	posDsp = SepsetProducer(graph, test, nullSepsets, threads);

	if (verbose) Rcpp::Rcout << "    Starting Orientations..." << std::endl;
    
	if (orientRule == ORIENT_SEPSETS) {

	  posDmapSp.setOrientRule(orientRule);
	  posDmapSp.setDepth(depth);
	  posDmapSp.setVerbose(verbose);
	  posDmapSp.setKnowledge(knowledge);
	  if (verbose) Rcpp::Rcout << "      Filling Triple Map..." << std::endl;
      
	  posDmapSp.fillMap();

	  fciorient_ = FciOrient(posDmapSp, whyOrient);
      
	} else {

	  posDsp.setOrientRule(orientRule);
	  posDsp.setDepth(depth);
	  posDsp.setVerbose(verbose);
	  posDsp.setKnowledge(knowledge);
	  if (verbose) Rcpp::Rcout << "      Filling Triple Map..." << std::endl;
      
	  posDsp.fillMap();

	  fciorient_ = FciOrient(posDsp, whyOrient);
	}

	fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
	fciorient_.setMaxPathLength(maxPathLength);
	fciorient_.setKnowledge(knowledge);

	if (verbose) Rcpp::Rcout << "    Orienting colliders..." << std::endl;
    
	fciorient_.ruleR0(graph);

	PossibleDsepFciConsumerProducer possibleDSep(graph, test, threads);
	possibleDSep.setKnowledge(knowledge);
	possibleDSep.setDepth(getDepth());
	possibleDSep.setMaxPathLength(maxPathLength);
	possibleDSep.setVerbose(verbose);

	if (verbose) Rcpp::Rcout << "    Checking Possible Dsep sets..." << std::endl;

	possDsepSepsets = possibleDSep.search();
	sepsets.addAll(possDsepSepsets);

	graph = possibleDSep.getGraph();

	// Step FCI D.

	// Reorient all edges as o-o.
	graph.reorientAllWith(ENDPOINT_CIRCLE);
    }

    // Step CI C (Zhang's step F3.)
    //fciOrientbk(getKnowledge(), graph, independenceTest.getVariables());    - Robert Tillman 2008
    //        fciOrientbk(getKnowledge(), graph, variables);
    //        new FciOrient(graph, new Sepsets(this.sepsets)).ruleR0(new Sepsets(this.sepsets));
    if (verbose) Rcpp::Rcout << "  Starting Orientations..." << std::endl;

    FciOrient fciorient_;
    mapSp = SepsetProducer(graph, sepsets, test);
    sp = SepsetProducer(graph, test, possDsepSepsets, threads);
    
    if (orientRule == ORIENT_SEPSETS) {

	mapSp.setOrientRule(orientRule);
	mapSp.setDepth(depth);
	mapSp.setVerbose(verbose);
	mapSp.setKnowledge(knowledge);
	if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
	mapSp.fillMap();

	fciorient_ = FciOrient(mapSp, whyOrient);
      
    } else {

	sp.setOrientRule(orientRule);
	sp.setDepth(depth);
	sp.setVerbose(verbose);
	sp.setKnowledge(knowledge);
	if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
	sp.fillMap();

	fciorient_ = FciOrient(sp, whyOrient);
    }

    fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
    fciorient_.setMaxPathLength(maxPathLength);
    fciorient_.setKnowledge(knowledge);

    if (verbose) Rcpp::Rcout << "  Orienting colliders..." << std::endl;
    
    fciorient_.ruleR0(graph);

    if (orientRule == ORIENT_MAJORITY || orientRule == ORIENT_CONSERVATIVE) {
	for (auto t : sp.getAmbiguousTriples())
	    graph.addAmbiguousTriple(t.x, t.y, t.z);
    }

    if (verbose) Rcpp::Rcout << "  Orienting implied edges..." << std::endl;
    
    fciorient_.doFinalOrientation(graph);

    // Rcpp::Rcout << "Checking for almost cycles\n";
    // for (Edge e : graph.getEdges()) {
    // 	if(e.isBidirected() || e.isPartiallyOriented()) {
    // 	    Rcpp::Rcout << "  Checking edge " << e << "...";
    // 	    Node n1 = e.getNode1();
    // 	    Node n2 = e.getNode2();
    // 	    bool failFlag = false;
    // 	    if (graph.isAncestorOf(n1, n2)) {
    // 		Rcpp::Rcout << " Node " << n1 << " is an ancestor of " << n2 << "\n";
    // 		failFlag = true;
    // 	    }
    // 	    if (graph.isAncestorOf(n2, n1)) {
    // 		Rcpp::Rcout << " Node " << n2 << " is an ancestor of " << n1 << "\n";
    // 		failFlag = true;
    // 	    }
    // 	    if (!failFlag) {
    // 		Rcpp::Rcout << " Passed\n";
    // 	    }
    // 	}
    // }

    // // Set algorithm and type
    // std::ostringstream alg;
    // alg << "FCI Stable: alpha = " << test->getAlpha();
    // graph.setAlgorithm(alg.str());
    // graph.setGraphType("partial ancestral graph");

    // Set algorithm and type
    std::string algString;
    if (orientRule == ORIENT_SEPSETS)
	algString = "FCI-Stable";
    else if (orientRule == ORIENT_MAXP)
	algString = "FCI-Max";
    else if (orientRule == ORIENT_MAJORITY)
	algString = "MFCI-Stable";
    else if (orientRule == ORIENT_CONSERVATIVE)
	algString = "CFCI-Stable";
    else
	algString = "FCI-Stable";
    
    std::ostringstream alg;
    if (initialGraph==NULL) {
	alg << algString;
    } else {
	alg << initialGraph->getAlgorithm() << "-" << algString;
    }
    graph.setAlgorithm(alg.str());
    graph.setGraphType("partial ancestral graph");
    

    // auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    // if (verbose) {
    // 	if (elapsedTime < 100*1000) {
    // 	    Rcpp::Rcout.precision(2);
    // 	} else {
    // 	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
    // 	}
    //     Rcpp::Rcout << "FCI-Stable Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    // }

    double elapsedTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();

    if (verbose) {
        double factor = (elapsedTime < 10) ? std::pow(10, 2 - std::ceil(std::log10(std::abs(elapsedTime)))) : 1.0;
	elapsedTime = std::round(elapsedTime * factor) / factor;
        Rcpp::Rcout << "  FCI-Stable Elapsed Time =  " << elapsedTime << " s" << std::endl;
    }
    
    return graph;
}

EdgeListGraph Fci::reorientWithRule(OrientRule rule) {

    auto startTime = std::chrono::high_resolution_clock::now();

    if (verbose) Rcpp::Rcout << "Starting New Orientations..." << std::endl;

    whyOrient = std::unordered_map<std::string, std::string>();
    
    orientRule = rule;

    if (isPossibleDsepSearchDone()) {

	graph = fasGraph;
	sepsets = fasSepsets;
	graph.reorientAllWith(ENDPOINT_CIRCLE);    
	
	if (verbose) Rcpp::Rcout << "  Starting Posssible DSep search" << std::endl;

	FciOrient fciorient_;
	// posDmapSp = SepsetProducer(graph, sepsets, test);
	// posDsp = SepsetProducer(graph, test, nullSepsets, threads);

	if (verbose) Rcpp::Rcout << "    Starting Orientations..." << std::endl;
    
	if (orientRule == ORIENT_SEPSETS) {

	    posDmapSp.setOrientRule(orientRule);
	    posDmapSp.setDepth(depth);
	    posDmapSp.setVerbose(verbose);
	    posDmapSp.setKnowledge(knowledge);
	    if (verbose) Rcpp::Rcout << "      Filling Triple Map..." << std::endl;
      
	    posDmapSp.fillMap();

	    fciorient_ = FciOrient(posDmapSp, whyOrient);
      
	} else {

	    posDsp.setOrientRule(orientRule);
	    posDsp.setDepth(depth);
	    posDsp.setVerbose(verbose);
	    posDsp.setKnowledge(knowledge);
	    if (verbose) Rcpp::Rcout << "      Filling Triple Map..." << std::endl;
      
	    posDsp.fillMap();

	    fciorient_ = FciOrient(posDsp, whyOrient);
	}

	fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
	fciorient_.setMaxPathLength(maxPathLength);
	fciorient_.setKnowledge(knowledge);

	if (verbose) Rcpp::Rcout << "    Orienting colliders..." << std::endl;
    
	fciorient_.ruleR0(graph);

	PossibleDsepFciConsumerProducer possibleDSep(graph, test, threads);
	possibleDSep.setKnowledge(knowledge);
	possibleDSep.setDepth(getDepth());
	possibleDSep.setMaxPathLength(maxPathLength);
	possibleDSep.setVerbose(verbose);

	if (verbose) Rcpp::Rcout << "    Checking Possible Dsep sets..." << std::endl;

	possDsepSepsets = possibleDSep.search();
	sepsets.addAll(possDsepSepsets);

	graph = possibleDSep.getGraph();

	// Step FCI D.

	// Reorient all edges as o-o.
	graph.reorientAllWith(ENDPOINT_CIRCLE);

	if (verbose) Rcpp::Rcout << "  Starting Orientations..." << std::endl;

	mapSp = SepsetProducer(graph, sepsets, test);
	sp = SepsetProducer(graph, test, possDsepSepsets, threads);
    
	if (orientRule == ORIENT_SEPSETS) {

	    mapSp.setOrientRule(orientRule);
	    mapSp.setDepth(depth);
	    mapSp.setVerbose(verbose);
	    mapSp.setKnowledge(knowledge);
	    if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
	    mapSp.fillMap();

	    fciorient_ = FciOrient(mapSp, whyOrient);
      
	} else {

	    sp.setOrientRule(orientRule);
	    sp.setDepth(depth);
	    sp.setVerbose(verbose);
	    sp.setKnowledge(knowledge);
	    if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
	    sp.fillMap();

	    fciorient_ = FciOrient(sp, whyOrient);
	}

	fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
	fciorient_.setMaxPathLength(maxPathLength);
	fciorient_.setKnowledge(knowledge);

	if (verbose) Rcpp::Rcout << "  Orienting colliders..." << std::endl;
    
	fciorient_.ruleR0(graph);

	if (orientRule == ORIENT_MAJORITY || orientRule == ORIENT_CONSERVATIVE) {
	    for (auto t : sp.getAmbiguousTriples())
		graph.addAmbiguousTriple(t.x, t.y, t.z);
	}

	if (verbose) Rcpp::Rcout << "  Orienting implied edges..." << std::endl;
    
	fciorient_.doFinalOrientation(graph);
	
    } else {

	graph = fasGraph;
	sepsets = fasSepsets;
	
	graph.reorientAllWith(ENDPOINT_CIRCLE);
	graph.clearAmbiguousTriples();
    
	FciOrient fciorient_;
    
	if (orientRule == ORIENT_SEPSETS) {
	    mapSp.setOrientRule(orientRule);
	    mapSp.setDepth(depth);
	    mapSp.setVerbose(verbose);
	    mapSp.setKnowledge(knowledge);
	    if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
	    mapSp.fillMap();

	    fciorient_ = FciOrient(mapSp, whyOrient);
      
	} else {
	    sp.setOrientRule(orientRule);
	    sp.setDepth(depth);
	    sp.setVerbose(verbose);
	    sp.setKnowledge(knowledge);
	    if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
	    sp.fillMap();

	    fciorient_ = FciOrient(sp, whyOrient);
	}
    
	// sp.setOrientRule(orientRule);
	// sp.setDepth(depth);
	// sp.setVerbose(verbose);
	// sp.setKnowledge(knowledge);
	// if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
    
	// sp.fillMap();
    
	// FciOrient fciorient_(sp, whyOrient);
	fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
	fciorient_.setMaxPathLength(maxPathLength);
	fciorient_.setKnowledge(knowledge);
	fciorient_.ruleR0(graph);

	if (orientRule == ORIENT_MAJORITY || orientRule == ORIENT_CONSERVATIVE) {
	    for (auto t : sp.getAmbiguousTriples())
		graph.addAmbiguousTriple(t.x, t.y, t.z);
	}
    
	fciorient_.doFinalOrientation(graph);

    }

    // // Set algorithm and type
    // std::ostringstream alg;
    // alg << "FCI Stable: alpha = " << test->getAlpha();
    // graph.setAlgorithm(alg.str());
    // graph.setGraphType("partial ancestral graph");

    // Set algorithm and type
    std::string algString;
    if (orientRule == ORIENT_SEPSETS)
	algString = "FCI-Stable";
    else if (orientRule == ORIENT_MAXP)
	algString = "FCI-Max";
    else if (orientRule == ORIENT_MAJORITY)
	algString = "MFCI-Stable";
    else if (orientRule == ORIENT_CONSERVATIVE)
	algString = "CFCI-Stable";
    else
	algString = "FCI-Stable";
    
    std::ostringstream alg;
    if (initialGraph==NULL) {
	alg << algString;
    } else {
	alg << initialGraph->getAlgorithm() << "-" << algString;
    }
    graph.setAlgorithm(alg.str());
    graph.setGraphType("partial ancestral graph");
    

    // auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    // if (verbose) {
    // 	if (elapsedTime < 100*1000) {
    // 	    Rcpp::Rcout.precision(2);
    // 	} else {
    // 	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
    // 	}
    //     Rcpp::Rcout << "Reorient Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    // }

    double elapsedTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - startTime).count();

    if (verbose) {
        double factor = (elapsedTime < 10) ? std::pow(10, 2 - std::ceil(std::log10(std::abs(elapsedTime)))) : 1.0;
	elapsedTime = std::round(elapsedTime * factor) / factor;
        Rcpp::Rcout << "  Reorient Elapsed Time =  " << elapsedTime << " s" << std::endl;
    }
    
    return graph;
}

/**
 * @param maxPathLength the maximum length of any discriminating path, or -1 if unlimited.
 */
void Fci::setMaxPathLength(int maxPathLength) {
    if (maxPathLength < -1) {
        throw std::invalid_argument("Max path length must be -1 (unlimited) or >= 0: ");
    }
    this->maxPathLength = maxPathLength;
}

//========================PRIVATE METHODS==========================//

void Fci::buildIndexing(std::vector<Node> nodes) {
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
void Fci::fciOrientbk(/*IKnowledge bk,*/ EdgeListGraph graph, std::vector<Node> variables) {

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
