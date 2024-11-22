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

    SepsetMap rfciSepsets = ruleR0_rfciPrune(graph, sepsets);

    // The original FCI, with or without JiJi Zhang's orientation rules
    // Optional step: Possible Dsep. (Needed for correctness but very time consuming.)
    if (isPossibleDsepSearchDone()) {
	if (verbose) Rcpp::Rcout << "  Starting Posssible DSep search" << std::endl;

	FciOrient fciorient_;
	posDsp = SepsetProducer(graph, test, nullSepsets, threads);

	if (verbose) Rcpp::Rcout << "    Starting Conservative Orientations..." << std::endl;
    
	posDsp.setOrientRule(ORIENT_CONSERVATIVE);
	posDsp.setDepth(depth);
	posDsp.setVerbose(verbose);
	posDsp.setKnowledge(knowledge);
	if (verbose) Rcpp::Rcout << "      Filling Triple Map..." << std::endl;
	
	posDsp.fillMap();
	
	fciorient_ = FciOrient(posDsp, whyOrient);
	
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

	if (verbose) Rcpp::Rcout << "    Checking Possible-Dsep sets..." << std::endl;

	possDsepSepsets = possibleDSep.search();
	possDsepSepsets.addAll(rfciSepsets);
	sepsets.addAll(possDsepSepsets);

	graph = possibleDSep.getGraph();

	// Step FCI D.

	// Reorient all edges as o-o.
	graph.reorientAllWith(ENDPOINT_CIRCLE);
    } else {
	possDsepSepsets = rfciSepsets;
	sepsets.addAll(possDsepSepsets);
    }

    // Step CI C (Zhang's step F3.)
    //fciOrientbk(getKnowledge(), graph, independenceTest.getVariables());    - Robert Tillman 2008
    //        fciOrientbk(getKnowledge(), graph, variables);
    //        new FciOrient(graph, new Sepsets(this.sepsets)).ruleR0(new Sepsets(this.sepsets));
    if (verbose) Rcpp::Rcout << "  Starting Final Orientations..." << std::endl;

    FciOrient fciorient_;
    mapSp = SepsetProducer(graph, sepsets, test);
    sp = SepsetProducer(graph, test, sepsets, threads);
    
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

    // if (isPossibleDsepSearchDone()) {

    // 	graph = fasGraph;
    // 	sepsets = fasSepsets;
    // 	graph.reorientAllWith(ENDPOINT_CIRCLE);    
	
    // 	if (verbose) Rcpp::Rcout << "  Starting Posssible DSep search" << std::endl;

    // 	FciOrient fciorient_;
    // 	// posDmapSp = SepsetProducer(graph, sepsets, test);
    // 	// posDsp = SepsetProducer(graph, test, nullSepsets, threads);

    // 	if (verbose) Rcpp::Rcout << "    Starting Orientations..." << std::endl;
    
    // 	if (orientRule == ORIENT_SEPSETS) {

    // 	    posDmapSp.setOrientRule(orientRule);
    // 	    posDmapSp.setDepth(depth);
    // 	    posDmapSp.setVerbose(verbose);
    // 	    posDmapSp.setKnowledge(knowledge);
    // 	    if (verbose) Rcpp::Rcout << "      Filling Triple Map..." << std::endl;
      
    // 	    posDmapSp.fillMap();

    // 	    fciorient_ = FciOrient(posDmapSp, whyOrient);
      
    // 	} else {

    // 	    posDsp.setOrientRule(orientRule);
    // 	    posDsp.setDepth(depth);
    // 	    posDsp.setVerbose(verbose);
    // 	    posDsp.setKnowledge(knowledge);
    // 	    if (verbose) Rcpp::Rcout << "      Filling Triple Map..." << std::endl;
      
    // 	    posDsp.fillMap();

    // 	    fciorient_ = FciOrient(posDsp, whyOrient);
    // 	}

    // 	fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
    // 	fciorient_.setMaxPathLength(maxPathLength);
    // 	fciorient_.setKnowledge(knowledge);

    // 	if (verbose) Rcpp::Rcout << "    Orienting colliders..." << std::endl;
    
    // 	fciorient_.ruleR0(graph);

    // 	PossibleDsepFciConsumerProducer possibleDSep(graph, test, threads);
    // 	possibleDSep.setKnowledge(knowledge);
    // 	possibleDSep.setDepth(getDepth());
    // 	possibleDSep.setMaxPathLength(maxPathLength);
    // 	possibleDSep.setVerbose(verbose);

    // 	if (verbose) Rcpp::Rcout << "    Checking Possible Dsep sets..." << std::endl;

    // 	possDsepSepsets = possibleDSep.search();
    // 	sepsets.addAll(possDsepSepsets);

    // 	graph = possibleDSep.getGraph();

    // 	// Step FCI D.

    // 	// Reorient all edges as o-o.
    // 	graph.reorientAllWith(ENDPOINT_CIRCLE);

    // 	if (verbose) Rcpp::Rcout << "  Starting Orientations..." << std::endl;

    // 	mapSp = SepsetProducer(graph, sepsets, test);
    // 	sp = SepsetProducer(graph, test, possDsepSepsets, threads);
    
    // 	if (orientRule == ORIENT_SEPSETS) {

    // 	    mapSp.setOrientRule(orientRule);
    // 	    mapSp.setDepth(depth);
    // 	    mapSp.setVerbose(verbose);
    // 	    mapSp.setKnowledge(knowledge);
    // 	    if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
    // 	    mapSp.fillMap();

    // 	    fciorient_ = FciOrient(mapSp, whyOrient);
      
    // 	} else {

    // 	    sp.setOrientRule(orientRule);
    // 	    sp.setDepth(depth);
    // 	    sp.setVerbose(verbose);
    // 	    sp.setKnowledge(knowledge);
    // 	    if (verbose) Rcpp::Rcout << "    Filling Triple Map..." << std::endl;
      
    // 	    sp.fillMap();

    // 	    fciorient_ = FciOrient(sp, whyOrient);
    // 	}

    // 	fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
    // 	fciorient_.setMaxPathLength(maxPathLength);
    // 	fciorient_.setKnowledge(knowledge);

    // 	if (verbose) Rcpp::Rcout << "  Orienting colliders..." << std::endl;
    
    // 	fciorient_.ruleR0(graph);

    // 	if (orientRule == ORIENT_MAJORITY || orientRule == ORIENT_CONSERVATIVE) {
    // 	    for (auto t : sp.getAmbiguousTriples())
    // 		graph.addAmbiguousTriple(t.x, t.y, t.z);
    // 	}

    // 	if (verbose) Rcpp::Rcout << "  Orienting implied edges..." << std::endl;
    
    // 	fciorient_.doFinalOrientation(graph);
	
    // } else {

    // 	graph = fasGraph;
    // 	sepsets = fasSepsets;
	
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

std::vector<Node> Fci::getSepset(const Node& a, const Node& b, SepsetMap& sepsets) {

    if (sepsets.isInSepsetMap(a, b))
	return sepsets.get(a, b);
    
    std::vector<Node> sepset;
    double pval = 0;

    // if (a.isCensored() && b.isCensored()) {
    // 	sepsets.set(a, b, sepset, 1.0);
    // 	return sepset;
    // }
    
    std::vector<Node> adja = graph.getAdjacentNodes(a);
    std::vector<Node> adjb = graph.getAdjacentNodes(b);

    std::vector<Node> ppa = possibleParents(a, adja, b);
    std::vector<Node> ppb = possibleParents(b, adjb, a);

    int maxDepth = std::max(ppa.size(), ppb.size());
    maxDepth = std::min(maxDepth, std::max(depth, 1000));

    for (int d = 0; d <= maxDepth; d++) {

	if (d <= ppa.size()) {
	    ChoiceGenerator cg1(ppa.size(), d);
	    std::vector<int> *comb2;
	    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb2, ppa);
		double score;
		bool indep = test->isIndependent(a, b, s, &score);
		if (indep) {
		    pval = score;
		    sepset = s;
		}
	    }
	}

	if (d <= ppb.size()) {
	    ChoiceGenerator cg2(ppb.size(), d);
	    std::vector<int> *comb3;
	    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb3, ppb);
		double score;
		bool indep = test->isIndependent(a, b, s, &score);
		if (indep && (score > pval)) {
		    pval = score;
		    sepset = s;
		}
	    }
	}

	if (pval > 0) break;
    }

    sepsets.set(a, b, sepset, pval);
    
    return sepset;
    
}

std::vector<Node> Fci::getMinSepset(const Node& a, const Node& b,
				    std::vector<Node>& sepset, double *pReturn) {
  
    std::vector<Node> minSepset;
    double pval = 0;
    
    for (int d = 0; d <= sepset.size(); d++) {
	
        ChoiceGenerator cg(sepset.size(), d);
	std::vector<int> *comb;
	
	for (comb = cg.next(); comb != NULL; comb = cg.next()) {
	    std::vector<Node> s = GraphUtils::asList(*comb, sepset);
	    double score;
	    bool indep = test->isIndependent(a, b, s, &score);
	    if (indep) {
		pval = score;
		minSepset = s;
	    }
	}
	
	if (pval > 0) break;
    }

    if (pReturn != NULL && pval > 0) {
	*pReturn = pval;
    }
    
    // sepsets.set(a, b, minSepset, pval);
    
    return minSepset;
}

std::list<Triple> Fci::getMTriples(EdgeListGraph& graph) {
    std::list<Triple> mTriples;

    std::vector<Node> nodes = graph.getNodes();

    for (const Node& b : nodes) {
	std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);

	if (adjacentNodes.size() < 2) {
	    continue;
	}

	ChoiceGenerator cg(adjacentNodes.size(), 2);
	std::vector<int>  *combination;

	for (combination = cg.next(); combination != NULL; combination = cg.next()) {
	    Node a = adjacentNodes.at(combination->at(0));
	    Node c = adjacentNodes.at(combination->at(1));

	    // Skip triples that are shielded or where a and c are censored.
	    if (graph.isAdjacentTo(a, c)) {
		continue;
	    }
	    
	    mTriples.push_back(Triple(a, b, c));
	}
    }
    return mTriples;
}

SepsetMap Fci::ruleR0_rfciPrune(EdgeListGraph& graph, SepsetMap& sepsets) {
    std::list<Triple> mTriples = getMTriples(graph);

    if (verbose) RcppThread::Rcout << "  RFCI adjacency pruning...\n";
    
    SepsetMap rfciSepsets;

    while (!mTriples.empty()) {
	Triple t = mTriples.front();
	mTriples.pop_front();

	std::vector<Node> sepset = getSepset(t.x, t.z, sepsets);

	sepset.erase(std::remove(sepset.begin(), sepset.end(), t.y), sepset.end());
	
	bool indep1 = test->isIndependent(t.x, t.y, sepset);
	bool indep2 = test->isIndependent(t.z, t.y, sepset);

	if (indep1) {
	    double pval = 0;
	    std::vector<Node> minSepset = getMinSepset(t.x, t.y, sepset, &pval);
	    rfciSepsets.set(t.x, t.y, minSepset, pval);
	    graph.removeEdge(t.x, t.y);

	    if(verbose) {
		RcppThread::Rcout << "      Removed: " << t.x << " --- " << t.y << std::endl;
	    }

	    // RcppThread::Rcout << "      Removing: " << t.x << " --- " << t.y << " | { ";
	    // for (Node n : minSepset) {
	    // 	RcppThread::Rcout << n << " ";
	    // }
	    // RcppThread::Rcout << "}\n";

	}
	if (indep2) {
	    double pval = 0;
	    std::vector<Node> minSepset = getMinSepset(t.z, t.y, sepset, &pval);
	    rfciSepsets.set(t.z, t.y, minSepset, pval);
	    graph.removeEdge(t.z, t.y);

	    if(verbose) {
		RcppThread::Rcout << "      Removed: " << t.z << " --- " << t.y << std::endl;
	    }

	    // RcppThread::Rcout << "      Removing: " << t.z << " --- " << t.y << " | { ";
	    // for (Node n : minSepset) {
	    // 	RcppThread::Rcout << n << " ";
	    // }
	    // RcppThread::Rcout << "}\n";
	}

	if (indep1 || indep2) {

	    for (Node node : variables) {
		std::vector<Node> adj = graph.getAdjacentNodes(node);
		if (indep1) {
		    if (std::find(adj.begin(), adj.end(), t.x) != adj.end() &&
			std::find(adj.begin(), adj.end(), t.y) != adj.end()) {
			mTriples.push_back(Triple(t.x, node, t.y));
		    }
		}
		if (indep2) {
		    if (std::find(adj.begin(), adj.end(), t.z) != adj.end() &&
			std::find(adj.begin(), adj.end(), t.y) != adj.end()) {
			mTriples.push_back(Triple(t.z, node, t.y));
		    }
		}
	    }

	    auto it = mTriples.begin();

	    while (it != mTriples.end()) {
		if ((indep1 && it->contains(t.x) && it->contains(t.y)) ||
		    (indep2 && it->contains(t.z) && it->contains(t.y))) {
		    // RcppThread::Rcout << "Removing Triple: " << *it << std::endl;
		    it = mTriples.erase(it);
		} else {
		    it++;
		}
	    }
	
	}
	
    }

    return rfciSepsets;
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

std::vector<Node> Fci::possibleParents(const Node& x,
				       const std::vector<Node>& adjx,
				       const Node& y) {
    std::vector<Node> possParents;

    for (const Node& z : adjx) {
	if (z==x) continue;
	if (z==y) continue;

	if (possibleParentOf(x, z)) {
	    possParents.push_back(z);
	}
    }

    return possParents;
}

bool Fci::possibleParentOf(const Node& x, const Node& z) {
    return !knowledge.isForbidden(z, x) && !knowledge.isRequired(x, z);
}


/**
 * Orients according to background knowledge
 */
void Fci::fciOrientbk(/*IKnowledge bk,*/ EdgeListGraph& graph, std::vector<Node> variables) {

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
