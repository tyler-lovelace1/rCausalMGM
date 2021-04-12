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
    std::vector<Variable*> vars = test->getVariables();
    this->variables.insert(variables.end(), vars.begin(), vars.end());
    buildIndexing(vars);
}

/**
 * Constructs a new FCI search for the given independence test and background knowledge and a list of variables to
 * search over.
 */
FciMax::FciMax(IndependenceTest* test, std::vector<Variable*> searchVars) {
    if (test == NULL /*|| knowledge == null*/) {
        throw std::invalid_argument("independenceTest cannot be null");
    }

    this->test = test;
    std::vector<Variable*> vars = test->getVariables();
    this->variables.insert(variables.end(), vars.begin(), vars.end());

    std::unordered_set<Variable*> remVars;
    for (Variable* node1 : this->variables) {
        bool search = false;
        for (Variable* node2 : searchVars) {
            if (node1->getName() == node2->getName()) {
                search = true;
            }
        }
        if (!search) {
            remVars.insert(node1);
        }
    }
    for (Variable* var: remVars) {
      std::remove(this->variables.begin(), this->variables.end(), var);
    }
}

//========================PUBLIC METHODS==========================//

EdgeListGraph FciMax::search() {
    return search(test->getVariables());
}

EdgeListGraph FciMax::search(const std::vector<Variable*>& nodes) {
    FasStableProducerConsumer fas(initialGraph, test);

    return search(fas, nodes);
}

EdgeListGraph FciMax::search(FasStableProducerConsumer& fas, const std::vector<Variable*>& nodes) {

    auto startTime = std::chrono::high_resolution_clock::now();
    
    if (verbose) Rcpp::Rcout << "Starting FCI-Max algorithm..." << std::endl;

    whyOrient = std::unordered_map<std::string, std::string>();

    if (test == NULL)
        throw std::invalid_argument("independenceTest of FCI-Max may not be NULL.");

    // fas->setKnowledge(getKnowledge());
    fas.setDepth(depth);
    fas.setVerbose(verbose);
    graph = fas.search();
    sepsets = fas.getSepsets();
    graph.reorientAllWith(ENDPOINT_CIRCLE);

    // SepsetsPossibleDsep sp(graph, test,/* knowledge,*/ depth, maxPathLength);

    // The original FCI, with or without JiJi Zhang's orientation rules
    // Optional step: Possible Dsep. (Needed for correctness but very time consuming.)
    if (isPossibleDsepSearchDone()) {
          if (verbose) Rcpp::Rcout << "Starting Posssible DSep search..." << std::endl;
          // SepsetsSet ssset(sepsets, test);
          // FciMaxOrient orienter(&ssset);

          PossibleDsepFciConsumerProducer possibleDSep(graph, test);
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

    sepsetsMaxP = (SepsetProducer*) new SepsetProducerMaxP(graph, test, sepsets);
    sepsetsMaxP->setDepth(depth);
    sepsetsMaxP->setVerbose(verbose);
    // sepsetsMaxP->fillMap();

    // Step CI C (Zhang's step F3.)
    //fciOrientbk(getKnowledge(), graph, independenceTest.getVariables());    - Robert Tillman 2008
    //        fciOrientbk(getKnowledge(), graph, variables);
    //        new FciMaxOrient(graph, new Sepsets(this.sepsets)).ruleR0(new Sepsets(this.sepsets));
    if (verbose) Rcpp::Rcout << "Starting Orientations..." << std::endl;
    //SepsetsSet sepsetsset_(sepsets, test);
    FciOrient fciorient_(sepsetsMaxP, whyOrient);
    fciorient_.setCompleteRuleSetUsed(completeRuleSetUsed);
    fciorient_.setMaxPathLength(maxPathLength);
    // fciOrient.setKnowledge(knowledge);
    // if (verbose) Rcpp::Rcout << "FCIOrient Set Up" << std::endl;
    fciorient_.ruleR0(graph);
    // if (verbose) Rcpp::Rcout << "Rule 0 finished" << std::endl;
    fciorient_.doFinalOrientation(graph);

    // Set algorithm and type
    std::ostringstream alg;
    alg << "FCI-Max: alpha = " << independenceTest->getAlpha();
    graph.setAlgorithm(alg.str());
    graph.setGraphType("partial ancestral graph");
    
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count() / 1000.0;
    
    if (verbose) Rcpp::Rcout << "FCI-Max finished:  " << std::round(elapsedTime) << " s" << std::endl;
    
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

void FciMax::buildIndexing(std::vector<Variable*> nodes) {
    this->hashIndices =  std::unordered_map<Variable*, int>();
    for (Variable* node : nodes) {
        auto itr = find(variables.begin(), variables.end(), node);
        int index = std::distance(variables.begin(), itr);
        this->hashIndices.insert(std::pair<Variable*, int>(node, index));
    }
}

/**
 * Orients according to background knowledge
 */
void FciMax::fciOrientbk(/*IKnowledge bk,*/ EdgeListGraph graph, std::vector<Variable*> variables) {

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
