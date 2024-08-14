#ifndef FCI_HPP_
#define FCI_HPP_

/**
 * Extends Erin Korber's implementation of the Fast Causal Inference algorithm (found in FCI.java) with Jiji Zhang's
 * Augmented FCI rules (found in sec. 4.1 of Zhang's 2006 PhD dissertation, "Causal Inference and Reasoning in Causally
 * Insufficient Systems").
 * <p>
 * This class is based off a copy of FCI.java taken from the repository on 2008/12/16, revision 7306. The extension is
 * done by extending doFinalOrientation() with methods for Zhang's rules R5-R10 which implements the augmented search.
 * (By a remark of Zhang's, the rule applications can be staged in this way.)
 *
 * @author Erin Korber, June 2004
 * @author Alex Smith, December 2008
 * @author Joseph Ramsey
 * @author Choh-Man Teng
 * @author Jack Fiore, November 2020
 */

#include "EdgeListGraph.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"
#include "ChoiceGenerator.hpp"
#include "Knowledge.hpp"
#include "SepsetProducer.hpp"
#include "FasStableProducerConsumer.hpp"
#include "PossibleDsepFciConsumerProducer.hpp"
// #include "FasStableConcurrent.hpp"
#include "FciOrient.hpp"
// #include "SepsetsPossibleDsep.hpp"
// #include "SepsetsSet.hpp"
#include <limits>


class Fci {

private:

    /**
     * The PAG being constructed.
     */
    EdgeListGraph graph;

    EdgeListGraph fasGraph;

    /**
     * The SepsetMap being constructed.
     */
    SepsetMap sepsets;

    /**
     * The SepsetMap being constructed.
     */
    SepsetMap fasSepsets;

    /**
     * The background knowledge.
     */
    Knowledge knowledge;

    /**
     * The variables to search over (optional)
     */
    std::vector<Node> variables;

    IndependenceTest *test;

    /**
     * flag for complete rule set, true if should use complete rule set, false otherwise.
     */
    bool completeRuleSetUsed = false;

    /**
     * True iff the possible dsep search is done.
     */
    bool possibleDsepSearchDone = true;

    /** 
     * Orientation rule for identifying colliders
     */
    OrientRule orientRule = ORIENT_MAJORITY;

    /**
     * The maximum length for any discriminating path. -1 if unlimited; otherwise, a positive integer.
     */
    int maxPathLength = -1;

    /**
     * The depth for the fast adjacency search.
     */
    int depth = -1;

    /**
     * The number of consumer threads to create for multi-threaded steps. -1 to set automatically
     */ 
    int threads = -1;

    /**
     * Elapsed time of last search.
     */
    long elapsedTime;

    std::unordered_map<std::string,std::string> whyOrient;

    /*
        Maps from "<Node>,<Node>" to integer, where integer describes which rule oriented node -> node as an arrowhead
        Useful for determining the relative efficacy of each rule in finding latent confounders
     */

    bool verbose = false;
    bool fdr = false;
    // EdgeListGraph truePag;
    std::unordered_map<Node, int> hashIndices;   // must lock
    // ICovarianceMatrix covarianceMatrix;
    double penaltyDiscount = 2;
    // SepsetMap possibleDsepSepsets, nullSepsets;
    EdgeListGraph *initialGraph = NULL;
    int possibleDsepDepth = -1;

    SepsetProducer sp, mapSp, fasMapSp, posDsp, posDmapSp;
    SepsetMap possDsepSepsets, nullSepsets;

        //========================PRIVATE METHODS==========================//

    void buildIndexing(std::vector<Node> nodes);

    /**
     * Orients according to background knowledge
     */
    void fciOrientbk(EdgeListGraph& graph, std::vector<Node> variables);

    SepsetMap ruleR0_rfciPrune(EdgeListGraph& graph, SepsetMap& sepsets);

    std::list<Triple> getMTriples(EdgeListGraph& graph);

    std::vector<Node> getSepset(const Node& a, const Node& b, SepsetMap& sepsets);

    std::vector<Node> getMinSepset(const Node& a, const Node& b, std::vector<Node>& sepset, double* pReturn = NULL);

    std::vector<Node> possibleParents(const Node& x,
				      const std::vector<Node>& adjx,
				      const Node& y);
  
    bool possibleParentOf(const Node& x, const Node& z);


public:
        //============================CONSTRUCTORS============================//
    /**
     * Constructs a new FCI search for the given independence test and background knowledge.
     */
    Fci(IndependenceTest *test);

    /**
     * Constructs a new FCI search for the given independence test and background knowledge and a list of variables to
     * search over.
     */
    Fci(IndependenceTest *test, std::vector<Node> searchVars);

       //========================PUBLIC METHODS==========================//

    int getDepth() { return depth; }

    void setDepth(int depth);

    long getElapsedTime() { return this->elapsedTime; }

    void setInitialGraph(EdgeListGraph *initialGraph) { this->initialGraph = initialGraph; }

    EdgeListGraph search();

    EdgeListGraph search(const std::vector<Node>& nodes);

    EdgeListGraph search(FasStableProducerConsumer& fas, const std::vector<Node>& nodes);

    EdgeListGraph reorientWithRule(OrientRule rule);

    SepsetMap getSepsets() { return this->sepsets; }

    Knowledge getKnowledge() { return this->knowledge; }

    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }

    OrientRule getOrientRule() { return orientRule; }
    
    void setOrientRule(OrientRule rule) {this->orientRule = rule; }

    /**
     * @return true if Zhang's complete rule set should be used, false if only R1-R4 (the rule set of the original FCI)
     * should be used. False by default.
     */
    bool isCompleteRuleSetUsed() { return completeRuleSetUsed; }

    /**
     * @param completeRuleSetUsed set to true if Zhang's complete rule set should be used, false if only R1-R4 (the rule
     *                            set of the original FCI) should be used. False by default.
     */
    void setCompleteRuleSetUsed(bool completeRuleSetUsed) { this->completeRuleSetUsed = completeRuleSetUsed; }

    bool isPossibleDsepSearchDone() { return possibleDsepSearchDone; }

    void setPossibleDsepSearchDone(bool possibleDsepSearchDone) { this->possibleDsepSearchDone = possibleDsepSearchDone; }

    /**
     * @return the maximum length of any discriminating path, or -1 of unlimited.
     */
    int getMaxPathLength() { return (maxPathLength == std::numeric_limits<int>::max()) ? -1 : maxPathLength; }

    /**
     * @param maxPathLength the maximum length of any discriminating path, or -1 if unlimited.
     */
    void setMaxPathLength(int maxPathLength);

    // /**
    //  * The independence test.
    //  */
    // IndependenceTest getIndependenceTest() { return test; }

    // void setTruePag(EdgeListGraph& truePag) { this->truePag = truePag; }

    double getPenaltyDiscount() { return penaltyDiscount; }

    void setPenaltyDiscount(double penaltyDiscount) { this->penaltyDiscount = penaltyDiscount; }

    int getPossibleDsepDepth() { return possibleDsepDepth; }

    void setPossibleDsepDepth(int possibleDsepDepth) { this->possibleDsepDepth = possibleDsepDepth; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setFDR(bool fdr) { this->fdr = fdr; }

    void setThreads(int threads) { this->threads = threads; }
};

#endif /* FCI_HPP_ */
