#ifndef SEARCHCV_HPP_
#define SEARCHCV_HPP_

#include "armaLapack.hpp"
#include "RcppThread.h"
#include "DataSet.hpp"
#include "MGM.hpp"
#include "CoxMGM.hpp"
#include "PcStable.hpp"
#include "Fci.hpp"
#include "IndTestMultiCox.hpp"
#include "LinearRegression.hpp"
#include "LogisticRegression.hpp"
#include "CoxRegression.hpp"

struct CvResult;

class SearchCV {

    DataSet originalData;
    DataSet internalData;

    LinearRegression regression;
    LogisticRegression logisticRegression;
    CoxRegression coxRegression;

    // Hyperparameter vectors
    arma::vec lambdas;                       // p values of lambda
    arma::vec alphas;                        // q values of alpha
    std::vector<OrientRule> orientRules;     // r values of orientation rule

    std::string alg = "pc"; // "pc" or "fci"

    bool fdr = false;

    uint nfolds; // number of folds K

    arma::uvec foldid; // fold ids in range 1:K

    std::vector<Node> scoreNodes; // nodes to include in scoring, all nodes by default
    std::map<Node, std::vector<Node>> variablesPerNode;

    bool verbose = false;

    int parallelism = -1;

    EdgeListGraph *initialGraph = NULL;

    Knowledge knowledge;

    // arma::cube llMean;  // p x q x r
    // arma::cube llSe;    // p x q x r

    std::set<CvResult> results;

    std::vector<Node> expandVariable(DataSet& dataSet, const Node& var);
    
    double multiTestLL(arma::mat &coeffs, const Node& dep,
		       std::vector<Node>& indep, arma::uvec testRows);
    
    double scoreTestLLTask(const Node& dep, std::vector<Node>& indep, int k);

    std::vector<double> scoreGraphTestLL(EdgeListGraph graph, int k);

    // bool checkFoldID(const arma::uvec& foldid);

public:

    SearchCV() {}
    SearchCV(DataSet& data, std::string alg, uint nfolds, int threads = -1);
    SearchCV(DataSet& data, std::string alg, const arma::uvec& foldid, int threads = -1);

    void setInitialGraph(EdgeListGraph *initialGraph) { this->initialGraph = initialGraph; }
    void setLambdas(arma::vec lambdas) { this->lambdas = lambdas; }
    void setAlphas(arma::vec alphas) { this->alphas = alphas; }
    void setOrientRules(std::vector<OrientRule> orientRules) { this->orientRules = orientRules; }
    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }
    void setVerbose(bool verbose) { this->verbose = verbose; }
    void setFDR(bool fdr) { this->fdr = fdr; }

    std::set<CvResult> getCVResults() { return results; };

    arma::uvec getFoldID() { return foldid; }

    // std::vector<EdgeListGraph> MgmCV(std::vector<double> lambdas, arma::mat& loglik);

    std::vector<EdgeListGraph> causalCV();

    std::vector<EdgeListGraph> causalMGMGridCV();

    std::vector<EdgeListGraph> causalMGMRandCV();

    static bool checkFoldID(const arma::uvec& foldid);

    // double scoreNode(EdgeListGraph graph, Node target, int k);

    // friend arma::vec scoreGraphTest1(const Rcpp::DataFrame& data, Rcpp::List graph, arma::uvec foldid);
    // friend arma::vec scoreGraphTest2(const Rcpp::DataFrame& data, Rcpp::List graph, uint nfolds);
    
    // friend Rcpp::List pcCvTest1(const Rcpp::DataFrame& data,
    // 				arma::vec alphas, arma::uvec foldid,
    // 				Rcpp::Nullable<Rcpp::List> initialGraph);
    
    // friend Rcpp::List pcCvTest2(const Rcpp::DataFrame& data,
    // 				arma::vec alphas, uint nfolds,
    // 				Rcpp::Nullable<Rcpp::List> initialGraph);

    // friend Rcpp::List pcGridCvTest2(const Rcpp::DataFrame& data, arma::vec lambdas,
    // 				    arma::vec alphas, uint nfolds);

    friend Rcpp::List causalMGMRandCvTest(const Rcpp::DataFrame& data,
					  std::string alg,
					  arma::vec lambdas,
					  arma::vec alphas,
					  uint nfolds,
					  bool fdr);
};


struct CvResult {
    double mean;
    double se;
    double mbSize;
    double alpha;
    std::vector<double> lambda;
    OrientRule rule;

    CvResult() {}

    CvResult(double mean, double se, double mbSize) : mean(mean), se(se), mbSize(mbSize) {}

    CvResult(double mean, double se, double mbSize, double alpha) : mean(mean), se(se), mbSize(mbSize), alpha(alpha) {}

    CvResult(double mean, double se, double mbSize, double alpha, OrientRule rule) : mean(mean), se(se), mbSize(mbSize), alpha(alpha), rule(rule) {}

    CvResult(double mean, double se, double mbSize, double alpha, std::vector<double> lambda, OrientRule rule) : mean(mean), se(se), mbSize(mbSize), alpha(alpha), lambda(lambda), rule(rule) {}

    bool operator==(const CvResult& rhs) const {
	return (mean == rhs.mean) && (se == rhs.se) && (mbSize == rhs.mbSize) &&
	    (alpha == rhs.alpha) && (lambda == rhs.lambda) && (rule == rhs.rule);
    }
	
    bool operator<(const CvResult& rhs) const {
	if (mbSize == rhs.mbSize) {
	    if (alpha == rhs.alpha) {
		    
		double lambdaSum = 0.0;
		double rhsLambdaSum = 0.0;
		for (int i = 0; i < lambda.size(); i++) {
		    lambdaSum += lambda.at(i);
		    rhsLambdaSum += rhs.lambda.at(i);
		}
		    
		if (lambdaSum == rhsLambdaSum) {
		    if (rule == rhs.rule) {
			if (se == rhs.se) {
			    return mean < rhs.mean;
			}
			return se < rhs.se;
		    }
			
		    std::vector<OrientRule> ruleOrder = { ORIENT_CONSERVATIVE,
							  ORIENT_MAJORITY,
							  ORIENT_MAXP,
							  ORIENT_SEPSETS };
		    int ruleIdx = 4, rhsRuleIdx = 4;
		    for (int i = 0; i < ruleOrder.size(); i++) {
			if (rule == ruleOrder[i]) ruleIdx = i;
			if (rhs.rule == ruleOrder[i]) rhsRuleIdx = i;
		    }
			    
		    return ruleIdx < rhsRuleIdx;
		}
		return lambdaSum > rhsLambdaSum;
	    }
	    return alpha < rhs.alpha;
	}
	return mbSize < rhs.mbSize;
    }
	
    bool operator!=(const CvResult& rhs) const {
	return !(*this == rhs);
    }
	
    bool operator> (const CvResult& rhs) const {
	return rhs < *this;
    }
	
    bool operator<=(const CvResult& rhs) const {
	return !(*this > rhs);
    }
	
    bool operator>=(const CvResult& rhs) const {
	return !(*this < rhs);
    }
};


#endif /* SEARCHCV_HPP_ */

