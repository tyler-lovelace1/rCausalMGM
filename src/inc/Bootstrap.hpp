#ifndef BOOTSTRAP_HPP_
#define BOOTSTRAP_HPP_

#include "MGM.hpp"
#include "CoxMGM.hpp"
#include "StabilityUtils.hpp"
#include "SepsetMap.hpp"
#include "IndTestMultiCox.hpp"
#include "PcStable.hpp"
#include "Knowledge.hpp"
#include "Fci.hpp"

class Bootstrap {

private:

    DataSet d;
    std::string alg;
    int B;
    int N;
    arma::umat samps;
    std::vector<double> lambda;
    double alpha;

    EdgeListGraph finalGraph;
    Rcpp::DataFrame stabs;
    
    OrientRule orientRule;
    Knowledge knowledge;
    
    bool replace = false;
    bool verbose = false;
    int threads = -1;

    EdgeListGraph makeEnsembleGraph(std::vector<EdgeListGraph>& graphVec);


public:
    Bootstrap(DataSet& dat, std::string alg, int numBoots, bool replace) :
        d(dat),
	alg(alg),
        B(numBoots),
	replace(replace) {
	if (replace)   N = dat.getNumRows();
	else           N = (int) std::ceil(0.632 * dat.getNumRows());
    }

    EdgeListGraph runBootstrap();

    void setVerbose(bool v) { verbose = v; }
    void setThreads(int threads) { this->threads = threads; }
    void setAlpha(double a) { alpha = a; }
    void setLambda(std::vector<double> l) { lambda = l; }
    void setAlgorithm(std::string alg) { this->alg = alg; }
    void setOrientRule(OrientRule rule) {this->orientRule = rule; }
    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }


    arma::umat getSubSamps() { return samps; }

    // void setFdr(bool fdr) { this->fdr = fdr; }
    // void setBaseGraph(EdgeListGraph* bg) {
    // 	baseGraph = bg;
    // 	alg = bg->getAlgorithm();
    // 	if (!bg->getHyperParam("lambda").isNull()) {
    // 	    lambda = Rcpp::as<std::vector<double>>(bg->getHyperParam("lambda"));
    // 	}
    // 	if (!bg->getHyperParam("alpha").isNull()) {
    // 	    alpha = Rcpp::as<std::vector<double>>(bg->getHyperParam("alpha"))[0];
    // 	}
    // }

    Rcpp::DataFrame getStabs() { return stabs; }

};

#endif /* BOOTSTRAP_HPP_ */
