#ifndef BOOTSTRAP_HPP_
#define BOOTSTRAP_HPP_

#include "StabilityUtils.hpp"
#include "armaLapack.hpp"

class Bootstrap {

private:

    DataSet d;
    std::string alg, ensemble;
    int B;
    int N;
    int subSize;
    std::vector<double> lambda;
    double alpha;
    double thresh;
    int iterLimit = 500;
    // EdgeListGraph* baseGraph = NULL;
    EdgeListGraph finalGraph;
    Rcpp::DataFrame stabs;
    arma::cube theta;
    // arma::umat subs;
    bool fdr;
    bool verbose = false;
    int threads = -1;

    arma::umat getBootstrapSamples();

    EdgeListGraph makeEnsembleGraph(std::vector<EdgeListGraph>& graphVec);


public:
    Bootstrap(DataSet& dat, std::string alg, std::string ensemble,
	      int numBoots, double thresh, double sampleFrac) :
        d(dat),
	alg(alg),
	ensemble(ensemble),
        B(numBoots),
	thresh(thresh),
        subSize(std::floor(sampleFrac * dat.getNumRows())),
	N(dat.getNumRows()) {}

    EdgeListGraph runBootstrap();

    // void setComputeStabs(bool cs) { computeStabs = cs; }
    // bool getComputeStabs() { return computeStabs; }

    // arma::mat getStabs() { return stabilities; }

    void setVerbose(bool v) { verbose = v; }
    void setThreads(int threads) { this->threads = threads; }
    void setAlpha(double a) { alpha = a; }
    void setLambda(std::vector<double> l) { lambda = l; }
    void setAlgorithm(std::string alg) { this->alg = alg; }
    void setFdr(bool fdr) { this->fdr = fdr; }
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
