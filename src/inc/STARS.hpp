#ifndef STARS_HPP_
#define STARS_HPP_

#include "StabilityUtils.hpp"
#include <RcppArmadillo.h>

class STARS {

private:

    DataSet d;
    std::string alg;
    int N;
    int b;
    std::vector<double> params;
    double gamma;
    // bool includeZeros = true;
    int iterLimit = 500;
    // Graph pdGraph;
    // Graph lastGraph;
    EdgeListGraph* initialGraph = NULL;
    double lastParam;
    bool adjacency = true;
    bool leaveOneOut = false;
    // arma::cube stabilities;
    // arma::umat subs;
    // bool computeStabs = false;
    bool verbose = false;
    int threads = -1;

public:
    STARS(DataSet& dat, std::string alg,
	  std::vector<double>& params, double g,
	  int numSub, bool adj = true, bool loo = false) :
        d(dat),
	alg(alg),
	adjacency(adj),
        leaveOneOut(loo),
        N(numSub),
        gamma(g),
        params(params),
        b(StabilityUtils::getSubSize(dat.getNumRows())) {}

    EdgeListGraph runStarsPar();

    // void setComputeStabs(bool cs) { computeStabs = cs; }
    // bool getComputeStabs() { return computeStabs; }

    // arma::mat getStabs() { return stabilities; }

    void setVerbose(bool v) { verbose = v; }
    void setThreads(int threads) { this->threads = threads; }
    void setInitialGraph(EdgeListGraph* ig) { initialGraph = ig; }

};

#endif /* STARS_HPP_ */
