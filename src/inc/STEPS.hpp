#ifndef STEPS_HPP_
#define STEPS_HPP_

#include "StabilityUtils.hpp"
#include "armaLapack.hpp"

class STEPS {

private:

    DataSet d;
    int N;
    int b;
    std::vector<double> lambda;
    double gamma;
    // bool includeZeros = true;
    int iterLimit = 500;
    double origLambda;
    // Graph pdGraph;
    // Graph lastGraph;
    std::vector<double> lastLambda;
    bool leaveOneOut = false;
    arma::mat stabilities;
    // arma::umat subs;
    bool computeStabs = false;
    bool verbose = false;
    int threads = -1;

public:
    STEPS() {}
    
    STEPS(DataSet& dat, std::vector<double>& lam, double g, int numSub, bool loo = false) :
        d(dat),
        leaveOneOut(loo),
        N(numSub),
        gamma(g),
        lambda(lam),
        b(StabilityUtils::getSubSize(dat.getNumRows())) {}

    STEPS(DataSet& dat, std::vector<double>& lam, double g, int numSub, int subSize, bool loo = false) :
        d(dat),
        leaveOneOut(loo),
        N(numSub),
        gamma(g),
        lambda(lam),
        b(subSize) {}

    EdgeListGraph runStepsPar();

    void setComputeStabs(bool cs) { computeStabs = cs; }
    bool getComputeStabs() { return computeStabs; }

    arma::mat getStabs() { return stabilities; }

    void setVerbose(bool v) { verbose = v; }
    void setThreads(int threads) { this->threads = threads; }

};

#endif /* STEPS_HPP_ */
