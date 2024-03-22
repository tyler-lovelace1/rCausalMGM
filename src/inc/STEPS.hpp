#ifndef STEPS_HPP_
#define STEPS_HPP_

#include "armaLapack.hpp"
#include "StabilityUtils.hpp"

class STEPS {

private:

    DataSet d;
    int N;
    int b;
    std::vector<double> lambda;
    double gamma;
    int iterLimit = 500;
    double origLambda;
    std::vector<double> lastLambda;
    bool leaveOneOut = false;
    arma::mat stabilities;
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

    std::vector<EdgeListGraph> runStepsPath(arma::mat& instabs, arma::umat& samps);

    void setComputeStabs(bool cs) { computeStabs = cs; }
    bool getComputeStabs() { return computeStabs; }

    arma::mat getStabs() { return stabilities; }

    void setVerbose(bool v) { verbose = v; }
    void setThreads(int threads) { this->threads = threads; }

};

#endif /* STEPS_HPP_ */
