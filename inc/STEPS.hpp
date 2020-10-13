#ifndef STEPS_HPP_
#define STEPS_HPP_

#include "StabilityUtils.hpp"
#include <RcppArmadillo.h>

class STEPS {

private:

    DataSet d;
    int N;
    int b;
    std::vector<double> lambda;
    double gamma;
    bool includeZeros = true;
    int iterLimit = 500;
    double origLambda;
    // Graph pdGraph;
    // Graph lastGraph;
    std::vector<double> lastLambda;
    bool leaveOneOut = false;
    arma::mat stabilities;
    arma::umat subs;
    bool computeStabs;

public:
    STEPS(DataSet& dat, std::vector<double>& lam, double g, int numSub, bool loo) :
        d(dat),
        leaveOneOut(loo),
        N(numSub),
        gamma(g),
        lambda(lam),
        b(StabilityUtils::getSubSize(dat.getNumRows())) {}

    EdgeListGraph runStepsPar();



};

#endif /* STEPS_HPP_ */