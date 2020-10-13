#ifndef STABILITYUTILS_HPP_
#define STABILITYUTILS_HPP_

#include "MGM.hpp"
#include <RcppArmadillo.h>

class StabilityUtils {

public:

    static arma::mat stabilitySearchPar(DataSet& data, MGM& gs, int N, int b);

    static arma::mat stabilitySearchPar(DataSet& data, MGM& gs); // LOO

    static arma::mat stabilitySearchPar(DataSet& data, MGM& gs, arma::umat subs);

    static arma::umat subSampleNoReplacement(int sampSize, int subSize, int numSub);

    static int getSubSize(int samplesize);

    static arma::mat skeletonToMatrix(EdgeListGraph& graph, DataSet& d);

    static int checkForVariance(DataSet& d, DataSet& full);

};

#endif /* STABILITYUTILS_HPP_ */