#ifndef STABILITYUTILS_HPP_
#define STABILITYUTILS_HPP_

#include "MGM.hpp"
#include <RcppArmadillo.h>
#include "BlockingQueue.hpp"
#include <thread>

class StabilityUtils {

public:

    static arma::mat stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int N, int b);

    static arma::mat stabilitySearchPar(DataSet& data, std::vector<double>& lambda); // LOO

    static arma::mat stabilitySearchPar(DataSet& data, std::vector<double>& lambda, arma::umat& subs);

    static arma::umat subSampleNoReplacement(int sampSize, int subSize, int numSub);

    static arma::urowvec subSampleIndices(int N, int subSize);

    static int getSubSize(int samplesize);

    static arma::mat skeletonToMatrix(EdgeListGraph& graph, DataSet& d);

    static int checkForVariance(DataSet& d, DataSet& full);

    static arma::vec standardizeData(const arma::vec& data);

};

#endif /* STABILITYUTILS_HPP_ */