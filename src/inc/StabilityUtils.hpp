#ifndef STABILITYUTILS_HPP_
#define STABILITYUTILS_HPP_

// [[Rcpp::depends(RcppThread)]]

#include "MGM.hpp"
#include "CoxMGM.hpp"
// #include "PcStable.hpp"
// #include "CpcStable.hpp"
// #include "PcMax.hpp"
// #include "Fci.hpp"
// #include "Cfci.hpp"
// #include "FciMax.hpp"
// #include "IndTestMulti.hpp"
#include "RcppThread.h"
#include "armaLapack.hpp"
#include "BlockingQueue.hpp"
#include <thread>

class StabilityUtils {

public:

    static arma::mat stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int num_threads, int N, int b);

    static arma::mat stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int num_threads); // LOO

    static arma::mat stabilitySearchPar(DataSet& data, std::vector<double>& lambda, int num_threads, arma::umat& subs);

    static double stabilitySearchStars(DataSet& data, std::string& alg, double param, EdgeListGraph* initialGraph, int num_threads, bool adjcacency, int N, int b);

    static double stabilitySearchStars(DataSet& data, std::string& alg, double param, EdgeListGraph* initialGraph, int num_threads, bool adjcacency); // LOO

    static double stabilitySearchStars(DataSet& data, std::string& alg, double param, EdgeListGraph* initialGraph, int num_threads, bool adjcacency, arma::umat& subs);

    static arma::umat subSampleNoReplacement(DataSet& data, int subSize, int numSub);

    static arma::umat subSampleWithReplacement(DataSet& data, int subSize, int numSub);

    static arma::umat subSampleLOO(DataSet& data);

    // static arma::urowvec subSampleIndices(int N, int subSize);

    static int getSubSize(int samplesize);

    static arma::mat skeletonToMatrix(EdgeListGraph& graph, DataSet& d);
    
    static arma::cube graphToCube(EdgeListGraph& graph, DataSet& d);
    
    static int checkForVariance(DataSet& d, DataSet& full);

    static int checkForVariance(DataSet& d);

    static arma::vec standardizeData(const arma::vec& data);

};

#endif /* STABILITYUTILS_HPP_ */
