#ifndef SEARCHCV_HPP_
#define SEARCHCV_HPP_

#include "armaLapack.hpp"
#include "LogisticRegression.hpp"
#include "LinearRegression.hpp"
#include "CoxIRLSRegression.hpp"
#include "DataSet.hpp"
#include "MGM.hpp"
#include "PcStable.hpp"
#include "FciStable.hpp"


class SearchCV {

    DataSet originalFull;
    DataSet internalFull;

    std::vector<DataSet> trainVec;
    std::vector<DataSet> testVec;

    arma::vec lambdas;
    arma::vec alphas;

    uint nfolds;

    arma::uvec foldid;

public:

    SearchCV() {}
    SearchCV(const DataSet& data, uint nfolds);
    SearchCV(const DataSet& data, const arma::uvec& foldid);

    std::vector<EdgeListGraph> MgmCV(std::vector<double> lambdas, arma::mat& loglik);

    std::vector<EdgeListGraph> PcStableCV(std::vector<double> alphas,
					  arma::mat& loglik);

    std::vector<EdgeListGraph> FciStableCV(std::vector<double> alphas,
					   arma::mat& loglik);

    std::vector<EdgeListGraph> PcStableGridCV(std::vector<double> lambdas,
					      std::vector<double> alphas,
					      arma::cube& loglik);

    std::vector<EdgeListGraph> FciStableGridCV(std::vector<double> lambdas,
					       std::vector<double> alphas,
					       arma::cube& loglik);


};

#endif /* SEARCHCV_HPP_ */
