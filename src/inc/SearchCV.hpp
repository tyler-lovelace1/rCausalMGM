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

    DataSet originalData;
    DataSet internalData;

    LinearRegression regression;
    LogisticRegression logisticRegression;
    CoxIRLSRegression coxRegression;

    arma::vec lambdas;
    arma::vec alphas;

    uint nfolds;

    arma::uvec foldid;

    std::vector<Node> scoreNodes;

    bool verbose = false;

    std::vector<Node> expandVariable(DataSet& dataSet, const Node& var);
    double multiLL(arma::mat &coeffs, const Node& dep, std::vector<Node>& indep, int k);

public:

    SearchCV() {}
    SearchCV(const DataSet& data, uint nfolds);
    SearchCV(const DataSet& data, const arma::uvec& foldid);

    std::vector<EdgeListGraph> MgmCV(std::vector<double> lambdas, arma::mat& loglik);

    std::vector<EdgeListGraph> PcStableCV(arma::mat& loglik);

    std::vector<EdgeListGraph> FciStableCV(arma::mat& loglik);

    double scoreNode(EdgeListGraph graph, Node target, int k);
};

#endif /* SEARCHCV_HPP_ */
