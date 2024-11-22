#ifndef REGRESSIONBICSCORE_HPP_
#define REGRESSIONBICSCORE_HPP_

#include "armaLapack.hpp"
#include "Score.hpp"
#include "DataSet.hpp"
#include "LogisticRegression.hpp"
#include "LinearRegression.hpp"
// #include "CoxRegression.hpp"

class RegressionBicScore : public Score {
private:

    DataSet originalData;
    std::vector<Node> searchVariables;
    DataSet internalData;
    double penalty;
    std::map<Node, std::vector<Node>> variablesPerNode;
    std::map<Node, double> nullLL;
    bool verbose = false;

    std::map<Node, double> logNevents;

    double N, logN;

    // CoxRegression coxRegression;
    LogisticRegression logisticRegression;
    LinearRegression regression;
    
    std::map<std::pair<Node,Node>, arma::vec> WZmap;

    // arma::uvec getIndices(const std::vector<Node>& nodes);

    std::vector<Node> expandVariable(DataSet& dataSet, const Node& var);

    arma::mat getSubsetData(DataSet& origData, std::vector<Node>& varSubset);

    // void resetWZ(Node target, std::vector<Node>& neighbors);

    double multiLL(arma::mat &coeffs, const Node& dep, std::vector<Node>& indep);

    double logLikLinearRegression(const Node& x, std::vector<Node>& regressors);

    double logLikMultinomialLogisticRegression(const Node& x, std::vector<Node>& regressors);

    // double logLikCoxRegression(const Node& x, std::vector<Node>& regressors);

public:

    RegressionBicScore() {}

    RegressionBicScore(DataSet& data, double penalty);

    RegressionBicScore(const RegressionBicScore& other) = default;
    RegressionBicScore& operator=(const RegressionBicScore& other) = default;

    RegressionBicScore(RegressionBicScore&& other) = default;
    RegressionBicScore& operator=(RegressionBicScore&& other) = default;

    std::vector<Node> getVariables() {return this->searchVariables; }

    double localScore(const Node& x, const std::vector<Node>& z);

    /**
     * @return the list of variable varNames.
     */
    std::vector<std::string> getVariableNames();

    Node getVariable(std::string name);

    double getPenalty() { return this->penalty; }

    void setPenalty(double penalty) { this->penalty = penalty; }

    DataSet getData() { return this->originalData; }

    bool isVerbose() { return this->verbose; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    int getSampleSize() { return originalData.getNumRows(); }

    std::vector<DataSet*> getDataSets() { return { &originalData }; }

    friend double DGScoreTest(const Rcpp::DataFrame& df,
			      std::string targetName,
			      std::vector<std::string>& regressorNames);

    
};


#endif /* REGRESSIONBICSCORE_HPP_ */
