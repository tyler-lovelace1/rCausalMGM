#ifndef DEGENERATEGAUSSIANSCORE_HPP_
#define DEGENERATEGAUSSIANSCORE_HPP_

#include "armaLapack.hpp"
#include "Score.hpp"
#include "DataSet.hpp"

class DegenerateGaussianScore : public Score {
private:

    // DataSet originalData;
    std::vector<Node> searchVariables;
    // DataSet internalData;
    double penalty;
    std::map<Node, std::vector<Node>> variablesPerNode;
    std::map<Node, int> indexMap; 
    bool verbose = false;

    double N, logN;

    arma::mat covMat;

    arma::uvec getIndices(const std::vector<Node>& nodes);

    std::vector<Node> expandVariable(DataSet& dataSet, const Node& var);

public:

    DegenerateGaussianScore() {}

    DegenerateGaussianScore(DataSet data, double penalty);

    DegenerateGaussianScore(const DegenerateGaussianScore& other) = default;
    DegenerateGaussianScore& operator=(const DegenerateGaussianScore& other) = default;

    DegenerateGaussianScore(DegenerateGaussianScore&& other) = default;
    DegenerateGaussianScore& operator=(DegenerateGaussianScore&& other) = default;

    std::vector<Node> getVariables() {return this->searchVariables; }

    double localScore(const Node& x, const std::vector<Node>& z);

    /**
     * @return the list of variable varNames.
     */
    std::vector<std::string> getVariableNames();

    Node getVariable(std::string name);

    double getPenalty() { return this->penalty; }

    void setPenalty(double penalty) { this->penalty = penalty; }

    // DataSet getData() { return DataSet(); }

    bool isVerbose() { return this->verbose; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    int getSampleSize() { return N; }

    // std::vector<DataSet*> getDataSets() { return { NULL }; }

    friend double DGScoreTest(const Rcpp::DataFrame& df,
			      std::string targetName,
			      std::vector<std::string>& regressorNames);

    
};


#endif /* DEGENERATEGAUSSIANSCORE_HPP_ */
