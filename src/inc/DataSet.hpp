/*
 * Represents a real-valued variable. The values are doubles, and the default
 * missing value marker for is numeric_limits<double>::quiet_NaN().
 *
 * @author Willie Wheeler 7/99
 * @author Joseph Ramsey modifications 12/00
`* @author Tyler Lovelace  1/20
 */

#ifndef DATASET_HPP_
#define DATASET_HPP_

class DataSet;

#include "armaLapack.hpp"
#include "Node.hpp"
#include "ContinuousVariable.hpp"
#include "DiscreteVariable.hpp"
#include <exception>
#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>
#include <typeinfo>
#include <iostream>

class DataSet {
private:
    std::vector<Node> variables;
    std::vector<std::string> variableNames;
    std::unordered_map<std::string, int> name2idx;
    std::unordered_map<Node, int> var2idx;
    arma::mat data;
    arma::umat missing;
    int maxDiscrete;
    int m, n;

    std::set<std::string> getUnique(const Rcpp::CharacterVector &col);
    // bool checkCensoring(const Rcpp::CharacterVector& col,
    // 			arma::vec& values,
    // 			arma::uvec& censor);

public:
    DataSet() {}
    DataSet(const int maxDiscrete) { this->maxDiscrete=maxDiscrete; }
    DataSet(const Rcpp::DataFrame& df);
    // DataSet(const Rcpp::DataFrame& df, const int maxDiscrete);
    // DataSet(const DataSet& ds, arma::urowvec rows); // subset rows
    DataSet(const DataSet& ds, const arma::urowvec& rows); // subset rows
    
    DataSet(const DataSet& ds) = default;
    DataSet& operator=(const DataSet& ds) = default;
    DataSet(DataSet&& ds) = default;
    DataSet& operator=(DataSet&& ds) = default;
    ~DataSet() = default;

    void dropMissing();
    void npnTransform();

    int getNumRows() { return n; }
    int getNumColumns() { return m; }

    void set(int i, int j, double val) { data(i, j) = val; }
    void set(int i, int j, int val) { data(i, j) = (double) val; }

    void addVariable(Node v);
    void addVariable(int i, Node v);

    void removeVariable(Node v);

    const Node& getVariable(int i) { return variables[i]; }
    const Node& getVariable(std::string name) const { return variables.at(name2idx.at(name)); }
    const std::vector<Node>& getVariables() { return variables; }
    std::vector<Node> getContinuousVariables();
    std::vector<Node> getDiscreteVariables();
    // std::vector<Node> getCensoredVariables();
    // std::vector<Node> copyVariables();
    // std::vector<Node> copyContinuousVariables();
    // std::vector<Node> copyDiscreteVariables();

    bool isMixed();
    bool isContinuous();
    bool isDiscrete();
    // bool isCensored();
    
    // void deleteVariables();

    int getInt(int row, int col);

    std::vector<std::string> getVariableNames() { return variableNames; }

    arma::mat getData() { return data; }

    arma::mat getSubsetData(arma::uvec colIdxs);
    arma::mat getSubsetData(arma::uvec rowIdxs, arma::uvec colIdxs);
    
    arma::mat getContinuousData();
    arma::mat getDiscreteData();
    // arma::mat getCensoredData();

    std::vector<int> getDiscLevels();

    int getColumn(Node v) { return var2idx.at(v); }

    /**
     * Removes a node from the dataset and replaces it with an updated version.
     */
    bool updateNode(const Node& node);

    friend void DataSetTest(const Rcpp::DataFrame &df, const int maxDiscrete);
    friend arma::mat DataSetNPNTest(const Rcpp::DataFrame &df);
    friend std::ostream &operator<<(std::ostream &os, DataSet &ds);

};

#endif /* DATASET_HPP_ */
