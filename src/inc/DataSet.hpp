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
#include "CensoredVariable.hpp"
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
    bool checkCensoring(const Rcpp::CharacterVector& col,
			arma::vec& values,
			arma::uvec& censor);

public:
    DataSet() {}
    DataSet(const int maxDiscrete) { this->maxDiscrete=maxDiscrete; }
    DataSet(const Rcpp::DataFrame& df);
    DataSet(const Rcpp::DataFrame& df, const int maxDiscrete);
    DataSet(const DataSet& ds, const arma::urowvec& rows); // subset rows

    // template<int RTYPE>
    // DataSet(const Rcpp::DataFrame& df) {
    // 	this->maxDiscrete = maxDiscrete;
    // 	const Rcpp::CharacterVector names = df.names();
    // 	this->m = names.length();
    // 	this->n = df.nrows();
    // 	// this->data.set_size(this->n, this->m);
    // 	this->data = arma::mat(n,m,arma::fill::zeros);
    // 	this->name2idx = std::unordered_map<std::string, int>();
    // 	this->var2idx = std::unordered_map<Node, int>();
    // 	int numUnique;
    // 	std::string val, curName;

    // 	for (int i = 0; i < m; i++) {
	
    // 	    curName = (std::string)names[i];

    // 	    Rcpp::Vector<RTYPE> colVec = df[i];
	
    // 	    if (colVec.attr("class") == "Surv") {
    // 		arma::vec values = Rcpp::as<arma::vec>(colVec[0]);
    // 		arma::uvec censor = Rcpp::as<arma::uvec>(colVec[1]);
    // 		variables.push_back(Node(new CensoredVariable(curName)));
    // 		variables[i].setCensor(values, censor);
    // 		data.col(i) = values;
    // 	    } else if (colVec.attr("class") == "factor") {
    // 		std::vector<std::string> levels = Rcpp::as<std::vector<std::string>>(colVec.attr("levels"));
    // 		arma::vec values = Rcpp::as<arma::vec>(colVec);
	    
    // 		variables.push_back(Node(new DiscreteVariable(curName, levels)));
    // 		data.col(i) = values;
    // 	    } else {
    // 		arma::vec values = Rcpp::as<arma::vec>(colVec);
    // 		variables.push_back(Node(new ContinuousVariable(curName)));
    // 		data.col(i) = values;
    // 	    }
    // 	}
    // }
    
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

    const Node& getVariable(int i) { return variables[i]; }
    const Node& getVariable(std::string name) const { return variables.at(name2idx.at(name)); }
    const std::vector<Node>& getVariables() { return variables; }
    std::vector<Node> getContinuousVariables();
    std::vector<Node> getDiscreteVariables();
    std::vector<Node> getCensoredVariables();
    // std::vector<Node> copyVariables();
    // std::vector<Node> copyContinuousVariables();
    // std::vector<Node> copyDiscreteVariables();

    bool isMixed();
    bool isContinuous();
    bool isDiscrete();
    bool isCensored();
    
    // void deleteVariables();

    int getInt(int row, int col);

    std::vector<std::string> getVariableNames() { return variableNames; }

    arma::mat getData() { return data; }
    
    arma::mat getContinuousData();
    arma::mat getDiscreteData();
    arma::mat getCensoredData();

    std::vector<int> getDiscLevels();

    int getColumn(Node v) { return var2idx.at(v); }

    friend void DataSetTest(const Rcpp::DataFrame &df, const int maxDiscrete);
    friend std::ostream &operator<<(std::ostream &os, DataSet &ds);

};

#endif /* DATASET_HPP_ */
