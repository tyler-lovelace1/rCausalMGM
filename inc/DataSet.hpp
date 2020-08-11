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
  std::vector<Variable*> variables;
  std::vector<std::string> variableNames;
  std::unordered_map<std::string, int> name2idx;
  std::unordered_map<Variable*, int> var2idx;
  arma::mat data;
  int maxDiscrete;
  int m, n;

  std::set<std::string> getUnique(const Rcpp::CharacterVector& col);

public:
  DataSet() {}
  DataSet(const int maxDiscrete) { this->maxDiscrete=maxDiscrete; }
  DataSet(const Rcpp::DataFrame& df, const int maxDiscrete);
  DataSet(DataSet& ds);
  DataSet& operator=(DataSet& ds);
  DataSet(DataSet&& ds);
  DataSet& operator=(DataSet&& ds);
  ~DataSet();

  int getNumRows() { return n; }
  int getNumColumns() {return m; }

  void set(int i, int j, double val) { data(i,j) = val; }
  void set(int i, int j, int val) { data(i,j) = (double) val; }

  void addVariable(Variable* v);
  void addVariable(int i, Variable* v);

  Variable* getVariable(int i) { return variables[i]; }
  Variable* getVariable(std::string name) { return variables[name2idx[name]]; }

  std::vector<Variable*> getVariables() { return variables; }


  int getColumn(Variable* v) { return var2idx[v]; }

  int getInt(int row, int col);

  arma::mat getData() { return data; }

  std::vector<std::string> getVariableNames() { return variableNames; }

  friend void DataSetTest(const Rcpp::DataFrame& df, const int maxDiscrete);
  friend std::ostream& operator<<(std::ostream& os, DataSet& ds);

};

#endif /* DATASET_HPP_ */
