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
#include <map>
#include <vector>
#include <typeinfo>

class DataSet {
private:
  std::vector<Variable*> variables;
  std::vector<std::string> variableNames;
  std::map<std::string, int> name2idx;
  arma::mat data;
  int maxDiscrete;
  int m, n;

  std::set<std::string> getUnique(const Rcpp::CharacterVector& col);

public:
  DataSet(const int& maxDiscrete) { this->maxDiscrete=maxDiscrete; }
  DataSet(const Rcpp::DataFrame& df, const int& maxDiscrete);
  ~DataSet();

  friend Rcpp::RObject rCausalMGMData(const Rcpp::DataFrame& df, const int& maxDiscrete);
  
};

#endif /* DATASET_HPP_ */
