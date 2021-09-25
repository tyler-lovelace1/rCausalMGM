/*
 * Represents a right censored variable with a real value and a censoring indicator. 
 * The values are doubles, and the default missing value marker for is 
 * numeric_limits<double>::quiet_NaN(). The censoring indicator is boolean value, 
 * false when the value is censored and true when uncensored.
 *
 * @author Tyler Lovelace conversion to C++ 7/21
 */

#ifndef CENSOREDVARIABLE_HPP_
#define CENSOREDVARIABLE_HPP_

#include "Variable.hpp"
#include <cmath>
#include <limits>
#include <exception>
#include <iostream>

class CensoredVariable : public Variable {
private:
  const double MISSING_VALUE = std::numeric_limits<double>::quiet_NaN();
  const std::string MISSING_VALUE_STRING = "NaN";
  uint n;
  uint unique;
  uint nEvents;
  arma::uvec order;
  arma::uvec censor;
  arma::uvec H;
  arma::uvec CC;

public:
  CensoredVariable(const std::string& name) {
    this->name = name;
    this->type = CENSORED;
  }

  CensoredVariable(const CensoredVariable& var) {
    name = var.name;
    type = var.type;
    n = var.n;
    unique = var.unique;
    nEvents = var.nEvents;
    order = arma::uvec(var.order);
    censor = arma::uvec(var.censor);
    H = arma::uvec(var.H);
    CC = arma::uvec(var.CC);
  }
  
  ~CensoredVariable() {}

  double getMissingValueMarker() { return this->MISSING_VALUE; }

  bool setCensor(arma::vec& values, arma::uvec& censor);
  
  bool setCensor(arma::vec&& values, arma::uvec&& censor);

  arma::uvec getOrder() { return this->order; }

  arma::uword getOrder(int i) { return this->order[i]; }

  arma::uvec getCensor() { return this->censor; }

  arma::uword getCensor(int i) { return this->censor[i]; }

  arma::uvec getH() { return this->H; }

  arma::uword getH(int i) { return this->H[i]; }

  arma::uvec getCC() { return this->CC; }

  arma::uword getCC(int i) { return this->CC[i]; }

  arma::uword getNEvents() { return nEvents; }

  bool isMissingValue(const double& val) { return std::isnan(val); }

  bool isMissingValue(const std::string& val) { return val==MISSING_VALUE_STRING; }

  // bool checkValue(const double& val) { return true; }

  bool checkValue(const std::string& val);
  
  CensoredVariable* like(const std::string& name) { return new CensoredVariable(name); }

  friend void test_censored(const std::string& name, arma::vec& values, arma::uvec& censor);
  friend std::ostream& operator<<(std::ostream& os, CensoredVariable& cv);
};

#endif /* CENSOREDVARIABLE_HPP_ */
