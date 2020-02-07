#include "LogisticRegressionResult.hpp"
#include <exception>
LogisticRegressionResult::LogisticRegressionResult(std::string target,
  std::vector<std::string> variableNames, arma::vec xMeans, arma::vec xStdDevs,
  int numRegressors, int ny0, int ny1, arma::vec coefs, arma::vec stdErrs,
  arma::vec probs, double intercept, std::string result, double logLikelihood) {

  if (variableNames.size() != numRegressors) {
    throw std::invalid_argument("Size of variableNames must be equal to numRegressors");
  }

  if (coefs.n_elem != numRegressors + 1) {
    throw std::invalid_argument("Number of coefficients must be equal to numRegressors + 1");
  }

  if (stdErrs.n_elem != numRegressors + 1) {
    throw std::invalid_argument("Number of stdErrs must be equal to numRegressors + 1");
  }

  if (probs.n_elem != numRegressors + 1) {
    throw std::invalid_argument("Number of probs must be equal to numRegressors + 1");
  }

  if (xMeans.n_elem != numRegressors + 1) {
    throw std::invalid_argument("Number of xMeans must be equal to numRegressors + 1");
  }

  if (xStdDevs.n_elem != numRegressors + 1) {
    throw std::invalid_argument("Number of xStdDevs must be equal to numRegressors + 1");
  }
  if (target.empty()) {
    throw std::invalid_argument("Target string is empty");
  }

  this->intercept = intercept;
  this->target = target;
  this->xMeans = xMeans;
  this->xStdDevs = xStdDevs;
  this->variableNames = variableNames;
  this->numRegressors = numRegressors;
  this->ny0 = ny0;
  this->ny1 = ny1;
  this->coefs = coefs;
  this->stdErrs = stdErrs;
  this->probs = probs;
  this->result = result;
  this->logLikelihood = logLikelihood;
}
