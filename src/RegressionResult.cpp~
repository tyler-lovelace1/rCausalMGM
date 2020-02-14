#include "RegressionResult.hpp"
#include <exception>

RegressionResult::RegressionResult(bool zeroInterceptAssumed, std::vector<std::string> regressorNames, int n, arma::vec b, arma::vec t, arma::vec p, arma::vec se, double r2, double rss, double alpha, arma::vec yHat, arma::vec res) {
  if (b.empty()) {
    throw std::invalid_argument("The array of Regression Coefficients cannot be null");
  }

  if (t.empty()) {
    throw std::invalid_argument("The array of T statistics cannot be null");
  }

  if (p.empty()) {
    throw std::invalid_argument("The array of p-values cannot be null");
  }

  if (se.empty()) {
    throw std::invalid_argument("The array of standard errors cannot be null");
  }

  this->zeroInterceptAssumed = zeroInterceptAssumed;

  // Need to set this one before calling getNumRegressors.
  this->regressorNames = regressorNames;

  this->n = n;
  this->b = b;
  this->t = t;
  this->p = p;
  this->se = se;
  this->r2 = r2;
  this->alpha = alpha;
  this->rss = rss;

  this->yHat = yHat;
  this->res = res;
}

double RegressionResult::getPredictedValue(arma::vec x) {
    double yHat = 0.0;

    int offset = zeroInterceptAssumed ? 0 : 1;


    for (int i = 0; i < getNumRegressors(); i++) {
        yHat += b[i + offset] * x[i];
    }

    if (!zeroInterceptAssumed) {
        yHat += b[0]; // error term.
    }

    return yHat;
}
