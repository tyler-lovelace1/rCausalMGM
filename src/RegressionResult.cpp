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

std::ostream& operator<<(std::ostream& os, RegressionResult& rr) {
  os << "Regression Result Summary:";
  os << "\n";
  os << "Does this model assume a zero intercept? ";
  os << rr.zeroInterceptAssumed;
  os << "\n";
  os << "List of Regressor Names: ";
  for (int i = 0; i < rr.regressorNames.size(); i++) {
    os << rr.regressorNames[i];
    os << "\t";
  }
  os << "\n";
  os << "Number of Data Points: ";
  os << rr.n;
  os << "\n";
  os << "Regression coefficients: ";
  for (int i = 0; i < rr.b.size(); i++) {
    os << rr.b[i];
    os << "\t";
  }
  os << "\n";
  os << "t-statistics: ";
  for (int i = 0; i < rr.t.size(); i++) {
    os << rr.t[i];
    os << "\t";
  }
  os << "\n";
  os << "P-Values: ";
  for (int i = 0; i < rr.p.size(); i++) {
    os << rr.p[i];
    os << "\t";
  }
  os << "\n";
  os << "Standard Errors: ";
  for (int i = 0; i < rr.se.size(); i++) {
    os << rr.se[i];
    os << "\t";
  }
  os << "\n";
  os << "R-Squared: ";
  os << rr.r2;
  os << "\n";
  os << "Alpha value: ";
  os << rr.alpha;
  os << "\n";
  os << "Residual Sum of Squares: ";
  os << rr.rss;
  os << "\n";
  // os << "Y-Hat values: ";
  // for (int i = 0; i < rr.yHat.size(); i++) {
  //   os << rr.yHat[i];
  //   os << "\n";
  // }
  // os << "\n";
  // os << "Residuals: ";
  // for (int i = 0; i < rr.res.size(); i++) {
  //   os << rr.res[i];
  //   os << "\n";
  // }
  os << "\n";
  return os;
}
