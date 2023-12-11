#include "CoxRegressionResult.hpp"
#include <exception>

CoxRegressionResult::CoxRegressionResult(std::vector<std::string> regressorNames, int n,
					 arma::vec b, arma::vec t, arma::vec p, arma::vec se,
					 arma::vec resid, double alpha, double loglikelihood) {
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

    this->regressorNames = regressorNames;

    this->n = n;
    this->b = b;
    this->t = t;
    this->p = p;
    this->se = se;
    this->alpha = alpha;
    this->loglikelihood = loglikelihood;
    this->resid = resid;
}

std::ostream& operator<<(std::ostream& os, const CoxRegressionResult& crr) {
    os << "Cox Regression Result Summary:";
    os << "\n";
    os << "List of Regressor Names: ";
    for (int i = 0; i < crr.regressorNames.size(); i++) {
	os << crr.regressorNames[i];
	os << "\t";
    }
    os << "\n";
    os << "Number of Data Points: ";
    os << crr.n;
    os << "\n";
    os << "Cox Regression coefficients: ";
    for (int i = 0; i < crr.b.size(); i++) {
	os << crr.b[i];
	os << "\t";
    }
    os << "\n";
    os << "Hazard Ratios: ";
    for (int i = 0; i < crr.b.size(); i++) {
	os << std::exp(crr.b[i]);
	os << "\t";
    }
    os << "\n";
    os << "t-statistics: ";
    for (int i = 0; i < crr.t.size(); i++) {
	os << crr.t[i];
	os << "\t";
    }
    os << "\n";
    os << "P-Values: ";
    for (int i = 0; i < crr.p.size(); i++) {
	os << crr.p[i];
	os << "\t";
    }
    os << "\n";
    os << "Standard Errors: ";
    for (int i = 0; i < crr.se.size(); i++) {
	os << crr.se[i];
	os << "\t";
    }
    os << "\n";
    os << "Alpha value: ";
    os << crr.alpha;
    os << "\n";
    os << "Log-likelihood: ";
    os << crr.loglikelihood;
    os << "\n";
    os << "\n";
    return os;
}
