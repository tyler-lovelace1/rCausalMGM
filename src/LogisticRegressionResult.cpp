#include "LogisticRegressionResult.hpp"
#include <exception>

LogisticRegressionResult::LogisticRegressionResult(std::string target, std::vector<std::string> regressorNames,
                                                   arma::vec xMeans, arma::vec xStdDevs, int numRegressors, int ny0,
                                                   int ny1, arma::vec coefs, arma::vec stdErrs, arma::vec probs,
                                                   double intercept, double logLikelihood, std::vector<std::string> sigmMarker,
                                                   double chiSq, double alpha)
{

    if (regressorNames.size() != numRegressors)
    {
        throw std::invalid_argument("Number of Regressor Names in vector does not equal numRegressors");
    }

    if (coefs.size() != numRegressors + 1)
    {
        throw std::invalid_argument("Invalid size for coefs");
    }

    if (stdErrs.size() != numRegressors + 1)
    {
        throw std::invalid_argument("Invalid size for stdErrs");
    }

    if (probs.size() != numRegressors + 1)
    {
        throw std::invalid_argument("Invalid size for probs");
    }

    if (xMeans.size() != numRegressors + 1)
    {
        throw std::invalid_argument("Invalid size for xMeans");
    }

    if (xStdDevs.size() != numRegressors + 1)
    {
        throw std::invalid_argument("Invalid size for xxStdDevs");
    }
    if (target.empty())
    {
        throw std::invalid_argument("invalid target (empty)");
    }

    this->intercept = intercept;
    this->target = target;
    this->xMeans = xMeans;
    this->xStdDevs = xStdDevs;
    this->regressorNames = regressorNames;
    this->numRegressors = numRegressors;
    this->ny0 = ny0;
    this->ny1 = ny1;
    this->coefs = coefs;
    this->stdErrs = stdErrs;
    this->probs = probs;
    this->logLikelihood = logLikelihood;
    this->sigMarker = sigmMarker;
    this->chiSq = chiSq;
    this->alpha = alpha;
}

std::ostream &operator<<(std::ostream &os, LogisticRegressionResult &lrr)
{
    os << "Regression Result Summary:";

    os << "\n";
    os << "Regressor Names ";
    std::vector<std::string> rnames = lrr.getRegressorNames();
    for (int i = 0; i < rnames.size(); i++)
    {
        os << rnames[i];
        os << "\t";
    }

    os << "\n";
    os << "Target: ";
    std::string targ = lrr.getTarget();
    os << targ;

    os << "\n";
    os << "The number of data points with target = 0: ";
    int ny0 = lrr.getNy0();
    os << ny0;

    os << "\n";
    os << "The number of data points with target = 1: ";
    int ny1 = lrr.getNy1();
    os << ny1;

    os << "\n";
    os << "The number of Regressors: ";
    int numR = lrr.getNumRegressors();
    os << numR;

    os << "\n";
    os << "Regression coefficients: ";
    arma::vec regCoefs = lrr.getCoefs();
    for (int j = 0; j < regCoefs.size(); j++)
    {
        os << regCoefs[j];
        os << "\t";
    }

    os << "\n";
    os << "Standard Errors: ";
    arma::vec stdErrs = lrr.getStdErrs();
    for (int l = 0; l < stdErrs.size(); l++)
    {
        os << stdErrs[l];
        os << "\t";
    }

    os << "\n";
    os << "P-Values: ";
    arma::vec pVals = lrr.getProbs();
    for (int m = 0; m < pVals.size(); m++)
    {
        os << pVals[m];
        os << "\t";
    }

    os << "\n";
    os << "X Means: ";
    arma::vec xMeans = lrr.getxMeans();
    for (int n = 0; n < xMeans.size(); n++)
    {
        os << xMeans[n];
        os << "\t";
    }

    os << "\n";
    os << "X Std Devs: ";
    arma::vec xStdDevs = lrr.getxStdDevs();
    for (int o = 0; o < xStdDevs.size(); o++)
    {
        os << xStdDevs[o];
        os << "\t";
    }

    os << "\n";
    os << "Intercept: ";
    double intercept = lrr.getIntercept();
    os << intercept;

    os << "\n";
    os << "Log likelihood of the Regression: ";
    double ll = lrr.getLogLikelihood();
    os << ll;

    os << "\n";
    return os;
}
