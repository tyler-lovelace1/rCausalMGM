#ifndef LOGISTICREGRESSIONRESULT_HPP_
#define LOGISTICREGRESSIONRESULT_HPP_

#include "armaLapack.hpp"
#include <string>
#include <vector>
#include <iostream>

/**
 * Stores the various components of a logisitc regression result so they can be passed
 * around together more easily.
 *
 * @author Joseph Ramsey
 * @author Jack Fiore Conversion to c++ 5/11
 */

class LogisticRegressionResult
{
private:
    std::vector<std::string> sigMarker;
    double chiSq;
    double alpha;
    std::vector<std::string> regressorNames;
    std::string target;
    int ny0;
    int ny1;
    int numRegressors;
    arma::vec coefs;
    arma::vec stdErrs;
    arma::vec probs;
    arma::vec xMeans;
    arma::vec xStdDevs;
    double intercept;
    double logLikelihood;


    /**
     * Constructs a new LogRegrResult.
     *
     * @param ny0           the number of cases with target = 0.
     * @param ny1           the number of cases with target = 1.
     * @param numRegressors the number of regressors
     * @param coefs         the array of regression coefficients.
     * @param stdErrs       the array of std errors of the coefficients.
     * @param probs         the array of P-values for the regression
     *                      coefficients.
     */
public:
    LogisticRegressionResult() {}
    
    LogisticRegressionResult(std::string target, std::vector<std::string> regressorNames,
			     arma::vec xMeans, arma::vec xStdDevs, int numRegressors, int ny0,
			     int ny1, arma::vec coefs, arma::vec stdErrs, arma::vec probs,
			     double intercept, double logLikelihood,
			     std::vector<std::string> sigmMarker, double chiSq, double alpha);

    /**
     * The variables.
     */
    std::vector<std::string> getRegressorNames() { return this->regressorNames; }

    /**
     * The target.
     */
    std::string getTarget() { return this->target; }

    /**
     * The number of data points with target = 0.
     */
    int getNy0() { return this->ny0; }

    /**
     * The number of data points with target = 1.
     */
    int getNy1() { return this->ny1; }

    /**
     * The number of regressors.
     */
    int getNumRegressors() { return this->numRegressors; }

    /**
     * The vector of regression coefficients.
     */
    arma::vec getCoefs() { return this->coefs; }

    /**
     * The array of standard errors for the regression coefficients.
     */
    arma::vec getStdErrs() { return this->stdErrs; }

    /**
     * The array of P-values for the regression coefficients.
     */
    arma::vec getProbs() { return this->probs; }

    /**
     * THe array of means.
     */
    arma::vec getxMeans() { return this->xMeans; }

    /**
     * The array of standard devs.
     */
    arma::vec getxStdDevs() { return this->xStdDevs; }

    double getIntercept() { return this->intercept; }

    /**
     * The log likelihood of the regression
     */
    double getLogLikelihood() { return this->logLikelihood; }

    friend std::ostream& operator<<(std::ostream& os, LogisticRegressionResult& lrr);

};


#endif /* LOGISTICREGRESSIONRESULT_HPP_ */
