#ifndef LOGISTICREGRESSIONRESULT_HPP_
#define LOGISTICREGRESSIONRESULT_HPP_

#include "armaLapack.hpp"
#include <string>
#include <vector>

/**
 * Stores the various components of a logistic regression result in a single class so
 * that they can be passed together as an argument or return value.
 *
 * @author Frank Wimberly
 * @author Jack Fiore Conversion to c++ 2/30
 */

class LogisticRegressionResult
{
  private:
    /**
     * String representation of the result
     */
    std::string result;

    /**
     * The variables.
     */
    std::vector<std::string> variableNames;


    /**
     * The target.
     */
    std::string target;

    /**
     * The number of data points with target = 0.
     */
    int ny0;

    /**
     * The number of data points with target = 1.
     */
    int ny1;


    /**
     * The number of regressors.
     */
    int numRegressors;

    /**
     * The array of regression coefficients.
     */
    arma::vec coefs;

    /**
     * The array of standard errors for the regression coefficients.
     */
    arma::vec  stdErrs;

    /**
     * The array of coefP-values for the regression coefficients.
     */
    arma::vec probs;


    /**
     * THe array of means.
     */
    arma::vec xMeans;


    /**
     * The array of standard devs.
     */
    arma::vec xStdDevs;


    double intercept;

    /**
     * The log likelyhood of the regression
     */
    double logLikelihood;

    /**
     * Constructs a new LinRegrResult.
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
    LogisticRegressionResult(std::string target, std::vector<std::string> variableNames,
                    arma::vec xMeans, arma::vec xStdDevs, int numRegressors, int ny0,
                    int ny1, arma::vec coefs, arma::vec stdErrs, arma::vec probs,
                    double intercept, std::string result, double logLikelihood);

    std::string getTarget() {return this->target;}

    double getIntercept() {return this->intercept;}

    /**
     * @return the number of regressors.
     */
    int getNumRegressors() {return numRegressors;}

    /**
     * @return the number of cases with target = 0.
     */
    int getNy0() {return ny0;}

    /**
     * @return the number of cases with target = 1.
     */
    int getNy1() {return ny1;}

    /**
     * @return the total number of cases.
     */
    int getnCases() {return ny0 + ny1;}

    /**
     * @return the array of strings containing the variable names.
     */
    std::vector<std::string> getVariableNames() {return variableNames;}

    /**
     * @return the array of regression coeffients.
     */
    arma::vec getCoefs() {return coefs;}

    /**
     * @return the array of coefT-statistics for the regression coefficients.
     */
     arma::vec getStdErrs() {return stdErrs;}

    /**
     * @return the array of coefP-values for the regression coefficients.
     */
    arma::vec getProbs() {return probs;}

    std::vector<std::string> getVarNames() {return variableNames;}

    public std::string toString() {return result;}

    arma::vec getxMeans() {return xMeans;}

    arma::vec getxStdDevs() {return xStdDevs;}

    /**
     * @return -2LogLiklihood
     */
    double getLogLikelihood() {return logLikelihood;}

}

#endif /* LOGISTICREGRESSIONRESULT_HPP_ */
