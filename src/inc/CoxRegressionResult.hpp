#ifndef COXREGRESSIONRESULT_HPP_
#define COXREGRESSIONRESULT_HPP_

#include "armaLapack.hpp"
#include <string>
#include <vector>
#include <iostream>


/**
 * Stores the various components of a Cox Regression result so they can be passed
 * around together more easily.
 *
 * @author Tyler Lovelace
 */

class CoxRegressionResult
{
  private:

    /**
     * Regressor names.
     *
     * @see #CoxRegressionResult
     */
    std::vector<std::string> regressorNames;

    /**
     * The number of data points.
     *
     * @see #CoxRegressionResult
     */
    int n;

    /**
     * The array of Cox regression coefficients.
     *
     * @see #CoxRegressionResult
     */
    arma::vec b;

    /**
     * The array of t-statistics for the Cox regression coefficients.
     *
     * @see #CoxRegressionResult
     */
    arma::vec t;

    /**
     * The array of p-values for the Cox regression coefficients.
     *
     * @see #CoxRegressionResult
     */
    arma::vec p;

    /**
     * Standard errors of the coefficients
     *
     * @see #CoxRegressionResult
     */
    arma::vec se;

    /**
     * Residuals for Cox IRLS model
     *
     * @see #CoxRegressionResult
     */
    arma::vec resid;

    /**
     * Alpha value to determine significance.
     *
     * @see #CoxRegressionResult
     */
    double alpha;
    
    /**
     * Log-likelihood of the Cox regression.
     *
     * @see #CoxRegressionResult
     */
    double loglikelihood;

    /**
     * A result for a variety of coxregression algorithm.
     *
     * @param regressorNames       The list of regressor variable names, in
     *                             order.
     * @param n                    The sample size.
     * @param b                    The list of coefficients, in order.
     * @param t                    The list of t-statistics for the
     *                             coefficients, in order.
     * @param p                    The p-values for the coefficients, in order.
     * @param se                   The standard errors for the coefficients, in
     *                             order. 
     * @param resid                The scaled residuals for IRLS Cox regression 
     * @param alpha                The alpha value for the Cox regression,
     *                             determining which regressors are taken to be
     * @param loglikelihood        The log-likelihood for the Cox regression,
     *                             determining goodness of fit of the model
     */

  public:
    CoxRegressionResult() {}
    CoxRegressionResult(std::vector<std::string> regressorNames, int n, arma::vec b,
			arma::vec t, arma::vec p, arma::vec se, arma::vec resid,
			double alpha, double loglikelihood);
    // ~CoxregressionResult();

    double getLoglikelihood() { return this->loglikelihood; }

    int getN() { return this->n; }

    int getNumRegressors() { return this->regressorNames.size(); }

    arma::vec getCoef() { return this->b; }

    arma::vec getT() { return this->t; }

    arma::vec getP() { return this->p; }

    arma::vec getSe() { return this->se; }

    arma::vec getResid() { return this->resid; }

    std::vector<std::string> getRegressorNames() { return this->regressorNames; }

    friend std::ostream& operator<<(std::ostream& os, const CoxRegressionResult& crr);

};

#endif /* COXREGRESSIONRESULT_HPP_ */
