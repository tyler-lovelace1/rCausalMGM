#ifndef REGRESSIONRESULT_HPP_
#define REGRESSIONRESULT_HPP_

#include "armaLapack.hpp"
#include <string>
#include <vector>
#include <iostream>


/**
 * Stores the various components of a regression result so they can be passed
 * around together more easily.
 *
 * @author Joseph Ramsey
 * @author Jack Fiore Conversion to c++ 2/30
 */

class RegressionResult
{
  private:
    /**
     * True iff this model assumes a zero intercept.
     *
     * @see #RegressionResult
     */
    bool zeroInterceptAssumed;

    /**
     * Regressor names.
     *
     * @see #RegressionResult
     */
    std::vector<std::string> regressorNames;

    /**
     * The number of data points.
     *
     * @see #RegressionResult
     */
    int n;

    /**
     * The array of regression coefficients.
     *
     * @see #RegressionResult
     */
    arma::vec b;

    /**
     * The array of t-statistics for the regression coefficients.
     *
     * @see #RegressionResult
     */
    arma::vec t;

    /**
     * The array of p-values for the regression coefficients.
     *
     * @see #RegressionResult
     */
    arma::vec p;

    /**
     * Standard errors of the coefficients
     *
     * @see #RegressionResult
     */
    arma::vec se;

    /**
     * R square value.
     *
     * @see #RegressionResult
     */
    double r2;

    /**
     * Residual sums of squares.
     *
     * @see #RegressionResult
     */
    double rss;

    /**
     * Alpha value to determine significance.
     *
     * @see #RegressionResult
     */
    double alpha;

    /**
     * The predicted values.
     */
    arma::vec yHat;

    /**
     * The residuals.
     */
    arma::vec res;

    /**
     * A result for a variety of regression algorithm.
     *
     * @param zeroInterceptAssumed True iff a zero intercept was assumed in
     *                             doing the regression, in which case this
     *                             coefficient is provided; otherwise, not.
     * @param regressorNames       The list of regressor variable names, in
     *                             order.
     * @param n                    The sample size.
     * @param b                    The list of coefficients, in order. If a zero
     *                             intercept was not assumed, this list begins
     *                             with the intercept.
     * @param t                    The list of t-statistics for the
     *                             coefficients, in order. If a zero intercept
     *                             was not assumed, this list begins with the t
     *                             statistic for the intercept.
     * @param p                    The p-values for the coefficients, in order.
     *                             If a zero intercept was not assumed, this
     *                             list begins with the p value for the
     *                             intercept.
     * @param se                   The standard errors for the coefficients, in
     *                             order. If a zero intercept was not assumed,
     *                             this list begins with the standard error of
     *                             the intercept.
     * @param r2                   The R squared statistic for the regression.
     * @param rss                  The residual sum of squares of the
     *                             regression.
     * @param alpha                The alpha value for the regression,
     *                             determining which regressors are taken to be
     */

  public:
    RegressionResult() {}
    RegressionResult(bool zeroInterceptAssumed, std::vector<std::string> regressorNames,
                     int n, arma::vec b, arma::vec t, arma::vec p, arma::vec se, double r2,
                     double rss, double alpha, arma::vec yHat, arma::vec res);
    // ~RegressionResult();

    bool isZeroInterceptAssumed() { return this->zeroInterceptAssumed; }

    double getRSquared() { return this->r2; }

    double getRSS() { return this->rss; }

    int getN() { return this->n; }

    int getNumRegressors() { return this->regressorNames.size(); }

    arma::vec getCoef() { return this->b; }

    arma::vec getT() { return this->t; }

    arma::vec getP() { return this->p; }

    arma::vec getSe() { return this->se; }

    std::vector<std::string> getRegressorNames() { return this->regressorNames; }

    double getPredictedValue(arma::vec x);

    arma::vec getYHat() { return this->yHat; }

    arma::vec getResiduals() { return this->res; }

    friend std::ostream& operator<<(std::ostream& os, RegressionResult& rr);

};

#endif /* REGRESSIONRESULT_HPP_ */
