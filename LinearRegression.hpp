// make a test function that takes in a data frame into the dataset class
// alll the distribution stuff we need is in the boost library   t distribution and chi squared
// p value and regression coefficients must be checked for correctness

#ifndef LINEARREGRESSION_HPP_
#define LINEARREGRESSION_HPP_

#include "armaLapack.hpp"
#include "DataSet.hpp"
#include "Variable.hpp"
#include "RegressionResult.hpp"
#include <list>

class LinearRegression
{
  private:
    /**
     * The data set.
     */
    DataSet data;

    /**
     * The variables.
     */
    std::vector<Variable*> variables;

    /**
     * The significance level for determining which regressors are significant
     * based on their p values.
     */
    double alpha = 0.05;

    /**
     * The graph of significant regressors into the target.
     */

    arma::uvec rows;
    arma::vec res2;

  public:
    //============================CONSTRUCTORS==========================//

    LinearRegression(DataSet data);
    // LinearRegression(arma::mat data, std::vector<Variable*>&  variables);

    //===========================PUBLIC METHODS========================//

    /**
     * Sets the alpha level for deciding which regressors are significant
     * based on their p values.
     */
    void setAlpha(double alpha) {this->alpha = alpha;}

    /**
     * Regresses the target on the given regressors.
     *
     * @param target     The target variable.
     * @param regressors The regressor variables.
     * @return The regression plane, specifying for each regressors its
     * coefficeint, se, t, and p values, and specifying the same for the
     * constant.
     */
     RegressionResult* regress(Variable* target, std::vector<Variable*>& regressors);

     // static RegressionResult regress(double[] target, double[][] regressors) = 0;

     double rss(arma::mat x, arma::vec y, arma::vec b); // check if only if 1 dimension   is a vector

     double tss(arma::vec y);

     arma::uvec getRows() {return rows;}

     void setRows(arma::uvec rows) {this->rows = rows;}

     arma::vec getResidualsWithoutFirstRegressor() {return res2;}

     friend void LinearRegressionTest(const Rcpp::DataFrame& df);

};
#endif /* LINEARREGRESSION_HPP_ */
