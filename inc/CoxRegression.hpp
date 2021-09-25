#ifndef COXREGRESSION_HPP_
#define COXREGRESSION_HPP_

#include "armaLapack.hpp"
#include "DataSet.hpp"
#include "CensoredVariable.hpp"
#include "CoxRegressionResult.hpp"
// #include <list>

class CoxRegression {
private:
    /**
     * The data set.
     */
    DataSet data;
    
    /**
     * The data matrix.
     */
    arma::mat dataMat;

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


public:
    //============================CONSTRUCTORS==========================//
    CoxRegression() {}
    CoxRegression(DataSet& data);
    CoxRegression(CoxRegression& lr);
    CoxRegression(CoxRegression&& lr);
    // Coxregression(arma::mat data, std::vector<Variable*>&  variables);

    CoxRegression& operator=(CoxRegression& lr);
    CoxRegression& operator=(CoxRegression&& lr);

    //===========================PUBLIC METHODS========================//

    /**
     * Sets the alpha level for deciding which regressors are significant
     * based on their p values.
     */
    void setAlpha(double alpha) {this->alpha = alpha;}

    /**
     * Cox regression of the target given the regressors.
     *
     * @param target     The target variable.
     * @param regressors The regressor variables.
     * @return The regression plane, specifying for each regressors its
     * coefficeint, se, t, and p values, and specifying the same for the
     * constant.
     */
    CoxRegressionResult regress(CensoredVariable* target, std::vector<Variable*>& regressors);

    /**
     * Cox regression of the target given the regressors.
     *
     * @param target     The target variable.
     * @param regressors The regressor variables.
     * @param _rows The samples to be used in regression.
     * @return The regression plane, specifying for each regressors its
     * coefficeint, se, t, and p values, and specifying the same for the
     * constant.
     */
    CoxRegressionResult regress(CensoredVariable* target, std::vector<Variable*>& regressors,
				arma::uvec& _rows);

    double loss(arma::vec&& beta, arma::mat& X, CensoredVariable* target);

    void gradHess(arma::vec& beta, arma::vec& grad, arma::mat& hess,
		  arma::mat& X, CensoredVariable* target);

    arma::uvec getRows() {return rows;}

    void setRows(arma::uvec rows) {this->rows = rows;}

    // arma::vec getResidualsWithoutFirstRegressor() {return res2;}

    friend void CoxRegressionTest(const Rcpp::DataFrame& df, std::string targetName,
				  std::vector<std::string>& regressorNames);

};
#endif /* COXREGRESSION_HPP_ */
