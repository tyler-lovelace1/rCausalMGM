#ifndef COXIRLSREGRESSION_HPP_
#define COXIRLSREGRESSION_HPP_

#include "armaLapack.hpp"
#include "DataSet.hpp"
#include "CoxRegressionResult.hpp"
// #include <list>

class CoxIRLSRegression {
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
    std::vector<Node> variables;

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
    CoxIRLSRegression() {}
    CoxIRLSRegression(DataSet& data);
    CoxIRLSRegression(const CoxIRLSRegression& other) = default;
    CoxIRLSRegression(CoxIRLSRegression&& other) = default;

    CoxIRLSRegression& operator=(const CoxIRLSRegression& other) = default;
    CoxIRLSRegression& operator=(CoxIRLSRegression&& other) = default;

    ~CoxIRLSRegression() = default;

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
    CoxRegressionResult regress(const Node& target,
				std::vector<Node>& regressors);

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
    CoxRegressionResult regress(const Node& target,
    				std::vector<Node>& regressors,
    				arma::uvec& _rows);

    double loss(arma::vec& beta, arma::mat& X, Node target);
    
    double gradHess(arma::vec& eta, arma::vec& grad, arma::vec& diagHess,
		    Node target);

    void infoMat(arma::vec& beta, arma::mat& hess, arma::mat& X, Node target);

    arma::uvec getRows() { return rows; }

    void setRows(arma::uvec rows) { this->rows = rows; }

    std::vector<Node> getVariables() { return variables; }

    // arma::vec getResidualsWithoutFirstRegressor() {return res2;}

    friend void CoxIRLSRegressionTest(const Rcpp::DataFrame& df, std::string targetName,
				      std::vector<std::string>& regressorNames, int repetitions);

};
#endif /* COXIRLSREGRESSION_HPP_ */
