#ifndef LINEARREGRESSION_HPP_
#define LINEARREGRESSION_HPP_

#include "armaLapack.hpp"
#include "DataSet.hpp"
#include "Node.hpp"
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
    // arma::vec res2;

    std::map<std::pair<Node,Node>, arma::vec> WZmap;

public:
    //============================CONSTRUCTORS==========================//
    LinearRegression() {}
    LinearRegression(DataSet& data);
    LinearRegression(const LinearRegression& lr) = default;
    LinearRegression(LinearRegression&& lr) = default;

    LinearRegression& operator=(const LinearRegression& lr) = default;
    LinearRegression& operator=(LinearRegression&& lr) = default;

    ~LinearRegression() = default;

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
    RegressionResult regress(const Node& target, std::vector<Node>& regressors);

    /**
     * Regresses the target on the given regressors.
     *
     * @param target     The target variable.
     * @param regressors The regressor variables.
     * @param _rows The samples to be used in regression.
     * @return The regression plane, specifying for each regressors its
     * coefficeint, se, t, and p values, and specifying the same for the
     * constant.
     */
    RegressionResult regress(const Node& target,
			     std::vector<Node>& regressors,
			     arma::uvec& _rows);

    // static RegressionResult regress(double[] target, double[][] regressors) = 0;

    double rss(const arma::mat& x, const arma::vec& y, const arma::vec& b); // check if only if 1 dimension   is a vector

    double tss(const arma::vec& y);

    arma::uvec getRows() {return rows;}

    void setRows(arma::uvec rows) {this->rows = rows;}

    void setWZmap(std::map<std::pair<Node,Node>, arma::vec>& WZmap) { this->WZmap = WZmap; }

    arma::vec getWZ(Node coxNode, Node target);

    // arma::vec getResidualsWithoutFirstRegressor() {return res2;}

    friend void LinearRegressionTest(const Rcpp::DataFrame& df);

};
#endif /* LINEARREGRESSION_HPP_ */
