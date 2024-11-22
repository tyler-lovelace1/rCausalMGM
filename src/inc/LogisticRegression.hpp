#ifndef LOGISTICREGRESSION_HPP_
#define LOGISTICREGRESSION_HPP_

#include "armaLapack.hpp"
#include "DataSet.hpp"
#include "Node.hpp"
#include "LogisticRegressionResult.hpp"

class LogisticRegression
{
  private:
    DataSet data;
    double alpha;
    // double dataCols[][]; ?
    arma::mat dataCols;
    arma::uvec rows;

    // std::map<std::pair<Node,Node>, arma::vec> WZmap;

    bool binary(const Node& x) { return x.getNumCategories()==2; }

  public:
    LogisticRegression() {}
    LogisticRegression(DataSet& data);
    LogisticRegression(const LogisticRegression& lr) = default;
    LogisticRegression(LogisticRegression&& lr) = default;

    LogisticRegression& operator=(const LogisticRegression& lr) = default;
    LogisticRegression& operator=(LogisticRegression&& lr) = default;

    ~LogisticRegression() = default;

    LogisticRegressionResult regress(const Node& x, std::vector<Node>& regressors); // double regressors[][]

    LogisticRegressionResult regress(const Node& x,
				     std::vector<Node>& regressors,
				     arma::uvec _rows); // double regressors[][]

    LogisticRegressionResult regress(arma::uvec& target,
				     std::string targetName,
				     arma::mat& regressors,
				     std::vector<std::string>& regressorNames); //double regressors[][]

    double getAlpha() { return this->alpha; }

    void setAlpha(double alpha) { this->alpha = alpha; }

    arma::uvec getRows() { return this->rows; }

    void setRows(arma::uvec rows) { this->rows = rows; }

    // void setWZmap(std::map<std::pair<Node,Node>, arma::vec>& WZmap) { this->WZmap = WZmap; }

    // arma::vec getWZ(Node coxNode, Node target);
    
    double norm(double z);

    friend void LinearRegressionTest(const Rcpp::DataFrame& df);
};

#endif /* LOGISTICREGRESSION_HPP_ */
