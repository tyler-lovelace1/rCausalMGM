#ifndef LOGISTICREGRESSION_HPP_
#define LOGISTICREGRESSION_HPP_

#include "armaLapack.hpp"
#include "DataSet.hpp"
#include "Variable.hpp"
#include "LogisticRegressionResult.hpp"

class LogisticRegression
{
  private:
    DataSet data;
    double alpha;
    // double dataCols[][]; ?
    arma::mat dataCols;
    arma::uvec rows;

  public:
    LogisticRegression() {}
    LogisticRegression(DataSet& data);
    LogisticRegression(LogisticRegression& lr);
    LogisticRegression(LogisticRegression&& lr);

    LogisticRegression& operator=(LogisticRegression& lr);
    LogisticRegression& operator=(LogisticRegression&& lr);

    LogisticRegressionResult regress(DiscreteVariable* x, std::vector<Variable*>& regressors); // double regressors[][]

  LogisticRegressionResult regress(DiscreteVariable* x, std::vector<Variable*>& regressors, arma::uvec _rows); // double regressors[][]

    LogisticRegressionResult regress(arma::uvec& target, std::string targetName, arma::mat& regressors, std::vector<std::string>& regressorNames); //double regressors[][]

    double getAlpha() { return this->alpha; }

    void setAlpha(double alpha) { this->alpha = alpha; }

    arma::uvec getRows() { return this->rows; }

    void setRows(arma::uvec rows) { this->rows = rows; }

    double norm(double z);

    friend void LinearRegressionTest(const Rcpp::DataFrame& df);
};

#endif /* LOGISTICREGRESSION_HPP_ */
