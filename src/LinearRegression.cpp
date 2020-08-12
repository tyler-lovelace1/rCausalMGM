#include "LinearRegression.hpp"
#include "RegressionResult.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>

LinearRegression::LinearRegression(DataSet& data){
  this->data = data;
  this->variables = data.getVariables();
  this->rows = arma::uvec(data.getNumRows());
  for (int i = 0; i < data.getNumRows(); i++) rows[i] = i;
}

LinearRegression::LinearRegression(LinearRegression& lr) {
  this->data = lr.data;
  this->variables = this->data.getVariables();
  this->rows = arma::uvec(this->data.getNumRows());
  for (int i = 0; i < data.getNumRows(); i++) rows[i] = i;
}

LinearRegression::LinearRegression(LinearRegression&& lr) {
  this->data = lr.data;
  this->variables = this->data.getVariables();
  this->rows = arma::uvec(this->data.getNumRows());
  for (int i = 0; i < data.getNumRows(); i++) rows[i] = i;
}

LinearRegression& LinearRegression::operator=(LinearRegression& lr) {
  this->data = lr.data;
  this->variables = this->data.getVariables();
  this->rows = arma::uvec(this->data.getNumRows());
  for (int i = 0; i < data.getNumRows(); i++) rows[i] = i;
}

LinearRegression& LinearRegression::operator=(LinearRegression&& lr) {
  this->data = lr.data;
  this->variables = this->data.getVariables();
  this->rows = arma::uvec(this->data.getNumRows());
  for (int i = 0; i < data.getNumRows(); i++) rows[i] = i;
}

RegressionResult* LinearRegression::regress(Variable* target, std::vector<Variable*>& regressors){
  int n = rows.size();
  int k = regressors.size() + 1;

  int target_ = data.getColumn(target);

  arma::uvec regressors_ = arma::uvec(regressors.size());

  for (int i = 0; i < regressors.size(); i++) {
    regressors_[i] = data.getColumn(regressors[i]);
  }

  if (target_ == -1) {
      Rcpp::Rcout <<"\n";
  }

  arma::uvec target_Vec(1);
  target_Vec.fill(target_);
  arma::mat y = data.getData().submat(rows, target_Vec);

  arma::mat xSub = data.getData().submat(rows, regressors_);

  arma::mat x;


  if (regressors.size() > 0) {
      x = arma::mat(xSub.n_rows, xSub.n_cols + 1);

      for (arma::uword i = 0; i < x.n_rows; i++) {
          for (arma::uword j = 0; j < x.n_cols; j++) {
              if (j == 0) {
                  x(i, j) = 1;
              } else {
                  x(i, j) = xSub(i, j - 1);
              }
          }
      }
  } else {
      x = arma::mat(xSub.n_rows, xSub.n_cols);

      for (arma::uword i = 0; i < x.n_rows; i++) {
          for (arma::uword j = 0; j < x.n_cols; j++) {
              x(i, j) = xSub(i, j);
          }
      }
  }


  arma::mat xT = x.t();
  arma::mat xTx = xT * x;
  arma::mat xTxInv = xTx.i();
  arma::mat xTy = xT * y;
  arma::mat b = xTxInv * xTy;

  arma::mat yHat = x * b;
  arma::mat res = y - yHat;

  arma::vec yHat_ = yHat.col(0);
  arma::vec res_ = res.col(0);

  arma::mat b2 = arma::mat(b);
  arma::mat yHat2 = x * b2;

  arma::mat res2 = y - yHat2;
  this->res2 = res2.col(0);

  double rss_ = LinearRegression::rss(x, y, b);
  double se = std::sqrt(rss_ / (n - k));
  double tss_ = LinearRegression::tss(y);
  double r2 = 1.0 - (rss_ / tss_);


  arma::vec sqErr = arma::vec(x.n_cols);
  arma::vec t = arma::vec(x.n_cols);
  arma::vec p = arma::vec(x.n_cols);

  boost::math::students_t dist(n-k);
  for (arma::uword i = 0; i < x.n_cols; i++) {
      double s_ = se * se * xTxInv(i, i);
      double se_ = std::sqrt(s_);
      double t_ = b(i, 0) / se_;
      double p_ = 2 * (1.0 - boost::math::cdf(dist, std::abs(t_)));

      sqErr[i] = se_;
      t[i] = t_;
      p[i] =  p_;
  }

  std::vector<std::string> vNames(regressors.size());

  for (int i = 0; i < regressors.size(); i++) {
    Rcpp::Rcout << regressors[i]->getName() << "\n";
    vNames[i] = regressors[i]->getName(); // getName Function may not be implemented
  }

  // arma::vec bArray = b.columns() == 0 ? new double[0] : b.getColumn(0).toArray(); // double check,
                                    //dealing with case where we dont give it anything to regess on

  return new RegressionResult(regressors.size() == 0, vNames, n, b, t, p, sqErr, r2, rss_, alpha, yHat_, res_); // MUST CONVERT B INTO A VECTOR
}

double LinearRegression::rss(arma::mat x, arma::vec y, arma::vec b) {
    double rss = 0.0;

    for (arma::uword i = 0; i < x.n_rows; i++) {
        double yH = 0.0;

        for (arma::uword j = 0; j < x.n_cols; j++) {
            yH += b(j) * x(i, j);
        }

        double d = y(i) - yH;

        rss += d * d;
    }

    return rss;
}

double LinearRegression::tss(arma::vec y) {
    double mean = 0.0;

    for (arma::uword i = 0; i < y.n_rows; i++) {
        mean += y(i);
    }

    mean /= (double) (y.n_rows);

    double ssm = 0.0;

    for (arma::uword i = 0; i < y.n_rows; i++) {
        double d = mean - y(i);
        ssm += d * d;
    }
    return ssm;
}

// [[Rcpp::export]]
void LinearRegressionTest(const Rcpp::DataFrame& df) {
  DataSet data = DataSet(df, 5);
  LinearRegression reg(data);

  for (int i = 0; i < data.getNumColumns(); i++) {
    Rcpp::Rcout << "-----START----- \n";
    Variable* target = data.getVariable(i);
    Rcpp::Rcout << target->getName() << std::endl;
    std::vector<Variable*> regressors(data.getVariables());
    regressors.erase(regressors.begin() + i);
    RegressionResult* result = reg.regress(target, regressors);
    Rcpp::Rcout << *result;
    Rcpp::Rcout << "-----END----- \n";
  }
}
