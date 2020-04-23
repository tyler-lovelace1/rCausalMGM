#include "LinearRegression.hpp"
#include "RegressionResult.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>

LinearRegression::LinearRegression(DataSet data){
  // BOth of these are stored as private member types in the dataset.hpp file
  this->data = data;
  this->variables = data.getVariables(); // this is stored as std::vector<Variable*> in dataset.hpp
  arma::uvec rows;
  this->rows = rows;
  for (int i = 0; i < data.getNumRows(); i++) rows[i] = (i); // data.n is a private variable stored in the dataset header file
  // std::cout << "Constructed!";
  Rcpp::Rcout << "Constructed!";
}

RegressionResult* LinearRegression::regress(Variable* target, std::vector<Variable*>& regressors){
  // int n = arma::uvec size(rows);
  Rcpp::Rcout << "Regressing!";
  int n = rows.size();
  int k = regressors.size() + 1;

  //UNFINISHED
  // int _target = variables.indexOf(target);
  int target_ = data.getColumn(target);

  // int[] _regressors = new int[regressors.size()];
  arma::uvec regressors_ = arma::uvec(regressors.size());

  // UNFINISHED
  // for (int i = 0; i < regressors.size(); i++) {
  //     _regressors[i] = variables.indexOf(regressors.get(i));
  //     if (_regressors[i] == -1) {
  //         System.out.println();
  //     }
  // }
  for (int i = 0; i < regressors.size(); i++) {
    regressors_[i] = data.getColumn(regressors[i]);
  }

  if (target_ == -1) {
      Rcpp::Rcout <<"\n";
  }

  arma::uvec target_Vec(1);
  target_Vec.fill(target_);
  // TetradMatrix y = data.getSelection(getRows(), new int[]{_target}).copy();
  arma::mat y = data.getData().submat(rows, target_Vec);

  // TetradMatrix xSub = data.getSelection(getRows(), _regressors);
  arma::mat xSub = data.getData().submat(rows, regressors_);

  // TetradMatrix x;
  arma::mat x;

  // if (regressors.size() > 0) {
  //     x = arma::mat(xSub.n_rows, xSub.n_cols + 1);
  //
  //     for (int i = 0; i < x.rows(); i++) {
  //         for (int j = 0; j < x.columns(); j++) {
  //             if (j == 0) {
  //                 x.set(i, j, 1);
  //             } else {
  //                 x.set(i, j, xSub.get(i, j - 1));
  //             }
  //         }
  //     }
  // } else {
  //     x = arma::mat(xSub.n_rows, xSub.n_cols);
  //
  //     for (int i = 0; i < x.rows(); i++) {
  //         for (int j = 0; j < x.columns(); j++) {
  //             x.set(i, j, xSub.get(i, j));
  //         }
  //     }
  // }

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

  // The indices of elements are specified via the uword type, which is a typedef
  // for an unsigned integer type. When using loops to access elements, it's best
  // to use uword instead of int. For example: for(uword i=0; i<X.n_elem; ++i) { X(i) = ... }

  // TetradMatrix xT = x.transpose();
  // TetradMatrix xTx = xT.times(x);
  // TetradMatrix xTxInv = xTx.inverse();
  // TetradMatrix xTy = xT.times(y);
  // TetradMatrix b = xTxInv.times(xTy);
  // double check that these operations are the correct operations to perform
  arma::mat xT = x.t();
  arma::mat xTx = xT * x;
  arma::mat xTxInv = xTx.i();
  arma::mat xTy = xT * y;
  arma::mat b = xTxInv * xTy;

  // TetradMatrix yHat = x.times(b);
  arma::mat yHat = x * b;  // product of two matrices


  // if (yHat.columns() == 0) yHat = y.like(); double check

  arma::mat res = y - yHat; //  y.copy().assign(yHat, PlusMult.plusMult(-1));

  arma::vec yHat_ = yHat.col(0);
  arma::vec res_ = res.col(0);

  arma::mat b2 = arma::mat(b); //double check
  arma::mat yHat2 = x * b2;

  // if (yHat.columns() == 0) yHat2 = y.like(); double check

  arma::mat res2 = y - yHat2; //  y.copy().assign(yHat, PlusMult.plusMult(-1));
  this->res2 = res2.col(0);

  double rss_ = LinearRegression::rss(x, y, b);
  double se = std::sqrt(rss_ / (n - k));
  double tss_ = LinearRegression::tss(y);
  double r2 = 1.0 - (rss_ / tss_);


  arma::vec sqErr = arma::vec(x.n_cols);
  arma::vec t = arma::vec(x.n_cols);
  arma::vec p = arma::vec(x.n_cols);

  boost::math::students_t dist(n-k);
  for (int i = 0; i < x.n_cols; i++) {
      double s_ = se * se * xTxInv(i, i);
      double se_ = std::sqrt(s_);
      double t_ = b(i, 0) / se_;
      double p_ = 2 * (1.0 - boost::math::cdf(dist, std::abs(t_))); //must use boost to get cdf for

      sqErr[i] = se_;
      t[i] = t_;
      p[i] =  p_;
  }

  // this.graph = createOutputGraph(target.getName(), x, regressors, p);

  // String[] vNames = new String[regressors.size()];
  std::vector<std::string> vNames;


  // for (int i = 0; i < regressors.size(); i++) {
  //     vNames[i] = regressors.get(i).getName(); // getName Function may not be implemented
  // }

  for (int i = 0; i < regressors.size(); i++) {
      vNames[i] = regressors[i]->getName(); // getName Function may not be implemented
  }

  // double[] bArray = b.columns() == 0 ? new double[0] : b.getColumn(0).toArray();
  // double[] tArray = t.toArray();
  // double[] pArray = p.toArray();
  // double[] seArray = sqErr.toArray();

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
    // first calculate the mean
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
  Rcpp::Rcout << "Hello World";
  LinearRegression reg(data);
  Rcpp::Rcout << "TINKER";
  for (int i = 0; i < data.getNumColumns(); i++) {
    Variable* target = data.getVariable(i);
    Rcpp::Rcout << target->getName() << std::endl;
    std::vector<Variable*> regressors(data.getVariables());
    regressors.erase(regressors.begin() + i);
    RegressionResult* result = reg.regress(target, regressors);
    Rcpp::Rcout << *result;
  }
}
