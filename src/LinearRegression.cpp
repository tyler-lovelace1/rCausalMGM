// [[Rcpp::depends(BH)]]

#include "LinearRegression.hpp"
#include "RegressionResult.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>

#include <fstream>

LinearRegression::LinearRegression(DataSet &data)
{
    this->data = data;
    dataMat = arma::mat(data.getData());
    variables = data.getVariables();
    rows = arma::uvec(data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

// LinearRegression::LinearRegression(LinearRegression &lr)
// {
//     this->data = lr.data;
//     dataMat = arma::mat(lr.dataMat);
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
// }

// LinearRegression::LinearRegression(LinearRegression &&lr)
// {
//     this->data = lr.data;
//     dataMat = lr.dataMat;
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
// }

// LinearRegression &LinearRegression::operator=(LinearRegression &lr)
// {
//     this->data = lr.data;
//     dataMat = arma::mat(lr.dataMat);
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
//     return *this;
// }

// LinearRegression &LinearRegression::operator=(LinearRegression &&lr)
// {
//     this->data = lr.data;
//     dataMat = lr.dataMat;
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
//     return *this;
// }

RegressionResult LinearRegression::regress(const Node& target, std::vector<Node>& regressors) {

    int n = rows.n_elem;
    int k = regressors.size() + 1;

    if (n < k)
	throw std::runtime_error("Linear regression ill-conditioned, samples less than regressors");

    int target_ = data.getColumn(target);

    arma::uvec regressors_ = arma::uvec(regressors.size());


    for (int i = 0; i < regressors.size(); i++)
        regressors_[i] = data.getColumn(regressors[i]);

    arma::uvec target_Vec(1);
    target_Vec.fill(target_);
    arma::mat y = dataMat.submat(rows, target_Vec);

    arma::mat xSub = dataMat.submat(rows, regressors_);

    arma::mat x = arma::mat(xSub.n_rows, xSub.n_cols + 1);
    
    for (arma::uword i = 0; i < x.n_rows; i++) {
	for (arma::uword j = 0; j < x.n_cols; j++) {
	    if (j == 0)	{
		x(i, j) = 1;
	    }
	    else {
		x(i, j) = xSub(i, j - 1);
	    }
	}
    }

    arma::vec b = arma::solve(x, y);
    arma::mat xTxInv = arma::inv_sympd(x.t() * x);
    
    arma::mat yHat = x * b;
    arma::mat res = y - yHat;

    arma::vec yHat_ = yHat.col(0);
    arma::vec res_ = res.col(0);

    /* NOT CURRENTLY IN USE, ONLY NEEDED FOR LRT
    
    // arma::mat b2 = arma::mat(b);
    // arma::mat yHat2 = x * b2;

    // arma::mat res2 = y - yHat2;
    // this->res2 = res2.col(0);

    */

    double rss_ = rss(x, y, b);
    double se = std::sqrt(rss_ / (n - k));
    double tss_ = tss(y);
    double r2 = 1.0 - (rss_ / tss_);

    arma::vec sqErr = arma::vec(x.n_cols);
    arma::vec t = arma::vec(x.n_cols);
    arma::vec p = arma::vec(x.n_cols);

    boost::math::students_t dist(n - k);
    
    for (arma::uword i = 0; i < x.n_cols; i++) {
	double s_ = se * se * xTxInv(i, i);
	double se_ = std::sqrt(s_);
	if (se_== 0.0) se_ = 1e-10;
	double t_ = b(i) / se_;
	double p_ = 2 * (1.0 - boost::math::cdf(dist, std::abs(t_)));

	// if (i == 1)
	//     {
	// 	// Rcpp::Rcout << "beta = " << b(i,0) << std::endl;
	// 	// Rcpp::Rcout << "SE = " << se_ << std::endl;
	// 	// Rcpp::Rcout << "t-statistic = " << t_ << std::endl;
	// 	// Rcpp::Rcout << "p-value = " << p_ << std::endl;
	//     }

	sqErr[i] = se_;
	t[i] = t_;
	p[i] = p_;
    }

    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)
        vNames[i] = regressors[i].getName(); // getName Function may not be implemented


    return RegressionResult(regressors.size() == 0, vNames, n, b, t, p, sqErr, r2, rss_, alpha, yHat_, res_); // MUST CONVERT B INTO A VECTOR
}

RegressionResult LinearRegression::regress(const Node& target,
					   std::vector<Node>& regressors,
					   arma::uvec& _rows)
{

    int n = _rows.n_elem;
    int k = regressors.size() + 1;

    if (n < k)
	throw std::runtime_error("Linear regression ill-conditioned, samples less than regressors");

    int target_ = data.getColumn(target);

    arma::uvec regressors_ = arma::uvec(regressors.size());


    for (int i = 0; i < regressors.size(); i++)
        regressors_[i] = data.getColumn(regressors[i]);

    arma::uvec target_Vec(1);
    target_Vec.fill(target_);
    arma::mat y = dataMat.submat(_rows, target_Vec);

    arma::mat xSub = dataMat.submat(_rows, regressors_);

    arma::mat x = arma::mat(xSub.n_rows, xSub.n_cols + 1);
    
    for (arma::uword i = 0; i < x.n_rows; i++) {
	for (arma::uword j = 0; j < x.n_cols; j++) {
	    if (j == 0)	{
		x(i, j) = 1;
	    }
	    else {
		x(i, j) = xSub(i, j - 1);
	    }
	}
    }

    arma::vec b = arma::solve(x, y);
    arma::mat xTxInv = arma::inv_sympd(x.t() * x);
    
    arma::mat yHat = x * b;
    arma::mat res = y - yHat;

    arma::vec yHat_ = yHat.col(0);
    arma::vec res_ = res.col(0);

    /* NOT CURRENTLY IN USE, ONLY NEEDED FOR LRT
    
    // arma::mat b2 = arma::mat(b);
    // arma::mat yHat2 = x * b2;

    // arma::mat res2 = y - yHat2;
    // this->res2 = res2.col(0);

    */

    double rss_ = rss(x, y, b);
    double se = std::sqrt(rss_ / (n - k));
    double tss_ = tss(y);
    double r2 = 1.0 - (rss_ / tss_);

    arma::vec sqErr = arma::vec(x.n_cols);
    arma::vec t = arma::vec(x.n_cols);
    arma::vec p = arma::vec(x.n_cols);

    boost::math::students_t dist(n - k);
    
    for (arma::uword i = 0; i < x.n_cols; i++) {
	double s_ = se * se * xTxInv(i, i);
	double se_ = std::sqrt(s_);
	if (se_== 0.0) se_ = 1e-10;
	double t_ = b(i) / se_;
	double p_ = 2 * (1.0 - boost::math::cdf(dist, std::abs(t_)));

	// if (i == 1)
	//     {
	// 	// Rcpp::Rcout << "beta = " << b(i,0) << std::endl;
	// 	// Rcpp::Rcout << "SE = " << se_ << std::endl;
	// 	// Rcpp::Rcout << "t-statistic = " << t_ << std::endl;
	// 	// Rcpp::Rcout << "p-value = " << p_ << std::endl;
	//     }

	sqErr[i] = se_;
	t[i] = t_;
	p[i] = p_;
    }

    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)
        vNames[i] = regressors[i].getName(); // getName Function may not be implemented

    return RegressionResult(regressors.size() == 0, vNames, n, b, t, p, sqErr, r2, rss_, alpha, yHat_, res_); // MUST CONVERT B INTO A VECTOR
}

double LinearRegression::rss(const arma::mat& x, const arma::vec& y, const arma::vec& b)
{
    double rss = 0.0;

    for (arma::uword i = 0; i < x.n_rows; i++)
    {
        double yH = 0.0;

        for (arma::uword j = 0; j < x.n_cols; j++)
        {
            yH += b(j) * x(i, j);
        }

        double d = y(i) - yH;

        rss += d * d;
    }

    return rss;
}

double LinearRegression::tss(const arma::vec& y)
{
    double mean = 0.0;

    for (arma::uword i = 0; i < y.n_rows; i++)
    {
        mean += y(i);
    }

    mean /= (double)(y.n_rows);

    double ssm = 0.0;

    for (arma::uword i = 0; i < y.n_rows; i++)
    {
        double d = mean - y(i);
        ssm += d * d;
    }
    return ssm;
}

// void LinearRegressionTest(const Rcpp::DataFrame& df) {
//   DataSet data = DataSet(df, 5);
//   LinearRegression reg(data);

//   for (int i = 0; i < data.getNumColumns(); i++) {
//     Rcpp::Rcout << "-----START----- \n";
//     const Node& target = data.getVariable(i);
//     Rcpp::Rcout << target->getName() << std::endl;
//     std::vector<Node> regressors(data.getVariables());
//     regressors.erase(regressors.begin() + i);
//     RegressionResult* result = reg.regress(target, regressors);
//     Rcpp::Rcout << *result;
//     Rcpp::Rcout << "-----END----- \n";
//   }
// }
