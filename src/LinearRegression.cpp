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

LinearRegression::LinearRegression(LinearRegression &lr)
{
    this->data = lr.data;
    dataMat = arma::mat(lr.dataMat);
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

LinearRegression::LinearRegression(LinearRegression &&lr)
{
    this->data = lr.data;
    dataMat = lr.dataMat;
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

LinearRegression &LinearRegression::operator=(LinearRegression &lr)
{
    this->data = lr.data;
    dataMat = arma::mat(lr.dataMat);
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
    return *this;
}

LinearRegression &LinearRegression::operator=(LinearRegression &&lr)
{
    this->data = lr.data;
    dataMat = lr.dataMat;
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
    return *this;
}

RegressionResult LinearRegression::regress(Variable *target, std::vector<Variable *>& regressors)
{
    // std::ofstream logfile;
    // logfile.open("lin_reg_debug.log", std::ios_base::app);

    int n = rows.size();
    int k = regressors.size() + 1;

    int target_ = data.getColumn(target);

    arma::uvec regressors_ = arma::uvec(regressors.size());

    for (int i = 0; i < regressors.size(); i++)
    {
        regressors_[i] = data.getColumn(regressors[i]);
    }

    // if (target_ == -1)
    // {
    //     // Rcpp::Rcout <<"\n";
    // }

    arma::uvec target_Vec(1);
    target_Vec.fill(target_);
    arma::mat y = dataMat.submat(rows, target_Vec);

    arma::mat xSub = dataMat.submat(rows, regressors_);

    arma::mat x;

    // if (regressors.size() > 0)
    // {
    x = arma::mat(xSub.n_rows, xSub.n_cols + 1);
    
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
    // }
    // else
    // {
    //     x = arma::mat(xSub.n_rows, xSub.n_cols);

    //     for (arma::uword i = 0; i < x.n_rows; i++)
    //     {
    //         for (arma::uword j = 0; j < x.n_cols; j++)
    //         {
    //             x(i, j) = xSub(i, j);
    //         }
    //     }
    // }

    arma::mat xT = x.t();
    arma::mat xTx = xT * x;
    arma::mat xTxInv = xTx.i();
    arma::mat xTy = xT * y;
    arma::mat b = xTxInv * xTy;

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

    double rss_ = LinearRegression::rss(x, y, b);
    double se = std::sqrt(rss_ / (n - k));
    double tss_ = LinearRegression::tss(y);
    double r2 = 1.0 - (rss_ / tss_);

    arma::vec sqErr = arma::vec(x.n_cols);
    arma::vec t = arma::vec(x.n_cols);
    arma::vec p = arma::vec(x.n_cols);

    boost::math::students_t dist(n - k);
    
    // std::unique_lock<std::mutex> linregLock(linregMutex);
    // logfile.open("lin_reg_debug.log", std::ios_base::app);
    // logfile.open("../debug.log", std::ios_base::app);
    // logfile << "LINEAR REGRESSION:\t" << target->getName() << " ? "
    // 	    << regressors[0]->getName() <<  " | [";
    // for (arma::uword i = 1; i < regressors.size(); i++) logfile << regressors[i]->getName() << ",";
    // logfile << "]\n  dist.df = " << dist.degrees_of_freedom() << std::endl;
    
    for (arma::uword i = 0; i < x.n_cols; i++) {
	double s_ = se * se * xTxInv(i, i);
	double se_ = std::sqrt(s_);
	double t_ = b(i, 0) / se_;

	// if (i == 0) logfile << "    var_ = intercept" << std::endl;
	// else logfile << "    var_ = " << regressors[i-1]->getName() << std::endl;
	
	// logfile << "    b_ = " << b(i, 0) << std::endl;
	// logfile << "    se_ = " << se_ << std::endl;
	// // logfile << "    s_ = " << s_ << std::endl;
	// logfile << "    t_ = " << t_ << std::endl;
	double p_ = 2 * (1.0 - boost::math::cdf(dist, std::abs(t_)));
	// logfile << "    p_ = " << p_ << std::endl;

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

    // logfile.close();

    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)
    {
        // Rcpp::Rcout << regressors[i]->getName() << "\n";
        vNames[i] = regressors[i]->getName(); // getName Function may not be implemented
    }

    // arma::vec bArray = b.columns() == 0 ? new double[0] : b.getColumn(0).toArray(); // double check,
    //dealing with case where we dont give it anything to regess on
    // logfile.close();

    return RegressionResult(regressors.size() == 0, vNames, n, b, t, p, sqErr, r2, rss_, alpha, yHat_, res_); // MUST CONVERT B INTO A VECTOR
}

RegressionResult LinearRegression::regress(Variable *target,
					   std::vector<Variable*>& regressors,
					   arma::uvec& _rows)
{
    // std::ofstream logfile;
    // logfile.open("lin_reg_debug.log", std::ios_base::app);

    int n = _rows.n_elem;
    int k = regressors.size() + 1;

    if (n < k)
	throw std::runtime_error("Linear regression ill-conditioned, samples less than regressors");

    // logfile << target->getName() << std::endl;

    int target_ = data.getColumn(target);

    // logfile << target_ << std::endl;

    // logfile << "regressors.size() = " << regressors.size() << std::endl;

    arma::uvec regressors_ = arma::uvec(regressors.size());

    // std::string names = "\t";

    for (int i = 0; i < regressors.size(); i++)
    {
	// names += regressors[i]->getName() + '\t';
	// logfile << i << " in loop" << std::endl;
	// logfile << regressors[i] << '\t';
	// logfile << regressors[i]->getName() + '\t';
        regressors_[i] = data.getColumn(regressors[i]);
	// logfile << regressors_[i] + '\t';
    }

    // logfile << "\n";

    // if (target_ == -1)
    // {
    //     // Rcpp::Rcout <<"\n";
    // }

    arma::uvec target_Vec(1);
    target_Vec.fill(target_);
    arma::mat y = dataMat.submat(_rows, target_Vec);

    arma::mat xSub = dataMat.submat(_rows, regressors_);

    // logfile << names << std::endl << xSub.rows(0,5) << std::endl << std::endl << xSub.rows(n-6,n-1) << std::endl;

    // arma::mat x;

    // if (regressors.size() > 0)
    // {
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
    // }
    // else
    // {
    //     x = arma::mat(xSub.n_rows, xSub.n_cols);

    //     for (arma::uword i = 0; i < x.n_rows; i++)
    //     {
    //         for (arma::uword j = 0; j < x.n_cols; j++)
    //         {
    //             x(i, j) = xSub(i, j);
    //         }
    //     }
    // }

    // const arma::mat xT = std::move(x.t());
    // const arma::mat xTx = std::move(xT * x);
    arma::vec b = arma::solve(x, y);
    arma::mat xTxInv = arma::inv_sympd(x.t() * x); // std::move(xTx.i());
    // const arma::mat xTy = std::move(xT * y);
    // const arma::mat b = std::move(xTxInv * xTy);

    arma::mat yHat = x * b;
    arma::mat res = y - yHat;

    arma::vec yHat_ = yHat.col(0);
    arma::vec res_ = res.col(0);
    // const arma::vec b_ = b.col(0);

    // logfile << '\t' << names << std::endl << b.t() << std::endl << std::endl;

    /* NOT CURRENTLY IN USE, ONLY NEEDED FOR LRT
    
    // arma::mat b2 = arma::mat(b);
    // arma::mat yHat2 = x * b2;

    // arma::mat res2 = y - yHat2;
    // this->res2 = res2.col(0);

    */

    double rss_ = LinearRegression::rss(x, y, b);
    double se = std::sqrt(rss_ / (n - k));
    double tss_ = LinearRegression::tss(y);
    double r2 = 1.0 - (rss_ / tss_);

    arma::vec sqErr = arma::vec(x.n_cols);
    arma::vec t = arma::vec(x.n_cols);
    arma::vec p = arma::vec(x.n_cols);

    boost::math::students_t dist(n - k);
    
    // std::unique_lock<std::mutex> linregLock(linregMutex);
    // logfile.open("debug.log", std::ios_base::app);
    // logfile << "LINEAR REGRESSION:\t" << target->getName() << " ? "
    // 	<< regressors[0]->getName() <<  " | [";
    // for (arma::uword i = 1; i < regressors.size(); i++) logfile << regressors[i]->getName() << ",";
    // logfile << "]\n  dist.df = " << dist.degrees_of_freedom() << std::endl;
    
    for (arma::uword i = 0; i < x.n_cols; i++) {
	double s_ = se * se * xTxInv(i, i);
	double se_ = std::sqrt(s_);
	if (se_== 0.0) se_ = 1e-10;
	double t_ = b(i) / se_;
	double p_ = 2 * (1.0 - boost::math::cdf(dist, std::abs(t_)));

	// if (i == 0) logfile << "    var_ = intercept" << std::endl;
	// else logfile << "    var_ = " << regressors[i-1]->getName() << std::endl;
	
	// logfile << "    b_ = " << b(i, 0) << std::endl;
	// logfile << "    se_ = " << se_ << std::endl;
	// // logfile << "    s_ = " << s_ << std::endl;
	// logfile << "    t_ = " << t_ << std::endl;	
	// logfile << "    p_ = " << p_ << std::endl;

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
    // logfile.close();
    

    // logfile << '\t' << names << std::endl << x.rows(0,5) << std::endl << std::endl << xSub.rows(n-6,n-1) << std::endl << std::endl;

    // logfile.close();

    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)
    {
        // Rcpp::Rcout << regressors[i]->getName() << "\n";
        vNames[i] = regressors[i]->getName(); // getName Function may not be implemented
    }

    // arma::vec bArray = b.columns() == 0 ? new double[0] : b.getColumn(0).toArray(); // double check,
    //dealing with case where we dont give it anything to regess on
    // logfile.close();

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
//     Variable* target = data.getVariable(i);
//     Rcpp::Rcout << target->getName() << std::endl;
//     std::vector<Variable*> regressors(data.getVariables());
//     regressors.erase(regressors.begin() + i);
//     RegressionResult* result = reg.regress(target, regressors);
//     Rcpp::Rcout << *result;
//     Rcpp::Rcout << "-----END----- \n";
//   }
// }
