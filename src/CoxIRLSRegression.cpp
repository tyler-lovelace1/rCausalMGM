// [[Rcpp::depends(BH)]]

#include "CoxIRLSRegression.hpp"
#include "CoxRegressionResult.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/normal.hpp>


CoxIRLSRegression::CoxIRLSRegression(DataSet& data) {
    this->data = data;
    dataMat = arma::mat(data.getData());
    variables = data.getVariables();
    rows = arma::uvec(data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

// CoxIRLSRegression::CoxIRLSRegression(CoxIRLSRegression &cr) {
//     this->data = cr.data;
//     dataMat = arma::mat(cr.dataMat);
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
// }

// CoxIRLSRegression::CoxIRLSRegression(CoxIRLSRegression &&cr) {
//     this->data = cr.data;
//     dataMat = cr.dataMat;
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
// }

// CoxIRLSRegression &CoxIRLSRegression::operator=(CoxIRLSRegression &cr) {
//     this->data = cr.data;
//     dataMat = arma::mat(cr.dataMat);
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
//     return *this;
// }

// CoxIRLSRegression &CoxIRLSRegression::operator=(CoxIRLSRegression &&cr) {
//     this->data = cr.data;
//     dataMat = cr.dataMat;
//     this->variables = this->data.getVariables();
//     this->rows = arma::uvec(this->data.getNumRows());
//     for (arma::uword i = 0; i < data.getNumRows(); i++)
//         rows[i] = i;
//     return *this;
// }

CoxRegressionResult CoxIRLSRegression::regress(const Node& target,
					       std::vector<Node>& regressors) {
    int n = this->rows.n_elem;
    int k = regressors.size();
    
    if (n < k)
	throw std::runtime_error("Cox regression ill-conditioned, samples less than regressors");

    int target_ = data.getColumn(target);

    arma::uvec regressors_ = arma::uvec(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	regressors_[i] = data.getColumn(regressors[i]);
    }

    // arma::vec target_vals = dataMat.col(target_);
    // target_vals = target_vals(this->rows);
    // arma::uvec censor = target.getCensorVec()(this->rows);
    // arma::uvec strata = target.getStrata()(this->rows);

    // Node Y(target);
    // Y.setCensor(target_vals, censor, strata);

    arma::mat X;

    if (regressors.empty()) {
	X = arma::mat(n, 1, arma::fill::randn);
    } else {
	X = dataMat.submat(this->rows, regressors_);
	X.insert_cols(0, arma::vec(n, arma::fill::ones));
    }

    double old_l, new_l;
    arma::vec beta(X.n_cols, arma::fill::zeros);
    arma::vec oldBeta(X.n_cols, arma::fill::zeros);
    arma::vec b = arma::vec((k > 0) ? k : 1);
    arma::vec se = arma::vec((k > 0) ? k : 1);
    arma::vec z = arma::vec((k > 0) ? k : 1);
    arma::vec pval = arma::vec((k > 0) ? k : 1);
    arma::vec res(X.n_rows, arma::fill::zeros);

    arma::vec diagHess(X.n_rows);
    arma::vec grad(X.n_rows);
    arma::vec eta(X.n_rows, arma::fill::zeros);

    if (regressors.size() == 0) {
	new_l = gradHess(eta, grad, diagHess, target);
	b[0] = 0.0;
	z[0] = 0.0;
	pval[0] = 1.0;
	se[0] = 1.0;
	return CoxRegressionResult({""}, n, b, z, pval, se, grad, this->alpha, new_l);
    }

    arma::vec Z(X.n_rows);
    arma::vec w(X.n_rows);

    arma::uvec nonzero;
    arma::mat Xprime;
    arma::vec Zprime;
    arma::vec sqrtW;
    
    // new_l = loss(beta, X, Y);

    double dbeta = 1;

    int iter = 0;

    auto start = std::chrono::high_resolution_clock::now();

    while (dbeta > 1e-8 || std::abs(new_l - old_l) > 1e-7) {

	double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	if (curr > 5 || iter > 100 || !beta.is_finite()) {
	    throw std::runtime_error("Cox Regression not converging");
	}
	    
	// old_l = new_l;
	oldBeta = beta;

	eta = X * beta;

	eta -= arma::mean(eta);

	old_l = gradHess(eta, grad, diagHess, target);

	// diagHess /= arma::mean(diagHess);
	
	w = -diagHess;

	diagHess.replace(0, -1e-10);

	Z = eta - grad / diagHess;

	nonzero = arma::find(w);
	sqrtW = arma::sqrt(w(nonzero));
	Xprime = arma::diagmat(sqrtW) * X.rows(nonzero);
	Zprime = sqrtW % Z(nonzero);

	// if (arma::rcond(Xprime.t() * Xprime) < 1e-15) {
	//     RcppThread::Rcout << "Cox regression of " << target.getName() << " given { ";
	//     for (const Node& node : regressors) {
	// 	RcppThread::Rcout << node.getName() << " ";
	//     }
	//     RcppThread::Rcout << "} is ill-conditioned\n";
	// }

	beta = arma::solve(Xprime, Zprime);

	// new_l = loss(beta, X, Y);

	new_l = 0.5 * arma::as_scalar((Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta)) - 0.5 * arma::sum(arma::square(grad) / diagHess) + old_l;

	dbeta = arma::norm(beta-oldBeta, 2) / arma::norm(beta, 2);
	iter++;

	// Rcpp::Rcout << "Iter:  " << iter << "    Loss:  " << new_l << "    |dx|/|x|:  "
	// 	    << dbeta << "    |dll|:  " << std::abs(new_l - old_l)
	// 	    << "\n    beta:  " << beta.t() << std::endl;
	
	// Rcpp::Rcout << "Iter:\t" << iter << "\n   Loss:\t" << new_l << "\n   beta:"
	// 	    << beta.t() << "\n   beta - betaOld:" << (beta-oldBeta).t()
	// 	    << "\n   dbeta:\t" << dbeta
	// 	    << "\n   IRLS Loss w/out constant:\t"
	// 	    << 0.5 * (Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta)
	// 	    << std::endl << std::endl;
		    // << "\n   Mean Weight:\t" << arma::mean(w)
		    // << "\n   IRLS Loss w/out constant:\t"
		    // << 0.5 * (Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta)
		    // << "\n   IRLS Loss w/ constant:\t"
		    // << 0.5 * (Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta) - 0.5 * arma::sum(arma::square(grad) / diagHess) + old_l
		    // << "\n   Taylor expansion Loss w/ constant:\t"
		    // << 0.5 * (X * beta - eta).t() * arma::diagmat(diagHess) * (X * beta - eta) + (X * beta - eta).t() * grad + old_l << std::endl;
    }

    // gradHess(beta, grad, hess, X, Y);

    // int df = n - k - 1;

    // w /= arma::mean(w);

    // censor = Y->getCensor();
    // arma::vec devResid = arma::sign(grad) % arma::sqrt(-2 * (grad + censor % arma::log(censor - grad)));

    // devResid.replace(arma::datum::nan, 0);

    // Rcpp::Rcout << "Martingale Residuals:" << grad.t() << std::endl;
    // Rcpp::Rcout << "Deviance Residuals:" << devResid.t() << std::endl;
    // Rcpp::Rcout << "Sum of Deviance / 2:    " << arma::sum(devResid)/2 << std::endl << std::endl;

    // arma::mat xTwxInv = arma::inv_sympd(X.t() * arma::diagmat(w) * X);

    // eta = X * beta;
    
    // res = Z - eta;

    // arma::vec scaledZ = eta + w % res;

    // Rcpp::Rcout << "eta" << eta.t()
    // 		<< "weight" << w.t()
    // 		<< "Z" << Z.t()
    // 		<< "scaledZ" << scaledZ.t() << std::endl;

    // double wrss = arma::as_scalar(res.t() * arma::diagmat(w) * res);

    // boost::math::students_t dist(df);
    
    // for (arma::uword i = 0; i < X.n_cols; i++) {
    // 	b[i] = beta[i];
    //     se[i] = std::sqrt(wrss * xTwxInv(i,i) / df);
    // 	if (se[i]== 0.0) se[i] = 1e-15;
    // 	z[i] = b[i] / se[i];
    // 	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    // }

    // std::vector<std::string> vNames(regressors.size());

    // for (int i = 0; i < regressors.size(); i++)	{
    // 	vNames[i] = regressors[i]->getName(); 
    // }

    beta.shed_row(0);

    X.shed_col(0);

    arma::mat info(k, k);
    
    infoMat(beta, info, X, target);

    arma::mat cov = -arma::inv(info);

    boost::math::normal dist;
    
    for (arma::uword i = 0; i < k; i++) {
	b[i] = beta[i];
	se[i] = std::sqrt(cov(i,i));
	if (se[i] == 0.0) se[i] = 1e-15;
	z[i] = b[i] / se[i];
	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    }

    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	vNames[i] = regressors[i].getName(); 
    }

    // res = w % res;

    arma::vec wz = w % Z;
    
    return CoxRegressionResult(vNames, n, b, z, pval, se, wz, alpha, new_l);
}

CoxRegressionResult CoxIRLSRegression::regress(const Node& target,
					       std::vector<Node>& regressors,
					       arma::uvec& _rows) {
    int n = _rows.n_elem;
    int k = regressors.size();
    
    if (n < k)
	throw std::runtime_error("Cox regression ill-conditioned, samples less than regressors");

    int target_ = data.getColumn(target);

    arma::uvec regressors_ = arma::uvec(regressors.size());

    // stringstream os;
    for (int i = 0; i < regressors.size(); i++)	{
	regressors_[i] = data.getColumn(regressors[i]);
	// os << regressors[i]->getName() << ","
    }
    // RcppThread::Rcout << os << std::endl; 

    arma::vec target_vals = dataMat.col(target_);
    target_vals = target_vals(_rows);
    arma::uvec censor = target.getCensorVec()(_rows);
    arma::uvec strata = target.getStrata()(_rows);

    // RcppThread::Rcout << "Values: " << target_vals.t();
    // RcppThread::Rcout << "Size of Values: " << target_vals.n_elem << std::endl;
    // RcppThread::Rcout << "Censor: " << censor.t();
    // RcppThread::Rcout << "Szie of Censor: " << censor.n_elem << std::endl;
    // RcppThread::Rcout << "Strata: " << strata.t();
    // RcppThread::Rcout << "Size of Strata: " << strata.n_elem << std::endl;
    
    Node Y(target);
    Y.setCensor(target_vals, censor, strata);

    arma::mat X;

    if (regressors.empty()) {
	X = arma::mat(n, 1, arma::fill::randn);
    } else {
	X = dataMat.submat(_rows, regressors_);
	X.insert_cols(0, arma::vec(n, arma::fill::ones));
    }

    double old_l, new_l;
    arma::vec beta(X.n_cols, arma::fill::zeros);
    arma::vec oldBeta(X.n_cols, arma::fill::zeros);
    arma::vec b = arma::vec((k > 0) ? k : 1);
    arma::vec se = arma::vec((k > 0) ? k : 1);
    arma::vec z = arma::vec((k > 0) ? k : 1);
    arma::vec pval = arma::vec((k > 0) ? k : 1);
    arma::vec res(X.n_rows, arma::fill::zeros);
    arma::vec diagHess(X.n_rows);
    arma::vec grad(X.n_rows);
    arma::vec eta(X.n_rows, arma::fill::zeros);

    if (regressors.size() == 0) {
	new_l = gradHess(eta, grad, diagHess, Y);
	b[0] = 0.0;
	z[0] = 0.0;
	pval[0] = 1.0;
	se[0] = 1.0;
	return CoxRegressionResult({""}, n, b, z, pval, se, grad, this->alpha, new_l);
    }

    arma::vec Z(X.n_rows);
    arma::vec w(X.n_rows);

    arma::uvec nonzero;
    arma::mat Xprime;
    arma::vec Zprime;
    arma::vec sqrtW;
    
    // new_l = loss(beta, X, Y);

    double dbeta = 1;

    int iter = 0;

    auto start = std::chrono::high_resolution_clock::now();

    while (dbeta > 1e-8 || std::abs(new_l - old_l) > 1e-7) {

	double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	if (curr > 5 || iter > 100 || !beta.is_finite()) {
	    throw std::runtime_error("Cox Regression not converging");
	}
	    
	// old_l = new_l;
	oldBeta = beta;

	eta = X * beta;

	eta -= arma::mean(eta);

	old_l = gradHess(eta, grad, diagHess, Y);

	// diagHess /= arma::mean(diagHess);
	
	w = -diagHess;

	diagHess.replace(0, -1e-10);

	Z = eta - grad / diagHess;

	nonzero = arma::find(w);
	sqrtW = arma::sqrt(w(nonzero));
	Xprime = arma::diagmat(sqrtW) * X.rows(nonzero);
	Zprime = sqrtW % Z(nonzero);

	// if (arma::rcond(Xprime.t() * Xprime) < 1e-15) {
	//     RcppThread::Rcout << "Cox regression of " << target.getName() << " given { ";
	//     for (const Node& node : regressors) {
	// 	RcppThread::Rcout << node.getName() << " ";
	//     }
	//     RcppThread::Rcout << "} is ill-conditioned\n";
	// }

	beta = arma::solve(Xprime, Zprime);

	// new_l = loss(beta, X, Y);

	new_l = 0.5 * arma::as_scalar((Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta)) - 0.5 * arma::sum(arma::square(grad) / diagHess) + old_l;

	dbeta = arma::norm(beta-oldBeta, 2) / arma::norm(beta, 2);
	iter++;

	// Rcpp::Rcout << "Iter:  " << iter << "    Loss:  " << new_l << "    |dx|/|x|:  "
	// 	    << dbeta << "    |dll|:  " << std::abs(new_l - old_l)
	// 	    << "\n    beta:  " << beta.t() << std::endl;
	// RcppThread::Rcout << "Iter:\t" << iter << "\n   Loss:\t" << new_l << "\n   beta:"
	// 		  << beta.t() << "\n   beta - betaOld:" << (beta-oldBeta).t()
	// 		  << "\n   dbeta:\t" << dbeta
	// 		  << "\n   IRLS Loss w/out constant:\t"
	// 		  << 0.5 * (Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta)
	// 		  << std::endl << std::endl;
		    // << "\n   Mean Weight:\t" << arma::mean(w)
		    // << "\n   IRLS Loss w/out constant:\t"
		    // << 0.5 * (Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta)
		    // << "\n   IRLS Loss w/ constant:\t"
		    // << 0.5 * (Z - X * beta).t() * arma::diagmat(diagHess) * (Z - X * beta) - 0.5 * arma::sum(arma::square(grad) / diagHess) + old_l
		    // << "\n   Taylor expansion Loss w/ constant:\t"
		    // << 0.5 * (X * beta - eta).t() * arma::diagmat(diagHess) * (X * beta - eta) + (X * beta - eta).t() * grad + old_l << std::endl;
    }

    // gradHess(beta, grad, hess, X, Y);

    // int df = n - k;

    // w /= arma::mean(w);

    // arma::mat xTwxInv = arma::inv_sympd(X.t() * arma::diagmat(w) * X);

    // devResid = arma::zeros(n);

    // eta = X * beta;
    // eta -= arma::mean(eta);
    
    // arma::vec cumHaz(n, arma::fill::zeros);

    // censor = Y->getCensor();
    // arma::vec devResid = arma::sign(grad) % arma::sqrt(-2 * (grad + censor % arma::log(censor - grad)));

    // devResid.replace(arma::datum::nan, 0);

    // Rcpp::Rcout << "Martingale Residuals:" << grad.t() << std::endl;
    // Rcpp::Rcout << "Deviance Residuals:" << devResid.t() << std::endl;
    // Rcpp::Rcout << "Sum of Deviance / 2:   " << arma::sum(devResid)/2 << std::endl << std::endl;

    // eta = X * beta;
    
    // res = Z - eta;

    // arma::vec scaledZ = eta + w % res;

    // Rcpp::Rcout << "eta:  " << eta.t() << std::endl;
    // Rcpp::Rcout << "Weight:  " << w.t() << std::endl;
    // Rcpp::Rcout << "Z:  " << Z.t() << std::endl;
    // Rcpp::Rcout << "Scaled Z:  " << scaledZ.t() << std::endl;

    // Rcpp::Rcout << "eta" << eta.t()
    // 		<< "weight" << w.t()
    // 		<< "Z" << Z.t()
    // 		<< "scaledZ" << scaledZ.t() << std::endl;


    // double wrss = arma::as_scalar(res.t() * arma::diagmat(w) * res);

    // beta = beta.subvec(1, k);

    // arma::vec(beta.subvec(1, k));

    // X.shed_col(0);

    // arma::mat info(X.n_cols, X.n_cols);
    
    // infoMat(beta, info, X, Y);

    // arma::mat cov = -arma::inv(info);

    // boost::math::normal dist;
    
    // for (arma::uword i = 0; i < X.n_cols; i++) {
    // 	b[i] = beta[i];
    // 	se[i] = std::sqrt(cov(i,i));
    // 	if (se[i] == 0.0) se[i] = 1e-15;
    // 	z[i] = b[i] / se[i];
    // 	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    // }

    beta.shed_row(0);

    X.shed_col(0);

    arma::mat info(k, k);
    
    infoMat(beta, info, X, Y);

    arma::mat cov = -arma::inv(info);

    // RcppThread::Rcout << "Covariance:\n" << cov << std::endl; 

    boost::math::normal dist;
    
    for (arma::uword i = 0; i < k; i++) {
	b[i] = beta[i];
	se[i] = std::sqrt(cov(i,i));
	if (se[i] == 0.0) se[i] = 1e-15;
	z[i] = b[i] / se[i];
	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    }

    
    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	vNames[i] = regressors[i].getName(); 
    }

    // res = w % res;

    arma::vec wz = w % Z;
    
    return CoxRegressionResult(vNames, n, b, z, pval, se, wz, alpha, new_l);
}



double CoxIRLSRegression::loss(arma::vec& beta, arma::mat& X, Node target) {

    int numStrata = target.getNumStrata();
    double loss = 0.0;

    for (int strat = 0; strat < numStrata; strat++) {

	arma::uvec index = target.getIndex(strat);
	arma::uvec order = target.getOrder(strat);
	arma::uvec H = target.getH(strat);
	arma::uvec censor = target.getCensor(strat);
	double HsumTheta, m, sub, d;

	int n = index.n_elem;

	arma::vec logtheta = X.rows(index) * beta;
  
	arma::vec theta = arma::exp(logtheta);
	double rs_sum = arma::accu(theta);

	int i = 0;
	for (int j = 0; j < H.n_elem; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;
	    for (int k = 0; k < H[j]; k++) {
		if (censor[order[i+k]]) {
		    m++;
		    loss += logtheta[order[i+k]];
		    HsumTheta += theta[order[i+k]];
		}
		sub += theta[order[i+k]];
	    }

	    if (HsumTheta - sub > 1e-5)
		throw std::runtime_error("Error in Cox regression loss, HsumTheta > sub: " +
					 std::to_string(HsumTheta) + " > " + std::to_string(sub));

	    if (sub - rs_sum > 1e-5) {
		if ((H[j] + i) != n) {
		    rs_sum = arma::accu(theta(order.subvec(i, n-1)));
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Cox regression loss, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }

	    for (int l = 0; l < m; l++) {
		d = l / ((double) m);
		loss -= std::log(rs_sum - d * HsumTheta);
	    }

	    i += H[j];
	    rs_sum -= sub;
	}
    }

    return loss;
}

double CoxIRLSRegression::gradHess(arma::vec& eta, arma::vec& grad, arma::vec& diagHess,
				   Node target) {

    int numStrata = target.getNumStrata();
    double loss = 0.0;

    grad.fill(0);
    diagHess.fill(0);

    for (int strat = 0; strat < numStrata; strat++) {

	arma::uvec index = target.getIndex(strat);
	arma::uvec order = target.getOrder(strat);
	arma::uvec H = target.getH(strat);
	arma::uvec censor = target.getCensor(strat);
	double HsumTheta, m, sub, d, phi, dSum, dSum2;
	double theta_weight_sum = 0, theta_weight2_sum = 0;

	int n = index.n_elem;

	// eta -= arma::mean(eta);

	arma::vec theta = arma::exp(eta(index));
	arma::vec theta_weight = arma::zeros(n);
	arma::vec theta_weight2 = arma::zeros(n);
	double rs_sum = arma::accu(theta);

	grad(index) += arma::conv_to<arma::vec>::from(censor);
    
	int i = 0;
	for (int j = 0; j < H.n_elem; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;
	    dSum = 0;
	    dSum2 = 0;

	    for (int k = 0; k < H[j]; k++) {
		if (censor[order[i+k]]) {
		    m++;
		    HsumTheta += theta[order[i+k]];
		    loss += eta[index[order[i+k]]];
		}
		sub += theta[order[i+k]];
	    }

	    if (m > 0) {
		if (sub - rs_sum > 1e-5) {
		    if ((H[j] + i) != n) {
			rs_sum = arma::accu(theta(order.subvec(i, n-1)));
			if (sub - rs_sum > 1e-5) {
			    throw std::runtime_error("Error in Cox IRLS Regression gradHess, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
			}
		    } else {
			sub = rs_sum;
		    }
		}

		for (int l = 0; l < m; l++) {
		    d = l / ((double) m);
		    phi = rs_sum - d * HsumTheta;
		    loss -= std::log(phi);
		    theta_weight_sum += 1 / phi;
		    dSum += d / phi;
		    theta_weight2_sum += 1 / std::pow(phi, 2);
		    dSum2 += (2*d - std::pow(d, 2)) / std::pow(phi, 2);
		}
	    }

	    for (int k = 0; k < H[j]; k++) {
		theta_weight[order[i+k]] = theta_weight_sum - censor[order[i+k]] * dSum;
		theta_weight2[order[i+k]] = theta_weight2_sum - censor[order[i+k]] * dSum2;
	    }
	
	    i += H[j];
	    rs_sum -= sub;
	}

	grad(index) -= theta % theta_weight;

	diagHess(index) = arma::square(theta) % theta_weight2 - theta % theta_weight;

    }

    // diagHess.replace(0, -1e-10);

    return loss;
}

void CoxIRLSRegression::infoMat(arma::vec& beta, arma::mat& hess,
				arma::mat& X, Node target) {

    int numStrata = target.getNumStrata();
    hess.fill(0);

    for (int strat = 0; strat < numStrata; strat++) {

	arma::uvec index = target.getIndex(strat);
	arma::uvec order = target.getOrder(strat);
	arma::uvec H = target.getH(strat);
	arma::uvec censor = target.getCensor(strat);
	double HsumTheta, m, sub, d, phi;

	arma::mat Xsub = X.rows(index);

	int n = index.n_elem;

	arma::vec theta = arma::exp(Xsub * beta);
	double rs_sum = arma::accu(theta);

	arma::mat temp(beta.n_elem, beta.n_elem, arma::fill::zeros);
	arma::vec sub_num(beta.n_elem), HsumThetaVec(beta.n_elem), Z(beta.n_elem);
	arma::mat sub_outer(beta.n_elem, beta.n_elem), HsumOuter(beta.n_elem, beta.n_elem);

	arma::vec num = arma::sum(arma::diagmat(theta) * Xsub, 0).t();
	arma::mat outer_num = Xsub.t() * arma::diagmat(theta) * Xsub;

	int i = 0;
	for (int j = 0; j < H.n_elem; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;

	    HsumThetaVec.fill(0);
	    sub_num.fill(0);

	    HsumOuter.fill(0);
	    sub_outer.fill(0);

	    for (int k = 0; k < H[j]; k++) {
		temp = theta[order[i+k]] * Xsub.row(order[i+k]).t() * Xsub.row(order[i+k]);
	    
		if (censor[order[i+k]]) {
		    m ++;

		    HsumTheta += theta[order[i+k]];
		    // grad += X.row(order[i+k]);
		    HsumThetaVec += theta[order[i+k]] * Xsub.row(order[i+k]);
		    HsumOuter += temp;
		}

		sub += theta[order[i+k]];
		sub_num += theta[order[i+k]] * Xsub.row(order[i+k]);
		sub_outer += temp;
	    }

	    if (HsumTheta - sub > 1e-5)
		throw std::runtime_error("Error in Cox regression gradHess, HsumTheta > sub: " +
					 std::to_string(HsumTheta) + " > " + std::to_string(sub));

	    if (sub - rs_sum > 1e-5) {
		if ((H[j] + i) != n) {
		    rs_sum = arma::accu(theta(order.subvec(i, n-1)));
		    num = arma::sum(arma::diagmat(theta(order.subvec(i, n-1)))
				    * Xsub.rows(order.subvec(i, n-1)), 0).t();
		    outer_num = (Xsub.rows(order.subvec(i, n-1)).t()
				 * arma::diagmat(theta(order.subvec(i, n-1)))
				 * Xsub.rows(order.subvec(i, n-1)));
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Cox regression gradHess, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }
	
	    // if (sub - rs_sum > 1e-5)
	    //     throw std::runtime_error("Error in Cox regression gradHess, sub > rs_sum: " +
	    // 			     std::to_string(sub) + " > " + std::to_string(rs_sum));

	    for (int l = 0; l < m; l++) {
		d = l / ((double) m);
		Z = num - d * HsumThetaVec;
		phi = rs_sum - d * HsumTheta;
		// grad -= Z / phi;
		hess -= ((outer_num - d * HsumOuter) / phi - Z * Z.t() / (phi*phi));
	    }

	    i += H[j];
	    rs_sum -= sub;
	    num -= sub_num;
	    outer_num -= sub_outer;
	}
    }
}

// no export // [[Rcpp::export]]
void CoxIRLSRegressionTest(const Rcpp::DataFrame& df, std::string targetName,
			   std::vector<std::string>& regressorNames, int repetitions) {
    DataSet data = DataSet(df);
    data.dropMissing();

    repetitions = std::max(1, repetitions);

    Rcpp::Rcout << "-----START----- \n";
    Node target = data.getVariable(targetName);
    std::vector<Node> inputRegressors;
    for (std::string varName : regressorNames) {
	inputRegressors.push_back(data.getVariable(varName));
    }

    std::vector<Node> regressors;
    for (Node var : inputRegressors) {
	// Rcpp::Rcout << var << std::endl;
	if (var.isContinuous()) {
	    regressors.push_back(var);
	    continue;
	}

	if (var.isDiscrete() && var.getNumCategories() < 3) {
	    regressors.push_back(var);
	    continue;
	}

	if (!var.isDiscrete()) {
	    throw std::invalid_argument("*Invalid variable type*");
	}

	std::vector<std::string> varCats = var.getCategories();

	std::vector<Node> variables;
	/*********************************************************************/
	std::string temp = var.getName();
	for (auto it = varCats.begin() + 1; it != varCats.end(); it++) {
	    Node newVar = Node(new DiscreteVariable(temp + "MULTINOM." + *it, 2));

	    /*********************************************************************/

	    variables.push_back(newVar);

	    // Rcpp::Rcout << newVar << std::endl;

	    data.addVariable(newVar);

	    int newVarIndex = data.getColumn(newVar);
	    int numCases = data.getNumRows();

	    for (int l = 0; l < numCases; l++) {
		int dataCellIndex = data.getInt(l, data.getColumn(var));
		if (dataCellIndex == var.getIndex(*it)) {
		    data.set(l, newVarIndex, 1);
		}
		else {
		    data.set(l, newVarIndex, 0);
		}
	    }
	}
	regressors.insert(regressors.end(), variables.begin(), variables.end());
    }

    CoxIRLSRegression reg(data);
    
    CoxRegressionResult result;
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < repetitions; i++) {
	result = reg.regress(target, regressors);
    }

    double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
    Rcpp::Rcout << result;
    Rcpp::Rcout << "-----END----- \n";
    
    Rcpp::Rcout << "-----START----- \n";
    arma::uvec rows = arma::regspace<arma::uvec>(0,data.getNumRows()/2-1);

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < repetitions; i++) {
	result = reg.regress(target, regressors, rows);
    }
    curr += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
    
    Rcpp::Rcout << result;
    Rcpp::Rcout << "-----END----- \n";

    Rcpp::Rcout << "Elapsed Time: " << curr << " seconds.\n";
    
}
