// [[Rcpp::depends(BH)]]

#include "CoxRegression.hpp"
#include "CoxRegressionResult.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>

#include <fstream>

CoxRegression::CoxRegression(DataSet& data) {
    this->data = data;
    dataMat = arma::mat(data.getData());
    variables = data.getVariables();
    rows = arma::uvec(data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

CoxRegression::CoxRegression(CoxRegression &cr) {
    this->data = cr.data;
    dataMat = arma::mat(cr.dataMat);
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

CoxRegression::CoxRegression(CoxRegression &&cr) {
    this->data = cr.data;
    dataMat = cr.dataMat;
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

CoxRegression &CoxRegression::operator=(CoxRegression &cr) {
    this->data = cr.data;
    dataMat = arma::mat(cr.dataMat);
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
    return *this;
}

CoxRegression &CoxRegression::operator=(CoxRegression &&cr) {
    this->data = cr.data;
    dataMat = cr.dataMat;
    this->variables = this->data.getVariables();
    this->rows = arma::uvec(this->data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
    return *this;
}

CoxRegressionResult CoxRegression::regress(CensoredVariable* target,
					   std::vector<Variable*>& regressors) {
    int n = this->rows.n_elem;
    int k = regressors.size();

    if (n < k)
	throw std::runtime_error("Cox regression ill-conditioned, samples less than regressors");

    int target_ = data.getColumn(target);

    arma::uvec regressors_ = arma::uvec(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	regressors_[i] = data.getColumn(regressors[i]);
    }

    arma::vec target_vals = dataMat.col(target_);
    target_vals = target_vals(this->rows);
    arma::uvec censor = target->getCensor()(this->rows);

    CensoredVariable* Y = new CensoredVariable(*target);
    Y->setCensor(target_vals, censor);

    arma::mat X;

    if (regressors.empty())
	X = arma::mat(n, 1, arma::fill::randn);
    else
	X = dataMat.submat(this->rows, regressors_);

    double old_l, new_l, a, m, t, c = 0.1;
    arma::vec beta(X.n_cols, arma::fill::zeros);
    arma::vec b = beta;
    arma::vec se = arma::vec(X.n_cols);
    arma::vec z = arma::vec(X.n_cols);
    arma::vec pval = arma::vec(X.n_cols);

    if (regressors.size() == 0) {
	new_l = loss(std::move(beta), X, Y);
	b[0] = 0.0;
	z[0] = 0.0;
	pval[0] = 1.0;
	se[0] = 1.0;
	return CoxRegressionResult({""}, n, b, z, pval, se, this->alpha, new_l);
    }

    arma::mat hess(X.n_cols, X.n_cols);
    arma::vec grad(X.n_cols), p(X.n_cols);
    
    new_l = loss(std::move(beta), X, Y);

    double dbeta = 1;

    int iter = 0;

    auto start = std::chrono::high_resolution_clock::now();

    while (dbeta > 1e-5) {

	double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	if (curr > 5) {
	    throw std::runtime_error("Cox Regression not converging");
	}
	    
	old_l = new_l;
	a = 1;

	gradHess(beta, grad, hess, X, Y);

	p = arma::inv(hess) * grad;
	p /= -arma::norm(p, 2);
	
	m = arma::dot(p, grad);

	if (m <= 0)
	    throw std::runtime_error("Error in Cox regression, m <= 0");

	t = -c * m;

	while (true) {
	    new_l = loss(beta + a*p, X, Y);
	    if (old_l - new_l <= a*t) {
		if (!std::isinf(new_l)) {
		    break;
		}
	    }
	    a /= 2;
	}

	beta += a*p;
	dbeta = arma::norm(a*p, 2) / arma::norm(beta, 2);
	iter++;
	// Rcpp::Rcout << "Iter:\t" << iter << "\n   Loss:\t" << new_l << "\n   beta:" <<
	//     beta.t() << std::endl;
    }

    gradHess(beta, grad, hess, X, Y);

    arma::mat cov = -arma::inv(hess);

    boost::math::students_t dist(n - k);
    
    for (arma::uword i = 0; i < X.n_cols; i++) {
	b[i] = beta[i];
        se[i] = std::sqrt(cov(i,i));
	if (se[i]== 0.0) se[i] = 1e-10;
	z[i] = b[i] / se[i];
	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    }

    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	vNames[i] = regressors[i]->getName(); 
    }

    return CoxRegressionResult(vNames, n, b, z, pval, se, alpha, new_l);
}

CoxRegressionResult CoxRegression::regress(CensoredVariable* target,
					   std::vector<Variable*>& regressors,
					   arma::uvec& _rows) {
    int n = _rows.n_elem;
    int k = regressors.size();

    if (n < k)
	throw std::runtime_error("Cox regression ill-conditioned, samples less than regressors");

    int target_ = data.getColumn(target);

    arma::uvec regressors_ = arma::uvec(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	regressors_[i] = data.getColumn(regressors[i]);
    }


    arma::vec target_vals = dataMat.col(target_);
    target_vals = target_vals(_rows);
    arma::uvec censor = target->getCensor()(_rows);

    CensoredVariable* Y = new CensoredVariable(*target);
    Y->setCensor(target_vals, censor);

    arma::mat X;

    if (regressors.empty())
	X = arma::mat(n, 1, arma::fill::randn);
    else
	X = dataMat.submat(_rows, regressors_);

    double old_l, new_l, a, m, t, c = 0.1;
    arma::vec beta(X.n_cols, arma::fill::zeros);
    arma::vec b = beta;
    arma::vec se = arma::vec(X.n_cols);
    arma::vec z = arma::vec(X.n_cols);
    arma::vec pval = arma::vec(X.n_cols);

    if (regressors.size() == 0) {
	new_l = loss(std::move(beta), X, Y);
	b[0] = 0.0;
	z[0] = 0.0;
	pval[0] = 1.0;
	se[0] = 1.0;
	return CoxRegressionResult({""}, n, b, z, pval, se, this->alpha, new_l);
    }

    arma::mat hess(X.n_cols, X.n_cols);
    arma::vec grad(X.n_cols), p(X.n_cols);
    
    new_l = loss(std::move(beta), X, Y);

    double dbeta = 1;

    int iter = 0;

    auto start = std::chrono::high_resolution_clock::now();

    while (dbeta > 1e-5) {

	double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	if (curr > 5) {
	    throw std::runtime_error("Cox Regression not converging");
	}
	
	old_l = new_l;
	a = 1;

	gradHess(beta, grad, hess, X, Y);

	p = arma::inv(hess) * grad;
	p /= -arma::norm(p, 2);
	
	m = arma::dot(p, grad);

	if (m <= 0)
	    throw std::runtime_error("Error in Cox regression, m <= 0");

	t = -c * m;

	while (true) {
	    new_l = loss(beta + a*p, X, Y);
	    if (old_l - new_l <= a*t) {
		if (!std::isinf(new_l)) {
		    break;
		}
	    }
	    a /= 2;
	}

	beta += a*p;
	dbeta = arma::norm(a*p, 2) / arma::norm(beta, 2);
	iter++;
	// Rcpp::Rcout << "Iter:\t" << iter << "\n   Loss:\t" << new_l << "\n   beta:" <<
	//     beta.t() << std::endl;
    }

    gradHess(beta, grad, hess, X, Y);

    arma::mat cov = -arma::inv(hess);

    boost::math::students_t dist(n - k);
    
    for (arma::uword i = 0; i < X.n_cols; i++) {
	b[i] = beta[i];
        se[i] = std::sqrt(cov(i,i));
	if (se[i]== 0.0) se[i] = 1e-10;
	z[i] = b[i] / se[i];
	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    }
    
    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	vNames[i] = regressors[i]->getName(); 
    }

    return CoxRegressionResult(vNames, n, b, z, pval, se, alpha, new_l);
}



double CoxRegression::loss(arma::vec&& beta, arma::mat& X, CensoredVariable* target) {
    arma::uvec order = target->getOrder();
    arma::uvec H = target->getH();
    arma::uvec censor = target->getCensor();
    double HsumTheta, m, sub;
    double loss = 0.0;

    int n = X.n_rows;

    arma::vec logtheta = X * beta;
  
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
	    loss -= std::log(rs_sum - ((double) l) / (m * HsumTheta));
	}

	i += H[j];
	rs_sum -= sub;
    }

    return loss;
}

void CoxRegression::gradHess(arma::vec& beta, arma::vec& grad, arma::mat& hess,
			     arma::mat& X, CensoredVariable* target) {
    arma::uvec order = target->getOrder();
    arma::uvec H = target->getH();
    arma::uvec censor = target->getCensor();
    double HsumTheta, m, sub, d, phi;

    int n = X.n_rows;

    grad.fill(0);
    hess.fill(0);

    arma::vec theta = arma::exp(X * beta);
    double rs_sum = arma::accu(theta);

    arma::mat temp(beta.n_elem, beta.n_elem, arma::fill::zeros);
    arma::vec sub_num(beta.n_elem), HsumThetaVec(beta.n_elem), Z(beta.n_elem);
    arma::mat sub_outer(beta.n_elem, beta.n_elem), HsumOuter(beta.n_elem, beta.n_elem);

    arma::vec num = arma::sum(arma::diagmat(theta) * X, 0).t();
    arma::mat outer_num = X.t() * arma::diagmat(theta) * X;

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
	    temp = theta[order[i+k]] * X.row(order[i+k]).t() * X.row(order[i+k]);
	    
	    if (censor[order[i+k]]) {
		m ++;

		HsumTheta += theta[order[i+k]];
		grad += X.row(order[i+k]);
		HsumThetaVec += theta[order[i+k]] * X.row(order[i+k]);
		HsumOuter += temp;
	    }

	    sub += theta[order[i+k]];
	    sub_num += theta[order[i+k]] * X.row(order[i+k]);
	    sub_outer += temp;
	}

	if (HsumTheta - sub > 1e-5)
	    throw std::runtime_error("Error in Cox regression gradHess, HsumTheta > sub: " +
				     std::to_string(HsumTheta) + " > " + std::to_string(sub));

	if (sub - rs_sum > 1e-5) {
	    if ((H[j] + i) != n) {
		rs_sum = arma::accu(theta(order.subvec(i, n-1)));
		num = arma::sum(arma::diagmat(theta(order.subvec(i, n-1)))
				* X.rows(order.subvec(i, n-1)), 0).t();
		outer_num = (X.rows(order.subvec(i, n-1)).t()
			     * arma::diagmat(theta(order.subvec(i, n-1)))
			     * X.rows(order.subvec(i, n-1)));
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
	    grad -= Z / phi;
	    hess -= ((outer_num - d * HsumOuter) / phi - Z * Z.t() / (phi*phi));
	}

	i += H[j];
	rs_sum -= sub;
	num -= sub_num;
	outer_num -= sub_outer;
    }

}

// [[Rcpp::export]]
void CoxRegressionTest(const Rcpp::DataFrame& df, std::string targetName,
		       std::vector<std::string>& regressorNames) {
    DataSet data = DataSet(df, 5);
    CoxRegression reg(data);

    Rcpp::Rcout << "-----START----- \n";
    Variable* target = data.getVariable(targetName);
    std::vector<Variable*> regressors;
    for (std::string varName : regressorNames) {
	regressors.push_back(data.getVariable(varName));
    }
    CoxRegressionResult result = reg.regress((CensoredVariable*) target, regressors);
    Rcpp::Rcout << result;
    Rcpp::Rcout << "-----END----- \n";

    Rcpp::Rcout << "-----START----- \n";
    arma::uvec rows = arma::regspace<arma::uvec>(0,499);
    result = reg.regress((CensoredVariable*) target, regressors, rows);
    Rcpp::Rcout << result;
    Rcpp::Rcout << "-----END----- \n";

}
