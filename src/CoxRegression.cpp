// [[Rcpp::depends(BH)]]

#include "CoxRegression.hpp"
#include "CoxRegressionResult.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/normal.hpp>


CoxRegression::CoxRegression(DataSet& data) {
    this->data = data;
    dataMat = arma::mat(data.getData());
    variables = data.getVariables();
    rows = arma::uvec(data.getNumRows());
    for (arma::uword i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}


CoxRegressionResult CoxRegression::regress(const Node& target,
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

    arma::vec target_vals = dataMat.col(target_);
    target_vals = target_vals(this->rows);
    arma::uvec censor = target.getCensorVec()(this->rows);
    arma::uvec strata = target.getStrata()(this->rows);

    Node Y(target);
    Y.setCensor(target_vals, censor, strata);

    arma::mat X;

    if (regressors.empty())
	X = arma::mat(n, 1, arma::fill::randn);
    else
	X = dataMat.submat(this->rows, regressors_);

    double old_l, new_l;
    arma::vec beta(X.n_cols, arma::fill::zeros);
    arma::vec b = beta;
    arma::vec se = arma::vec(X.n_cols);
    arma::vec z = arma::vec(X.n_cols);
    arma::vec pval = arma::vec(X.n_cols);
    arma::vec res(X.n_rows, arma::fill::zeros);
    arma::vec w(X.n_rows, arma::fill::zeros);
    arma::vec eta(X.n_rows, arma::fill::zeros);

    if (regressors.size() == 0) {
	new_l = loss(beta, X, Y);
	etaGradHess(eta, res, w, Y);
	b[0] = 0.0;
	z[0] = 0.0;
	pval[0] = 1.0;
	se[0] = 1.0;
	return CoxRegressionResult({""}, n, b, z, pval, se, res, this->alpha, -new_l);
    }
    
    arma::mat hess(X.n_cols, X.n_cols);
    arma::vec grad(X.n_cols), p(X.n_cols);
    
    new_l = 1e20;

    double dbeta = 1;

    int iter = 0;

    auto start = std::chrono::high_resolution_clock::now();

    while (dbeta > 1e-8 || std::abs(new_l - old_l) > 1e-7) {

	double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	if (curr > 5) {
	    throw std::runtime_error("Cox Regression not converging");
	}
	    
	old_l = new_l;

	new_l = gradHess(beta, grad, hess, X, Y);

	// Rcpp::Rcout << "Hessian:\n" << hess << std::endl;
	// Rcpp::Rcout << "Gradient:\n" << grad << std::endl;
		
	p = arma::solve(hess, grad, arma::solve_opts::likely_sympd);

	// Rcpp::Rcout << "Update:\n" << -p << std::endl;

	beta -= p;
	dbeta = arma::norm(p, 2) / arma::norm(beta, 2);

	// p = arma::inv(hess) * grad;
	// p /= -arma::norm(p, 2);
	
	// m = arma::dot(p, grad);

	// if (m <= 0)
	//     throw std::runtime_error("Error in Cox regression, m <= 0");

	// t = -c * m;

	// while (true) {
	//     new_l = loss(beta + a*p, X, Y);
	//     if (old_l - new_l <= a*t) {
	// 	if (!std::isinf(new_l)) {
	// 	    break;
	// 	}
	//     }
	//     a /= 2;
	// }

	// beta += a*p;
	// dbeta = arma::norm(a*p, 2) / arma::norm(beta, 2);
	
	iter++;
	// Rcpp::Rcout << "Iter:\t" << iter << "\n   Loss:\t" << new_l << "\n   beta:" <<
	//     beta.t() << std::endl;
	// Rcpp::Rcout << "Iter:  " << iter << "    Loss:  " << new_l << "    |dx|/|x|:  "
	// 	    << dbeta << "    |dll|:  " << std::abs(new_l - old_l)
	// 	    << "\n    beta:  " << beta.t() << std::endl;
    }

    new_l = gradHess(beta, grad, hess, X, Y);

    arma::mat cov = arma::solve(hess, arma::eye(arma::size(hess)), arma::solve_opts::likely_sympd);

    boost::math::normal dist(0, 1);
    
    for (arma::uword i = 0; i < X.n_cols; i++) {
	b[i] = beta[i];
        se[i] = std::sqrt(cov(i,i));
	if (se[i]== 0.0) se[i] = 1e-10;
	z[i] = b[i] / se[i];
	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    }

    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	vNames[i] = regressors[i].getName(); 
    }

    eta = X * beta;
    eta -= arma::mean(eta);
    etaGradHess(eta, res, w, Y);
    arma::vec Z = eta - res / w;
    arma::vec wz = -w % Z;

    return CoxRegressionResult(vNames, n, b, z, pval, se, wz, alpha, -new_l);
}

CoxRegressionResult CoxRegression::regress(const Node& target,
					   std::vector<Node>& regressors,
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
    arma::uvec censor = target.getCensorVec()(_rows);
    arma::uvec strata = target.getStrata()(_rows);

    Node Y(target);
    Y.setCensor(target_vals, censor, strata);

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
    arma::vec res(X.n_rows, arma::fill::zeros);
    arma::vec w(X.n_rows, arma::fill::zeros);
    arma::vec eta(X.n_rows, arma::fill::zeros);

    if (regressors.size() == 0) {
	new_l = loss(beta, X, Y);
	etaGradHess(eta, res, w, Y);
	b[0] = 0.0;
	z[0] = 0.0;
	pval[0] = 1.0;
	se[0] = 1.0;
	return CoxRegressionResult({""}, n, b, z, pval, se, res, this->alpha, -new_l);
    }

    arma::mat hess(X.n_cols, X.n_cols);
    arma::vec grad(X.n_cols), p(X.n_cols);
    
    new_l = 1e20; // loss(std::move(beta), X, Y);

    double dbeta = 1;

    int iter = 0;

    auto start = std::chrono::high_resolution_clock::now();

    while (dbeta > 1e-8 || std::abs(new_l - old_l) > 1e-7) {

	double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	if (curr > 5) {
	    throw std::runtime_error("Cox Regression not converging");
	}
	
	old_l = new_l;

	new_l = gradHess(beta, grad, hess, X, Y);

	// Rcpp::Rcout << "Hessian:\n" << hess << std::endl;
	// Rcpp::Rcout << "Gradient:\n" << grad << std::endl;
		
	p = arma::solve(hess, grad, arma::solve_opts::likely_sympd);

	// Rcpp::Rcout << "Update:\n" << -p << std::endl;

	beta -= p;
	dbeta = arma::norm(p, 2) / arma::norm(beta, 2);

	// a = 1;

	// gradHess(beta, grad, hess, X, Y);

	// p = arma::inv(hess) * grad;
	// p /= -arma::norm(p, 2);
	
	// m = arma::dot(p, grad);

	// if (m <= 0)
	//     throw std::runtime_error("Error in Cox regression, m <= 0");

	// t = -c * m;

	// while (true) {
	//     new_l = loss(beta + a*p, X, Y);
	//     if (old_l - new_l <= a*t) {
	// 	if (!std::isinf(new_l)) {
	// 	    break;
	// 	}
	//     }
	//     a /= 2;
	// }

	// beta += a*p;
	// dbeta = arma::norm(a*p, 2) / arma::norm(beta, 2);
	iter++;
	// Rcpp::Rcout << "Iter:  " << iter << "    Loss:  " << new_l << "    |dx|/|x|:  "
	// 	    << dbeta << "    |dll|:  " << std::abs(new_l - old_l)
	// 	    << "\n    beta:  " << beta.t() << std::endl;
    }

    new_l = gradHess(beta, grad, hess, X, Y);

    arma::mat cov = arma::solve(hess, arma::eye(arma::size(hess)), arma::solve_opts::likely_sympd);

    boost::math::normal dist(0, 1);
    
    for (arma::uword i = 0; i < X.n_cols; i++) {
	b[i] = beta[i];
        se[i] = std::sqrt(cov(i,i));
	if (se[i]== 0.0) se[i] = 1e-10;
	z[i] = b[i] / se[i];
	pval[i] = 2 * (1.0 - boost::math::cdf(dist, std::abs(z[i])));
    }
    
    std::vector<std::string> vNames(regressors.size());

    for (int i = 0; i < regressors.size(); i++)	{
	vNames[i] = regressors[i].getName(); 
    }

    eta = X * beta;
    eta -= arma::mean(eta);
    etaGradHess(eta, res, w, Y);
    arma::vec Z = eta - res / w;
    arma::vec wz = -w % Z;

    return CoxRegressionResult(vNames, n, b, z, pval, se, wz, alpha, -new_l);
}



double CoxRegression::loss(arma::vec& beta, arma::mat& X, const Node& target) {

    double loss = 0.0;
    int numStrata = target.getNumStrata();

    arma::vec eta(X.n_rows);
    arma::vec theta(X.n_rows);
    double rs_sum = 0.0;
        
    for (int strat = 0; strat < numStrata; strat++) {
	arma::uvec index = target.getIndex(strat);
	arma::uvec order = target.getOrder(strat);
	arma::uvec H = target.getH(strat);
	arma::uvec censor = target.getCensor(strat);
	
	double HsumTheta, m, sub, d, phi;

	int n = index.n_elem;

	// initiallize risk set sum
	rs_sum = 0.0;

	
	for (uint rowIdx = 0; rowIdx < n; rowIdx++) {
	    uint xRowIdx = index[rowIdx];
	    eta[xRowIdx] = arma::dot(X.row(xRowIdx), beta);
	    theta[xRowIdx] = std::exp(eta[xRowIdx]);
	    rs_sum += theta[xRowIdx];
	}


	// Calculate loss for this strata
	int i = 0;
	for (int j = 0; j < H.n_elem; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;

	    for (int k = 0; k < H[j]; k++) {
		uint rowIdx = order[i+k];
		uint xRowIdx = index[rowIdx];
		uint delta = censor[rowIdx];

		m += delta;
		loss -= delta * eta[xRowIdx];
		HsumTheta += delta * theta[xRowIdx];
		sub += theta[xRowIdx];	        
	    }

	    if (HsumTheta - sub > 1e-5)
		throw std::runtime_error("Error in Cox regression gradHess, HsumTheta > sub: " +
					 std::to_string(HsumTheta) + " > " + std::to_string(sub));

	    if (sub - rs_sum > 1e-5) {
		if ((H[j] + i) != n) {
		    arma::uvec remainingXIdx = index(order.subvec(i, n-1));
		    arma::vec remainingTheta = theta(remainingXIdx);
		    rs_sum = arma::accu(remainingTheta);
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Cox regression gradHess, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }
	
	    for (int l = 0; l < m; l++) {
		loss += std::log(rs_sum - l / ((double) m) * HsumTheta);
	    }

	    i += H[j];
	    rs_sum -= sub;
	}
    }
    
    return loss;
}

double CoxRegression::gradHess(arma::vec& beta, arma::vec& grad, arma::mat& hess,
			     arma::mat& X, const Node& target) {

    double loss = 0.0;
    int numStrata = target.getNumStrata();

    arma::vec eta(X.n_rows);
    arma::vec theta(X.n_rows);
    double rs_sum = 0.0;

    // arma::mat temp(beta.n_elem, beta.n_elem, arma::fill::zeros);
    arma::vec sub_num(beta.n_elem), HsumThetaVec(beta.n_elem), Z(beta.n_elem);
    arma::mat sub_outer(beta.n_elem, beta.n_elem), HsumOuter(beta.n_elem, beta.n_elem);
    arma::vec num(beta.n_elem);
    arma::mat outer_num(beta.n_elem, beta.n_elem);

    for (uint ii = 0; ii < X.n_cols; ii++) {
	grad[ii] = 0;
	for (uint jj = ii; jj < X.n_cols; jj++) {
	    hess(ii,jj) = 0;
	}
    }
        
    for (int strat = 0; strat < numStrata; strat++) {
	arma::uvec index = target.getIndex(strat);
	arma::uvec order = target.getOrder(strat);
	arma::uvec H = target.getH(strat);
	arma::uvec censor = target.getCensor(strat);
	
	double HsumTheta, m, sub, d, phi;

	int n = index.n_elem;

	// initiallize risk set sum, numerator sum, and outer product
	// numerator sum for this strata
	rs_sum = 0.0;

	for (uint ii = 0; ii < X.n_cols; ii++) {
	    num[ii] = 0;
	    for (uint jj = ii; jj < X.n_cols; jj++) {
		outer_num(ii,jj) = 0;
	    }
	}

	for (uint rowIdx = 0; rowIdx < n; rowIdx++) {
	    uint xRowIdx = index[rowIdx];
	    eta[xRowIdx] = arma::dot(X.row(xRowIdx), beta);
	    theta[xRowIdx] = std::exp(eta[xRowIdx]);
	    rs_sum += theta[xRowIdx];
	    for (uint ii = 0; ii < X.n_cols; ii++) {
		double wxii = theta[xRowIdx] * X(xRowIdx, ii);
		num[ii] += wxii;
		for (uint jj = ii; jj < X.n_cols; jj++) {
		    outer_num(ii, jj) += wxii * X(xRowIdx, jj);
		}
	    }
	}


	// Calculate loss, gradient, and hessian for this strata
	int i = 0;
	for (int j = 0; j < H.n_elem; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;

	    for (uint ii = 0; ii < X.n_cols; ii++) {
		HsumThetaVec[ii] = 0;
		sub_num[ii] = 0;
		for (uint jj = ii; jj < X.n_cols; jj++) {
		    HsumOuter(ii,jj) = 0;
		    sub_outer(ii,jj) = 0;
		}
	    }

	    for (int k = 0; k < H[j]; k++) {
		uint rowIdx = order[i+k];
		uint xRowIdx = index[rowIdx];
		uint delta = censor[rowIdx];

		m += delta;
		loss -= delta * eta[xRowIdx];
		HsumTheta += delta * theta[xRowIdx];
		sub += theta[xRowIdx];

		for (uint ii = 0; ii < X.n_cols; ii++) {
		    double xii = X(xRowIdx, ii);
		    grad[ii] -= delta * xii;
		    double wxii = theta[xRowIdx] * xii;
		    HsumThetaVec[ii] += delta * wxii;
		    sub_num[ii] += wxii;
		    for (uint jj = ii; jj < X.n_cols; jj++) {
			double xjj = X(xRowIdx, jj);
			HsumOuter(ii,jj) += delta * wxii * xjj;
			sub_outer(ii,jj) += wxii * xjj;
		    }
		}
	        
	    }

	    if (HsumTheta - sub > 1e-5)
		throw std::runtime_error("Error in Cox regression gradHess, HsumTheta > sub: " +
					 std::to_string(HsumTheta) + " > " + std::to_string(sub));

	    if (sub - rs_sum > 1e-5) {
		if ((H[j] + i) != n) {
		    arma::uvec remainingXIdx = index(order.subvec(i, n-1));
		    arma::vec remainingTheta = theta(remainingXIdx);
		    rs_sum = arma::accu(remainingTheta);
		    num = arma::sum(arma::diagmat(remainingTheta)
				    * X.rows(remainingXIdx), 0).t();
		    outer_num = (X.rows(remainingXIdx).t()
				 * arma::diagmat(remainingTheta)
				 * X.rows(remainingXIdx));
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Cox regression gradHess, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }
	
	    for (int l = 0; l < m; l++) {
		d = l / ((double) m);
		phi = rs_sum - d * HsumTheta;
		loss += std::log(phi);

		Z = num - d * HsumThetaVec;

		for (uint ii = 0; ii < X.n_cols; ii++) {
		    grad[ii] += Z[ii] / phi;
		    for (uint jj = ii; jj < X.n_cols; jj++) {
			hess(ii, jj) += (outer_num(ii,jj) - d * HsumOuter(ii, jj)) / phi - Z[ii] * Z[jj] / (phi * phi);
		    }
		}
		
	    }

	    i += H[j];
	    rs_sum -= sub;
	    num -= sub_num;
	    outer_num -= sub_outer;
	}
    }

    // set the lower triangle of the hessian equal to the upper
    // triangle

    for (uint ii = 0; ii < X.n_cols; ii++) {
	for (uint jj = ii+1; jj < X.n_cols; jj++) {
	    hess(jj, ii) = hess(ii, jj);
	}
    }
    
    return loss;
}


double CoxRegression::etaGradHess(arma::vec& eta, arma::vec& grad, arma::vec& diagHess,
				  const Node& target) {

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

    diagHess.replace(0, -1e-10);

    return loss;
}


// no export // [[Rcpp::export]]
void CoxRegressionTest(const Rcpp::DataFrame& df, std::string targetName,
		       std::vector<std::string>& regressorNames, int repetitions) {
    DataSet data = DataSet(df);
    
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

    CoxRegression reg(data);

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
