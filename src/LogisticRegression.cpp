// [[Rcpp::depends(BH, RcppThread)]]

#include "LogisticRegression.hpp"
#include "LogisticRegressionResult.hpp"
// #include "Variable.hpp"
// #include "DiscreteVariable.hpp"
#include "Node.hpp"
#include "RcppThread.h"
#include <iostream>
#include <algorithm>
#include <math.h>
// #include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <fstream>
#include <thread>

LogisticRegression::LogisticRegression(DataSet &data) {
    this->data = data;
    this->dataCols = data.getData().t();
    this->rows = arma::uvec(data.getNumRows());
    for (int i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}


LogisticRegressionResult LogisticRegression::regress(const Node& x,
						     std::vector<Node>& regressors)
{
    if (!binary(x)) {
        throw std::runtime_error("Target must be binary.");
    }

    for (const Node& var : regressors) {
	// if (var.isCensored()) RcppThread::Rcout << "LogisticRegression::Censored Regressor " + var.getName() + "\n";
        // if (var.isCensored()) continue;
        if (!var.isContinuous() && !binary(var)) {
	    // RcppThread::Rcout << "LogisticRegression::Regressors must be continuous, binary, or censored.\n";
            throw std::runtime_error("Regressors must be continuous or binary.");
        }
    }

    // std::ofstream logfile;
    // logfile.open("log_reg_debug.log", std::ios_base::app);

    arma::mat regressors_ = arma::mat(regressors.size(), rows.size());

    for (arma::uword j = 0; j < regressors.size(); j++)
    {
	// logfile << regressors[j]->getName() + '\t';
        int col = data.getColumn(regressors[j]);
	// logfile << col + '\t';

        for (arma::uword i = 0; i < rows.size(); i++) {
            regressors_(j, i) = dataCols(col, rows[i]);
        }
    }

    // for (arma::uword j = 0; j < regressors.size(); j++) {
    // 	if (regressors[j].isCensored()) {
    // 	    // RcppThread::Rcout << "LogisticRegression::Censored Regressor " + regressors[j].getName() + "\n";
    // 	    regressors_.row(j) = getWZ(regressors[j], x).elem(rows).t();
    // 	}
    // }

    // logfile << "\n";

    arma::uvec target = arma::uvec(rows.size());
    // logfile << x->getName() << std::endl;
    int col = data.getColumn(x);
    // logfile << col + '\n';

    for (arma::uword k = 0; k < rows.size(); k++)
    {
        target[k] = data.getInt(rows[k], col);
    }

    std::vector<std::string> regressorNames(regressors.size());

    for (int l = 0; l < regressors.size(); l++)
    {
        regressorNames[l] = regressors[l].getName();
    }

    // logfile.close();

    return regress(target, x.getName(), regressors_, regressorNames);
}


LogisticRegressionResult LogisticRegression::regress(const Node& x,
						     std::vector<Node>& regressors,
						     arma::uvec _rows)
{
    if (!binary(x)) {
        throw std::runtime_error("Target must be binary.");
    }

    for (const Node& var : regressors) {
	// if (var.isCensored()) RcppThread::Rcout << "LogisticRegression::Censored Regressor " + var.getName() + "\n";
	// if (var.isCensored()) continue;
        if (!var.isContinuous() && !binary(var)) {
	    // RcppThread::Rcout << "LogisticRegression::Regressors must be continuous, binary, or censored.\n";
            throw std::runtime_error("Regressors must be continuous or binary.");
        }
    }

    arma::mat regressors_ = arma::mat(regressors.size(), _rows.size());

    for (arma::uword j = 0; j < regressors.size(); j++)
    {
        int col = data.getColumn(regressors[j]);

        for (arma::uword i = 0; i < _rows.size(); i++)
        {
            regressors_(j, i) = dataCols(col, _rows[i]);
        }
    }

    // for (arma::uword j = 0; j < regressors.size(); j++) {
    // 	// RcppThread::Rcout << (regressors[j].isCensored() ? "Censored: " : "Not Censored: ") + regressors[j].getName() + "\n";
    // 	if (regressors.at(j).isCensored()) {
    // 	    // RcppThread::Rcout << "LogisticRegression::Censored Regressor " + regressors[j].getName() + "\n";
    // 	    regressors_.row(j) = getWZ(regressors[j], x).elem(_rows).t();
    // 	}
    // }

    arma::uvec target = arma::uvec(_rows.size());
    int col = data.getColumn(x);

    for (arma::uword k = 0; k < _rows.size(); k++)
    {
        target[k] = data.getInt(_rows[k], col);
    }

    std::vector<std::string> regressorNames(regressors.size());

    for (int l = 0; l < regressors.size(); l++)
    {
        regressorNames[l] = regressors[l].getName();
    }

    return regress(target, x.getName(), regressors_, regressorNames);
}


LogisticRegressionResult LogisticRegression::regress(arma::uvec& target,
						     std::string targetName,
						     arma::mat& regressors,
						     std::vector<std::string>& regressorNames)
{
    arma::mat x;
    int numRegressors = regressors.n_rows;
    int numCases = target.size();

    // RcppThread::Rcout << "Num regressors: " << numRegressors << std::endl;

    // RcppThread::Rcout << "regressors mat:\n" << regressors << std::endl;

    // make a new matrix x with all the columns of regressorsh
    // but with first column all 1.0's.
    x = arma::ones(numRegressors + 1, numCases);

    // copy numRegressors of values from regressors into x
    if (numRegressors > 0) {
	x.submat(1, 0, size(regressors)) = regressors;
    }
    
    // RcppThread::Rcout << "x mat:\n" << x << std::endl;

    arma::vec xMeans = arma::zeros(numRegressors + 1);
    arma::vec xStdDevs = arma::zeros(numRegressors + 1);

    arma::vec y0 = arma::vec(numCases, arma::fill::zeros);
    arma::vec y1 = arma::vec(numCases, arma::fill::zeros);

    // std::ofstream logfile;

    // logfile.open("../debug.log", std::ios_base::app);

    for (arma::uword i = 0; i < numCases; i++)
    {
        y0[i] = 0;
        y1[i] = 0;
    }

    // RcppThread::Rcout << "target: " << target.subvec(0,5).t();
    
    int ny0 = 0;
    int ny1 = 0;
    int nc = 0;

    for (arma::uword i = 0; i < numCases; i++)
    {
        if (target[i] == 0)
        {
            y0[i] = 1;
            ny0++;
        }
        else
        {
            y1[i] = 1;
            ny1++;
        }
        nc += y0[i] + y1[i];
        for (arma::uword j = 1; j <= numRegressors; j++)
        {
            xMeans[j] += (y0[i] + y1[i]) * x(j, i);
            xStdDevs[j] += (y0[i] + y1[i]) * x(j, i) * x(j, i);
        }
    }

    for (arma::uword j = 1; j <= numRegressors; j++)
    {
        xMeans[j] /= nc;
        xStdDevs[j] /= nc;
        xStdDevs[j] = std::sqrt(std::abs(xStdDevs[j] - xMeans[j] * xMeans[j]));
    }

    xMeans[0] = 0.0;
    xStdDevs[0] = 1.0;

    for (arma::uword i = 0; i < nc; i++)
    {
        for (arma::uword j = 1; j <= numRegressors; j++)
        {
            x(j, i) = (x(j, i) - xMeans[j]) / xStdDevs[j];
        }
    }

    // arma::mat xTx = x * x.t(); // numRegressors + 1 x numRegressors + 1 

    arma::vec par = arma::vec(numRegressors + 1, arma::fill::zeros);
    arma::vec parStdErr = arma::vec(numRegressors + 1, arma::fill::zeros);
    arma::mat arr = arma::mat(numRegressors + 1, numRegressors + 2, arma::fill::zeros);
    arma::vec coefficients;
    std::vector<std::string> sigMarker;
    sigMarker.reserve(numRegressors);
    arma::vec pValues = arma::vec(numRegressors + 1, arma::fill::zeros);
    arma::vec zScores = arma::vec(numRegressors + 1, arma::fill::zeros);
    double lnV;
    double ln1mV;
    double llP = 2e+10;
    double ll = 1e+10;
    double llN = 0.0;
    double lam = 0.888;
    double chiSq;

    // auto start = std::chrono::high_resolution_clock::now();
    
    par[0] = std::log((double)ny1 / (double)ny0);
    for (arma::uword j = 1; j <= numRegressors; j++) {
	par[j] = 0.0;
    }

    llP = 2e+10;
    ll = 1e+10;
    llN = 0.0;

    int iter = 0;

    auto start = std::chrono::high_resolution_clock::now();

    while (std::abs(llP - ll) > 1e-7) {
	double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	iter++;
	if (curr > 5 || iter > 100 || !par.is_finite()) {
	    // RcppThread::Rcout << "Time out : Logistic Regression not converging" << std::endl
	    //  		      << "Loglikelihood: " << ll << std::endl;
	    if (!par.is_finite()) {
		std::stringstream ss;
		ss << "Logistic Regression not converging: Non-finite values in coefficient vector" << std::endl;
		ss << "Regressing " << targetName << " on { ";
		for (std::string varName : regressorNames) {
		    ss << varName << " ";
		}
		ss << "}\n";
		ss << "   Iter: " << iter << std::endl
		   << "   Loglikelihood: " << ll << std::endl
		   << "   |dx|: " << arma::mean(arma::abs(arr.col(numRegressors + 1))) << std::endl
		   << "   updates: " << arr.col(numRegressors + 1) << std::endl
		   << "   Parameters = " << par.t() << std::endl;
		throw std::runtime_error(ss.str());
	    }
	    // RcppThread::Rcout << "Time out : Logistic Regression not converging" << std::endl
	    // 		      << "   Iter: " << iter << std::endl
	    // 		      << "   Loglikelihood: " << ll << std::endl
	    // 		      << "   Null Loglikelihood: " << llN << std::endl;
	    // throw std::runtime_error("Logistic Regression not converging");
	    break;
	}

	llP = ll;
	ll = 0.0;

	for (arma::uword j = 0; j <= numRegressors; j++) {
	    for (arma::uword k = j; k <= numRegressors + 1; k++) {
		arr(j, k) = 0.0;
	    }
	}

	for (arma::uword i = 0; i < nc; i++) {
	    double q;
	    double v = par[0];

	    for (arma::uword j = 1; j <= numRegressors; j++) {
		v += par[j] * x(j, i);
	    }

	    if (v > 15.0) {
		lnV = -std::exp(-v);
		ln1mV = -v;
		q = std::exp(-v);
		v = std::exp(lnV);
	    } else {
		if (v < -15.0) {
		    lnV = v;
		    ln1mV = -std::exp(v);
		    q = std::exp(v);
		    v = std::exp(lnV);
		} else {
		    v = 1.0 / (1 + std::exp(-v));
		    lnV = std::log(v);
		    ln1mV = std::log(1.0 - v);
		    q = v * (1.0 - v);
		}
	    }

	    ll -= (2.0 * y1[i] * lnV + 2.0 * y0[i] * ln1mV);

	    for (arma::uword j = 0; j <= numRegressors; j++) {
		double xij = x(j, i);
		arr(j, numRegressors + 1) += xij * (y1[i] * (1.0 - v) + y0[i] * (-v));

		for (arma::uword k = j; k <= numRegressors; k++) {
		    arr(j, k) += xij * x(k, i) * q * (y0[i] + y1[i]);
		}
	    }
	}

	// for (arma::uword j = 1; j <= numRegressors; j++) {
	//     ll += lam / 2 * std::pow(par[j], 2);
	//     arr(j, numRegressors + 1) += lam * par[j];
	//     arr(j, j) += lam;
	// }

	if (llP == 1e+10) {
	    llN = ll;
	}

	// RcppThread::Rcout << arr << std::endl;

	for (arma::uword j = 1; j <= numRegressors; j++) {
	    for (arma::uword k = 0; k < j; k++) {
		arr(j, k) = arr(k, j);
	    }
	}

	arr.col(numRegressors + 1) = arma::solve(arr(arma::span(0,numRegressors),
						     arma::span(0,numRegressors)),
						 arr.col(numRegressors + 1),
						 arma::solve_opts::likely_sympd);

	// for (arma::uword i = 0; i <= numRegressors; i++) {
	//     double s = arr(i, i);
	//     arr(i, i) = 1.0;
	//     for (arma::uword k = 0; k <= numRegressors + 1; k++) {
	// 	arr(i, k) = arr(i, k) / s;
	//     }

	//     for (arma::uword j = 0; j <= numRegressors; j++) {
	// 	if (i != j) {
	// 	    s = arr(j, i);
	// 	    arr(j, i) = 0.0;
	// 	    for (arma::uword k = 0; k <= numRegressors + 1; k++) {
	// 		arr(j, k) = arr(j, k) - s * arr(i, k);
	// 	    }
	// 	}
	//     }
	// }

	// RcppThread::Rcout << arr << std::endl;

	for (arma::uword j = 0; j <= numRegressors; j++) {
	    par[j] += arr(j, numRegressors + 1);
	}

	// RcppThread::Rcout << "    Iter " << iter << ":    ll = " << ll << "    |dx| = " << arma::mean(arma::abs(arr.col(numRegressors + 1))) << "\n";
    }

    // if (ll > llN) ll = llN;
    
    chiSq = llN - ll;
    
    // if (chiSq < 0 || std::isnan(chiSq)) {
    // 	// RcppThread::Rcout << "Logistic Regression not converging" << std::endl
    // 	// 		  << "Loglikelihood: " << ll << std::endl
    // 	// 		  << "ChiSq: " << chiSq << std::endl;
    // 	// throw std::runtime_error("Logistic Regression did not converge");
    // }
    
    double zScore;

    // arma::mat cov(numRegressors+1, numRegressors+1, arma::fill::zeros);

    arma::mat cov = arma::solve(arr(arma::span(0,numRegressors),
				    arma::span(0,numRegressors)),
				arma::eye(numRegressors+1, numRegressors+1),
				arma::solve_opts::likely_sympd);

    // bool success = arma::inv_sympd(cov,
    // 				   arr(arma::span(0,numRegressors),
    // 				       arma::span(0,numRegressors)),
    // 				   arma::inv_opts::allow_approx);

    // if (!success) {
    // 	RcppThread::Rcout << "arr:\n" << arr << "\ncov:\n" << cov << "\npar:\n" << par.t() << std::endl;
    // }

    for (arma::uword j = 1; j <= numRegressors; j++) {
	par[j] = par[j] / xStdDevs[j];
	parStdErr[j] = std::sqrt(cov(j, j)) / xStdDevs[j];
	par[0] = par[0] - par[j] * xMeans[j];
	zScore = par[j] / parStdErr[j];

	if (std::isnan(zScore)) {
	    // RcppThread::Rcout << "Logistic Regression not converging" << std::endl
	    // 		      << "Loglikelihood: " << ll << std::endl
	    // 		      << "NaN zScore" << std::endl;
	    throw std::runtime_error("Logistic Regression not converging, NaN coefficient");
	}

	// double prob = norm(std::abs(zScore));
	pValues[j] = norm(std::abs(zScore));
	zScores[j] = zScore;
    }

    parStdErr[0] = std::sqrt(arr(0, 0));
    zScore = par[0] / parStdErr[0];

    if (std::isnan(zScore)) {
    	// RcppThread::Rcout << "Logistic Regression not converging" << std::endl
    	// 		  << "Loglikelihood: " << ll << std::endl
    	// 		  << "NaN zScore" << std::endl;;
    	throw std::runtime_error("Logistic Regression not converging, NaN intercept");
    }

    // RcppThread::Rcout << "Logistic Regression converged" << std::endl
    // 		      << "   Iter: " << iter << std::endl
    // 		      << "   Loglikelihood: " << ll << std::endl
    // 		      << "   Null Loglikelihood: " << llN << std::endl;

    pValues[0] = norm(zScore);
    zScores[0] = zScore;

    double intercept = par[0];
    coefficients = par;

    return LogisticRegressionResult(targetName, regressorNames, xMeans, xStdDevs,
				    numRegressors, ny0, ny1, coefficients, parStdErr,
				    pValues, intercept, ll, sigMarker, chiSq, alpha);
}

double LogisticRegression::norm(double z) {
    // std::ofstream logfile;
    // logfile.open("../debug.log", std::ios_base::app);

    double q = z * z;
    // const double pi = M_PI;
    double piOver2 = M_PI / 2.0;

    // Rcpp::Rcout << "chisq = " << q << std::endl;
    
    // logfile << "LOGISTIC REGRESSION" << std::endl;
    // logfile << "q = " << q << std::endl << std::endl;
    // logfile.close();

    if (std::abs(q) > 7.0)
    {
        return (1.0 - 1.0 / q + 3.0 / (q * q)) * std::exp(-q / 2.0) /
               (std::abs(z) * std::sqrt(piOver2));
    }
    else
    {
        boost::math::chi_squared dist(1);

        double p = cdf(dist, q);
        return (p);
    }
}

// arma::vec LogisticRegression::getWZ(Node coxNode, Node target) {
//     if (!coxNode.isCensored())
// 	throw std::runtime_error("Node " + coxNode.getName() + " is not censored");

//     coxNode = data.getVariable(coxNode.getName());
    
//     // std::vector<std::string> _neighbors(coxNode.getNeighbors());
//     // std::vector<Node> neighbors;

//     // for (const std::string& name : _neighbors) {
//     // 	neighbors.push_back(data.getVariable(name));
//     // }

//     // Node targetNode;
//     // auto multinomLoc = target.getName().find("MULTINOM");
//     // if (multinomLoc != std::string::npos) {
//     // 	std::string varName = target.getName().substr(0, multinomLoc);
//     // 	targetNode = data.getVariable(varName);
//     // } else {
//     // 	targetNode = target;
//     // }
    
//     // auto targetLoc = std::find(neighbors.begin(), neighbors.end(), targetNode);
    
//     // if (targetLoc != neighbors.end()) {
//     // 	// RcppThread::Rcout << "LogisticRegression::getWZ:: Erasing " + targetLoc->getName() + "\n";
//     // 	// auto multinomLoc = target.getName().find("MULTINOM"); 
//     // 	// if (multinomLoc != std::string::npos) {
//     // 	//     std::string varName = target.getName().substr(0, multinomLoc);
//     // 	//     for (auto it = neighbors.begin(); it != neighbors.end(); it++) {
//     // 	// 	if (it->getName().find(varName) != std::string::npos) {
//     // 	// 	    // RcppThread::Rcout << "LogisticRegression::getWZ:: Erasing " + it->getName() + "\n";
//     // 	// 	    it = neighbors.erase(it);
//     // 	// 	    it--;
//     // 	// 	}
//     // 	//     }
//     // 	// } else {
//     // 	//     // RcppThread::Rcout << "LogisticRegression::getWZ:: Erasing " + targetLoc->getName() + "\n";
//     // 	//     neighbors.erase(targetLoc);
//     // 	// }
//     // 	// CoxRegressionResult result = coxRegression.regress(coxNode, neighbors);
//     // 	// RcppThread::Rcout << "WZ:\n" << result.getResid().t() << std::endl;
//     // 	return WZmap.at(std::minmax(coxNode, targetNode));
//     // }
    
//     // RcppThread::Rcout << "WZ:\n" << coxNode.getWZ().t() << std::endl;
//     return coxNode.getWZ();
// }

// void LogisticRegressionTest(const Rcpp::DataFrame& df) {
//     Rcpp::Rcout << "Start** \n";
//     DataSet data = DataSet(df, 5);
//     Rcpp::Rcout << "Dataset Constructed \n";
//     LogisticRegression reg(data);
//     Rcpp::Rcout << "Finished Constructing \n";

//     for (int i = 0; i < data.getNumColumns(); i++) {
// 	Rcpp::Rcout << "-----START----- \n";
// 	if(data.getVariable(i)->isDiscrete())
// 	    {
// 		const Node& x = (Node)data.getVariable(i);
// 		Rcpp::Rcout << x->getName() << std::endl;
// 		std::vector<Node> regressors(data.getVariables());

// 		regressors.erase(regressors.begin() + i);

// 		std::string xName = x->getName();
// 		xName = xName.substr(0, xName.find('.'));

// 		for (auto it = regressors.begin(); it != regressors.end(); it++) {
// 		    std::string tempName = (*it)->getName();
// 		    tempName = tempName.substr(0, tempName.find('.'));
// 		    if (tempName == xName) {
// 			it = regressors.erase(it);
// 			it--;
// 		    }
// 		}

// 		LogisticRegressionResult* result = reg.regress(x, regressors);
// 		Rcpp::Rcout << *result;
// 	    }
// 	Rcpp::Rcout << "-----END----- \n";
// 	Rcpp::Rcout << "\n";
// 	Rcpp::Rcout << "\n";
// 	Rcpp::Rcout << "\n";
// 	Rcpp::Rcout << "\n";
// 	Rcpp::Rcout << "\n";
//     }
// }
