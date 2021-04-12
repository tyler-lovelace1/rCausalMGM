// [[Rcpp::depends(BH)]]

#include "LogisticRegression.hpp"
#include "LogisticRegressionResult.hpp"
#include "Variable.hpp"
#include "DiscreteVariable.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <fstream>
#include <thread>

LogisticRegression::LogisticRegression(DataSet &data)
{
    this->data = data;
    this->dataCols = data.getData().t();
    this->rows = arma::uvec(data.getNumRows());
    for (int i = 0; i < data.getNumRows(); i++)
        rows[i] = i;
}

LogisticRegression::LogisticRegression(LogisticRegression &lr)
{
    this->data = lr.data;
    this->dataCols = this->data.getData().t();
    this->rows = arma::uvec(this->data.getNumRows());
    for (int i = 0; i < this->data.getNumRows(); i++)
        rows[i] = i;
}

LogisticRegression::LogisticRegression(LogisticRegression &&lr)
{
    this->data = lr.data;
    this->dataCols = this->data.getData().t();
    this->rows = arma::uvec(this->data.getNumRows());
    for (int i = 0; i < this->data.getNumRows(); i++)
        rows[i] = i;
}

LogisticRegression &LogisticRegression::operator=(LogisticRegression &lr)
{
    this->data = lr.data;
    this->dataCols = this->data.getData().t();
    this->rows = arma::uvec(this->data.getNumRows());
    for (int i = 0; i < this->data.getNumRows(); i++)
        rows[i] = i;
    return *this;
}

LogisticRegression &LogisticRegression::operator=(LogisticRegression &&lr)
{
    this->data = lr.data;
    this->dataCols = this->data.getData().t();
    this->rows = arma::uvec(this->data.getNumRows());
    for (int i = 0; i < this->data.getNumRows(); i++)
        rows[i] = i;
    return *this;
}

LogisticRegressionResult LogisticRegression::regress(DiscreteVariable *x, std::vector<Variable*>& regressors)
{
    // if (!binary(x)) {
    //     throw new IllegalArgumentException("Target must be binary.");
    // }

    // for (Variable var : regressors) {
    //     if (!(node instanceof ContinuousVariable || binary(node))) {
    //         throw new IllegalArgumentException("Regressors must be continuous or binary.");
    //     }
    // }

    arma::mat regressors_ = arma::mat(regressors.size(), rows.size());

    for (arma::uword j = 0; j < regressors.size(); j++)
    {
        int col = data.getColumn(regressors[j]);

        for (arma::uword i = 0; i < getRows().size(); i++)
        {
            regressors_(j, i) = dataCols(col, rows[i]);
        }
    }

    arma::uvec target = arma::uvec(rows.size());
    int col = data.getColumn(data.getVariable(x->getName()));

    for (arma::uword k = 0; k < rows.size(); k++)
    {
        target[k] = data.getInt(rows[k], col);
    }

    std::vector<std::string> regressorNames(regressors.size());

    for (int l = 0; l < regressors.size(); l++)
    {
        Variable *var = regressors[l];
        regressorNames[l] = var->getName();
    }

    return regress(target, x->getName(), regressors_, regressorNames);
}


LogisticRegressionResult LogisticRegression::regress(DiscreteVariable *x, std::vector<Variable*>& regressors, arma::uvec _rows)
{
    // if (!binary(x)) {
    //     throw new IllegalArgumentException("Target must be binary.");
    // }

    // for (Variable var : regressors) {
    //     if (!(node instanceof ContinuousVariable || binary(node))) {
    //         throw new IllegalArgumentException("Regressors must be continuous or binary.");
    //     }
    // }

    arma::mat regressors_ = arma::mat(regressors.size(), _rows.size());

    for (arma::uword j = 0; j < regressors.size(); j++)
    {
        int col = data.getColumn(regressors[j]);

        for (arma::uword i = 0; i < _rows.size(); i++)
        {
            regressors_(j, i) = dataCols(col, _rows[i]);
        }
    }

    arma::uvec target = arma::uvec(_rows.size());
    int col = data.getColumn(data.getVariable(x->getName()));

    for (arma::uword k = 0; k < _rows.size(); k++)
    {
        target[k] = data.getInt(_rows[k], col);
    }

    std::vector<std::string> regressorNames(regressors.size());

    for (int l = 0; l < regressors.size(); l++)
    {
        Variable *var = regressors[l];
        regressorNames[l] = var->getName();
    }

    return regress(target, x->getName(), regressors_, regressorNames);
}


LogisticRegressionResult LogisticRegression::regress(arma::uvec& target,
                                                      std::string targetName, arma::mat& regressors,
                                                      std::vector<std::string>& regressorNames)
{
    arma::mat x;
    int numRegressors = regressors.n_rows;
    int numCases = target.size();

    // make a new matrix x with all the columns of regressorsh
    // but with first column all 1.0's.
    x = arma::ones(numRegressors + 1, numCases);

    // copy numRegressors of values from regressors into x
    x.submat(1, 0, size(regressors)) = regressors;

    arma::vec xMeans = arma::zeros(numRegressors + 1);
    arma::vec xStdDevs = arma::zeros(numRegressors + 1);

    arma::vec y0 = arma::vec(numCases, arma::fill::zeros);
    arma::vec y1 = arma::vec(numCases, arma::fill::zeros);

    std::ofstream logfile;

    // logfile.open("../debug.log", std::ios_base::app);

    for (arma::uword i = 0; i < numCases; i++)
    {
        y0[i] = 0;
        y1[i] = 0;
    }

    int ny0 = 0;
    int ny1 = 0;
    int nc = 0;

    for (arma::uword i = 0; i < numCases; i++)
    {
        if (target[i] == 0.0)
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

    arma::vec par = arma::vec(numRegressors + 1, arma::fill::zeros);
    arma::vec parStdErr = arma::vec(numRegressors + 1, arma::fill::zeros);
    arma::mat arr = arma::mat(numRegressors + 1, numRegressors + 2, arma::fill::zeros);
    arma::vec coefficients;
    std::vector<std::string> sigMarker;
    sigMarker.reserve(numRegressors);
    arma::vec pValues = arma::vec(numRegressors + 1, arma::fill::ones);
    arma::vec zScores = arma::vec(numRegressors + 1, arma::fill::zeros);
    double lnV;
    double ln1mV;
    double llP = 2e+10;
    double ll = 1e+10;
    double llN = 0.0;
    double lam = 0.1;
    double chiSq;

    while (true)
    {

        par[0] = std::log((double)ny1 / (double)ny0);
        for (arma::uword j = 1; j <= numRegressors; j++)
        {
            par[j] = 0.0;
        }

        llP = 2e+10;
        ll = 1e+10;
        llN = 0.0;

        int iter = 0;

        auto start = std::chrono::high_resolution_clock::now();
        while (std::abs(llP - ll) > 1e-7)
        {
            double curr = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
            if (curr > 5)
            {
                // logfile << "timed out" << std::endl;
                // for (arma::uword j = 0; j <= numRegressors; j++) {
                //     par[j] += std::numeric_limits<double>::quiet_NaN();
                // }
                // break;
                throw new std::runtime_error("Log Reg taking forever");
            }

            llP = ll;
            ll = 0.0;

            for (arma::uword j = 0; j <= numRegressors; j++)
            {
                for (arma::uword k = j; k <= numRegressors + 1; k++)
                {
                    arr(j, k) = 0.0;
                }
            }

            for (arma::uword i = 0; i < nc; i++)
            {
                double q;
                double v = par[0];

                for (arma::uword j = 1; j <= numRegressors; j++)
                {
                    v += par[j] * x(j, i);
                }

                // if (numRegressors == 3 && regressorNames[0] == "X311MULTINOM.1" &&
                //     regressorNames[1] == "X311MULTINOM.2" && regressorNames[2] == "X311MULTINOM.3") {
                //     logfile << std::this_thread::get_id() << "\titer = " << iter << std::endl;
                //     logfile << std::this_thread::get_id() << "\t|v| > 15 = "
                // 	     << (std::abs(v) > 15) << std::endl;
                //     logfile << std::this_thread::get_id() << "\tv = " << v << std::endl;
                // }

                if (v > 15.0)
                {
                    // v = v < expBound ? v : expBound;
                    lnV = -std::exp(-v);
                    ln1mV = -v;
                    q = std::exp(-v);
                    v = std::exp(lnV);
                }
                else
                {
                    if (v < -15.0)
                    {
                        // v = v > -expBound ? v : -expBound;
                        lnV = v;
                        ln1mV = -std::exp(v);
                        q = std::exp(v);
                        v = std::exp(lnV);
                    }
                    else
                    {
                        v = 1.0 / (1 + std::exp(-v));
                        lnV = std::log(v);
                        ln1mV = std::log(1.0 - v);
                        q = v * (1.0 - v);
                    }
                }

                ll -= 2.0 * y1[i] * lnV - 2.0 * y0[i] * ln1mV;

                // if (numRegressors == 3 && regressorNames[0] == "X3MULTINOM.1" &&
                //     regressorNames[1] == "X3MULTINOM.2" &&regressorNames[2] == "X3MULTINOM.3") {
                //     logfile << "v after = " << v << std::endl;
                //     logfile << "q = " << q << std::endl;
                // }

                for (arma::uword j = 0; j <= numRegressors; j++) {
                    double xij = x(j, i);
                    arr(j, numRegressors + 1) += xij * (y1[i] * (1.0 - v) + y0[i] * (-v));

                    for (arma::uword k = j; k <= numRegressors; k++)
                    {
                        arr(j, k) += xij * x(k, i) * q * (y0[i] + y1[i]);
                    }
                }
            }

	    for (arma::uword j = 1; j <= numRegressors; j++) {
		// double nextPar = (par[j] + arr(j, numRegressors + 1));
		double nextPar = par[j];
		ll -= lam / 2 * std::pow(nextPar, 2);
		arr(j, numRegressors + 1) -= lam * nextPar;
		arr(j, j) -= lam * ((nextPar > 0) - (nextPar < 0));
	    }

            if (llP == 1e+10)
            {
                llN = ll;
            }

            for (arma::uword j = 1; j <= numRegressors; j++)
            {
                for (arma::uword k = 0; k < j; k++)
                {
                    arr(j, k) = arr(k, j);
                }
            }

            for (arma::uword i = 0; i <= numRegressors; i++)
            {
                double s = arr(i, i);
                arr(i, i) = 1.0;
                for (arma::uword k = 0; k <= numRegressors + 1; k++)
                {
                    arr(i, k) = arr(i, k) / s;
                }

                for (arma::uword j = 0; j <= numRegressors; j++)
                {
                    if (i != j)
                    {
                        s = arr(j, i);
                        arr(j, i) = 0.0;
                        for (arma::uword k = 0; k <= numRegressors + 1; k++)
                        {
                            arr(j, k) = arr(j, k) - s * arr(i, k);
                        }
                    }
                }
            }

            for (arma::uword j = 0; j <= numRegressors; j++)
            {
                par[j] += arr(j, numRegressors + 1);
            }

	    // logfile << std::endl;
	    // logfile << "Params:" <<std::endl;
	    // for (int i = 0; i <= numRegressors; i++) {
	    // 	if (i==0) logfile << "intercept\t";
	    // 	else logfile << regressorNames[i-1] << "\t";
	    // }
	    // logfile << std::endl;
	    // for (int i = 0; i <= numRegressors; i++) {
	    // 	logfile << par[i] << "\t";
	    // }
	    // logfile << std::endl;
	    // logfile << "ll = " << ll << std::endl;
	    // logfile << std::endl;
	    
            iter++;
        }

        chiSq = llN - ll;

        //Indicates whether each coefficient is significant at the alpha level.
        // std::vector<std::string> sigMarker;
        // sigMarker.reserve(numRegressors);
        // arma::vec pValues = arma::vec(numRegressors + 1);
        // arma::vec zScores = arma::vec(numRegressors + 1);

        double zScore;

        for (arma::uword j = 1; j <= numRegressors; j++)
        {
            par[j] = par[j] / xStdDevs[j];
            parStdErr[j] = std::sqrt(arr(j, j)) / xStdDevs[j];
            par[0] = par[0] - par[j] * xMeans[j];
            zScore = par[j] / parStdErr[j];

            // logfile.open("../debug.log", std::ios_base::app);
            // if (numRegressors == 3 && regressorNames[0] == "X311MULTINOM.1" &&
            //  	 regressorNames[1] == "X311MULTINOM.2" &&regressorNames[2] == "X311MULTINOM.3") {
            if (std::isnan(zScore))
            {
                logfile << "Regressors:" << std::endl;
                for (int i = 0; i < numRegressors; i++) {
                    logfile << regressorNames[i] << "\t";
                }
                logfile << std::endl;
                for (int i = 0; i < numRegressors; i++) {
                    logfile << par[i] << "\t";
                }
                logfile << std::endl;
                logfile << "xStdDevs[" << j << "] = " << xStdDevs[j] << std::endl;
                logfile << "xMeans[" << j << "] = " << xMeans[j] << std::endl;
                logfile << "arr(" << j << "," << j << ") = " << arr(j,j) << std::endl;
                logfile << "par[" << j << "] = " << par[j] << std::endl;
                logfile << "parStdErr[" << j << "] = " << parStdErr[j] << std::endl;
                logfile << "zScore = " << zScore << std::endl << std::endl;
                break;
            }

            double prob = norm(std::abs(zScore));
            pValues[j] = prob;
            zScores[j] = zScore;
        }

        // if (std::isnan(zScore))
        // {
        //     lam = std::max(0.01, lam * 10);
        //     logfile << "Re-run, lam = " << lam <<  "\n\n";
        //     continue;
        // }

        parStdErr[0] = std::sqrt(arr(0, 0));
        zScore = par[0] / parStdErr[0];

        // if (std::isnan(zScore))
        // {
        //     lam = std::max(0.01, lam * 10);
        //     logfile << "Re-run, lam = " << lam <<  "\n\n";
        //     continue;
        // }

        pValues[0] = norm(zScore);
        zScores[0] = zScore;

        break;
    }

    // logfile.open("../test_results/debug.log", std::ios_base::app);
    // logfile << "zScore = " << zScore << std::endl;
    // logfile << "parStdErr[" << 0 << "] = " << parStdErr[0] << std::endl << std::endl;
    // logfile << "par[" << 0 << "] = " << par[0] << std::endl << std::endl;
    // logfile.close();

    // logfile.close();

    double intercept = par[0];
    coefficients = par;

    return LogisticRegressionResult(targetName, regressorNames, xMeans, xStdDevs,
				    numRegressors, ny0, ny1, coefficients, parStdErr,
				    pValues, intercept, ll, sigMarker, chiSq, alpha);
}

double LogisticRegression::norm(double z)
{
    // std::ofstream logfile;
    // logfile.open("../test_results/debug.log", std::ios_base::app);

    double q = z * z;
    const double pi = boost::math::constants::pi<double>();
    double piOver2 = pi / 2.0;

    // Rcpp::Rcout << "chisq = " << q << std::endl;

    if (std::abs(q) > 7.0)
    {
        return (1.0 - 1.0 / q + 3.0 / (q * q)) * std::exp(-q / 2.0) /
               (std::abs(z) * std::sqrt(piOver2));
    }
    else
    {
        boost::math::chi_squared dist(1);

        // logfile << "LOGISTIC REGRESSION" << std::endl;
        // logfile << "dist.df = " << dist.degrees_of_freedom() << std::endl;
        // logfile << "q = " << q << std::endl << std::endl;
        // logfile.close();
        double p = cdf(dist, q);
        return (p);
    }
}

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
// 		DiscreteVariable* x = (DiscreteVariable*)data.getVariable(i);
// 		Rcpp::Rcout << x->getName() << std::endl;
// 		std::vector<Variable*> regressors(data.getVariables());

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
