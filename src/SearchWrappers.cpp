#include "MGM.hpp"
#include "CoxMGM.hpp"
#include "PcStable.hpp"
#include "Fci.hpp"
#include "STEPS.hpp"
#include "STARS.hpp"
#include "StabilityUtils.hpp"
#include "Bootstrap.hpp"
#include "IndTestMultiCox.hpp"
#include "SearchCV.hpp"
#include "GrowShrink.hpp"
#include "DegenerateGaussianScore.hpp"
#include "RegressionBicScore.hpp"
#include "Grasp.hpp"
#include "Boss.hpp"
#include "Knowledge.hpp"

//' Calculate the MGM graph on a dataset
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param lambda A numeric vector of three values for the regularization parameter lambda: the first for continuous-continuous edges, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print updates on the progress of opti-\mizing MGM. The default is FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' sim <- rCausalMGM::simRandomDAG(200, 25)
//' g <- rCausalMGM::mgm(sim$data)
//` print(g)
// [[Rcpp::export]]
Rcpp::List mgm(
    const Rcpp::DataFrame &data, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }
    
    std::vector<double> l(lambda.begin(), lambda.end());

    int lamLength = 3;

    if (l.size() == 1) {
	for (int i = 1; i < lamLength; i++) {
	    l.push_back(l[0]);
	}
    } else if (l.size() != lamLength) {
	throw std::runtime_error("The regularization parameter lambda should be either a vector of length " + std::to_string(lamLength) + " or a single value for this dataset.");
    }

    MGM mgm(ds, l);
    mgm.setVerbose(verbose);
    EdgeListGraph mgmGraph = mgm.search();
    
    RcppThread::checkUserInterrupt();

    Rcpp::List result = mgmGraph.toList();

    return result;
}


//' Calculate the CoxMGM graph on a dataset. The dataset must contain at least one censored variable formatted as Surv object from the survival package.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. All censored variables must be a survival::Surv object. Any rows with missing values will be dropped.
//' @param lambda A numeric vector of five values for the regularization parameter lambda: the first for continuous-continuous edges, the second for continuous-discrete, the third for discrete-discrete, the fourth for continuous-survival, and the fifth for discrete-survival. Defaults to c(0.2, 0.2, 0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print updates on the progress of opti-\mizing MGM. The default is FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' sim <- rCausalMGM::simRandomDAG(200, 25)
//' time1 <- exp(0.5 * sim$data$X1 - 0.5 * sim$data$X2 + rnorm(nrow(sim$data)))
//' censtime1 <- sample(time1)
//' event1 <- as.integer(time1 < censtime1)
//' time1 <- pmin(time1, censtime1)
//' sim$data$Survival1 <- survival::Surv(time1, event1)
//' ig <- rCausalMGM::coxmgm(sim$data)
//` print(ig)
// [[Rcpp::export]]
Rcpp::List coxmgm(
    const Rcpp::DataFrame &data, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2, 0.2, 0.2),
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }
    
    bool v = verbose; 

    std::vector<double> l(lambda.begin(), lambda.end());

    int lamLength = 5;

    if (l.size() == 1) {
	for (int i = 1; i < lamLength; i++) {
	    l.push_back(l[0]);
	}
    } else if (l.size() != lamLength) {
	throw std::runtime_error("The regularization parameter lambda should be either a vector of length " + std::to_string(lamLength) + " or a single value for this dataset.");
    }

    CoxMGM mgm(ds, l);
    mgm.setVerbose(v);
    // mgm.calcLambdaMax();
    EdgeListGraph mgmGraph = mgm.search();
    
    RcppThread::checkUserInterrupt();

    // auto elapsedTime = mgm.getElapsedTime();

    // if (v) {
    // 	if (elapsedTime < 100*1000) {
    // 	    Rcpp::Rcout.precision(2);
    // 	} else {
    // 	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
    // 	}
    //     Rcpp::Rcout << "CoxMGM Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    // }

    Rcpp::List result = mgmGraph.toList();

    return result;
}


//' Calculate the solution path for an MGM graph on a dataset
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param lambdas A numeric vector containing the values of lambda to learn an MGM with. The default value is NULL, in which case a log-spaced vector of nLambda values for lambda will be supplied instead.
//' @param nLambda A numeric value indicating the number of lambda values to test when the lambdas vector is NULL.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return The a graphPath object that contains MGM graphs learned by the solution path, as well as the BIC and AIC selected models
//' @export
//' @examples
//' sim <- rCausalMGM::simRandomDAG(200, 25)
//' ig <- rCausalMGM::mgm(sim$data)
//` print(ig.path)
// [[Rcpp::export]]
Rcpp::List mgmPath(
    const Rcpp::DataFrame &data,
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 30,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    arma::vec _lambda;
    std::vector<double> l; 

    int n = ds.getNumRows();
    int p = ds.getNumColumns();

    MGM mgm(ds);
    mgm.setVerbose(verbose);

    double logLambdaMax = std::log10(mgm.calcLambdaMax());

    logLambdaMax = std::min(logLambdaMax,
    			    std::log10(10 * std::sqrt(std::log10(ds.getNumColumns()) /
						      ((double) ds.getNumRows()))));

    if (lambdas.isNotNull()) {
        _lambda = arma::vec(Rcpp::as<arma::vec>(lambdas)); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
	if (ds.getNumRows() > ds.getNumColumns()) {
	    _lambda = arma::logspace(logLambdaMax+std::log10(0.05), logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	} else {
	    _lambda = arma::logspace(logLambdaMax-1, logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	}
    }

    arma::vec loglik(l.size(), arma::fill::zeros);
    arma::vec nParams(l.size(), arma::fill::zeros);
    std::vector<EdgeListGraph> mgmGraphs = mgm.searchPath(l, loglik, nParams);

    RcppThread::checkUserInterrupt();

    Rcpp::List graphList;
    
    for (int i = 0; i < l.size(); i++) {
        graphList.push_back(mgmGraphs[i].toList());
    }

    arma::vec aic = 2*nParams - 2*loglik;
    arma::vec bic = std::log(n)*nParams - 2*loglik;

    arma::uword aicIdx = arma::index_min(aic);
    arma::uword bicIdx = arma::index_min(bic);
    
    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.bic"]=mgmGraphs[bicIdx].toList(),
					   Rcpp::_["graph.aic"]=mgmGraphs[aicIdx].toList(),
					   Rcpp::_["graphs"]=graphList,
					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
					   Rcpp::_["alphas"]=R_NilValue,
					   Rcpp::_["AIC"] = aic,
					   Rcpp::_["BIC"] = bic,
					   Rcpp::_["loglik"] = loglik,
					   Rcpp::_["nParams"] = nParams,
					   Rcpp::_["n"] = n);

    result.attr("class") = "graphPath";

    return result;
}


//' Calculate the solution path for a CoxMGM graph on a dataset
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. All censored variables must be a survival::Surv object. Any rows with missing values will be dropped.
//' @param lambdas A numeric vector containing the values of lambda to learn an MGM with. The default value is NULL, in which case a log-spaced vector of nLambda values for lambda will be supplied instead.
//' @param nLambda A numeric value indicating the number of lambda values to test when the lambdas vector is NULL.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return The a graphPath object that contains MGM graphs learned by the solution path, as well as the BIC and AIC selected models
//' @export
//' @examples
//' sim <- rCausalMGM::simRandomDAG(200, 25)
//' time1 <- exp(0.5 * sim$data$X1 - 0.5 * sim$data$X2 + rnorm(nrow(sim$data)))
//' censtime1 <- sample(time1)
//' event1 <- as.integer(time1 < censtime1)
//' time1 <- pmin(time1, censtime1)
//' sim$data$Survival1 <- survival::Surv(time1, event1)
//' ig.path <- rCausalMGM::coxmgmPath(sim$data)
//' print(ig.path)
// [[Rcpp::export]]
Rcpp::List coxmgmPath(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 30,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    bool v = verbose; // Rcpp::is_true(Rcpp::all(verbose));

    arma::vec _lambda;
    std::vector<double> l; // = {0.2, 0.2, 0.2};

    int n = ds.getNumRows();

    CoxMGM mgm(ds);
    mgm.setVerbose(v);

    double logLambdaMax = std::log10(mgm.calcLambdaMax());

    Rcpp::Rcout << "LambdaMax: " << std::pow(10, logLambdaMax) << std::endl;

    // return Rcpp::List::create();

    logLambdaMax = std::min(logLambdaMax,
    			    std::log10(10 * std::sqrt(std::log10(ds.getNumColumns()) /
    						      ((double) ds.getNumRows()))));

    if (lambdas.isNotNull()) {
        _lambda = arma::vec(Rcpp::as<arma::vec>(lambdas)); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
	if (ds.getNumRows() > ds.getNumColumns()) {
	    _lambda = arma::logspace(logLambdaMax+std::log10(0.05), logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	} else {
	    _lambda = arma::logspace(logLambdaMax-1, logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	}
    }

    // Rcpp::Rcout << "Beginning path search for lambdas " << _lambda.t() << std::endl;
    arma::vec loglik(l.size(), arma::fill::zeros);
    arma::vec nParams(l.size(), arma::fill::zeros);
    std::vector<EdgeListGraph> mgmGraphs = mgm.searchPath(l, loglik, nParams);

    RcppThread::checkUserInterrupt();

    // auto elapsedTime = mgm.getElapsedTime();

    // if (v) {
    // 	if (elapsedTime < 100*1000) {
    // 	    Rcpp::Rcout.precision(2);
    // 	} else {
    // 	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
    // 	}
    //     Rcpp::Rcout << "CoxMGM Path Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    // }

    Rcpp::List graphList;
    
    for (int i = 0; i < l.size(); i++) {
        graphList.push_back(mgmGraphs[i].toList());
    }

    arma::vec aic = 2*nParams - 2*loglik;
    arma::vec bic = std::log(n)*nParams - 2*loglik;

    arma::uword aicIdx = arma::index_min(aic);
    arma::uword bicIdx = arma::index_min(bic);
    
    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.bic"]=mgmGraphs[bicIdx].toList(),
					   Rcpp::_["graph.aic"]=mgmGraphs[aicIdx].toList(),
					   Rcpp::_["graphs"]=graphList,
					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
					   Rcpp::_["alphas"]=R_NilValue,
					   Rcpp::_["AIC"] = aic,
					   Rcpp::_["BIC"] = bic,
					   Rcpp::_["loglik"] = loglik,
					   Rcpp::_["nParams"] = nParams,
					   Rcpp::_["n"] = n);

    result.attr("class") = "graphPath";
    
    // Rcpp::List result = Rcpp::List::create(Rcpp::_["graphs"]=graphList,
    // 					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
    // 					   Rcpp::_["loglik"] = loglik,
    // 					   Rcpp::_["nParams"] = nParams,
    // 					   Rcpp::_["AIC"] = 2*nParams - 2*loglik,
    // 					   Rcpp::_["BIC"] = std::log(n)*nParams - 2*loglik);
    

    return result;
}


//' Calculate the solution path for an MGM graph on a dataset with k-fold cross-validation
//'
//' @param data The dataframe
//' @param lambdas A range of lambda values used to calculate a solution path for MGM. If NULL, lambdas is set to nLambda logarithmically spaced values from 10*sqrt(log10(p)/n) to sqrt(log10(p)/n). Defaults to NULL.
//' @param nLambda The number of lambda values to fit an MGM for when lambdas is NULL
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::mgmCV(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List mgmCV(
    const Rcpp::DataFrame &data,
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 30,
    const int nfolds = 5,
    Rcpp::Nullable<Rcpp::NumericVector> foldid = R_NilValue,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    bool v = verbose; // Rcpp::is_true(Rcpp::all(verbose));

    arma::vec _lambda;
    std::vector<double> l; // = {0.2, 0.2, 0.2};

    int n = ds.getNumRows();
    int p = ds.getNumColumns();

    MGM mgm(ds);
    mgm.setVerbose(v);

    double logLambdaMax = std::log10(mgm.calcLambdaMax());

    logLambdaMax = std::min(logLambdaMax,
    			    std::log10(10 * std::sqrt(std::log10(p) / ((double) n))));

    if (lambdas.isNotNull()) {
        _lambda = arma::vec(Rcpp::as<arma::vec>(lambdas)); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
	if (n > p) {
	    _lambda = arma::logspace(logLambdaMax+std::log10(0.05), logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	} else {
	    _lambda = arma::logspace(logLambdaMax-1, logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	}
    }
    
    arma::uvec _foldid;
    
    if (foldid.isNull()) {
	// Rcpp::NumericVector folds = Rcpp::as<Rcpp::NumericVector>(arma::linspace(1, nfolds));
	// folds = Rcpp::rep_len(folds, n);
	// arma::vec folds = arma::linspace(1, nfolds, nfolds);
	// folds = arma::repmat(folds, 1, std::ceil(n/((double)nfolds)));
	// _foldid = arma::conv_to<arma::uvec>::from(folds(Rcpp::as<arma::uvec>(Rcpp::sample(n,n))-1));
	_foldid = arma::linspace<arma::uvec>(1, n, n) - 1;
	_foldid.transform([nfolds](arma::uword val) { return val % nfolds + 1; });
	_foldid = arma::shuffle(_foldid);
    } else {
	_foldid = Rcpp::as<arma::uvec>(foldid);
    }

    if (_foldid.n_elem != n) {
	throw std::invalid_argument("foldid has a length that does not equal number of samples in the dataset.");
    }

    if (!SearchCV::checkFoldID(_foldid)) {
	throw std::invalid_argument("Invalid input for foldid. Values must be in the range 1:nfolds, with folds of approximately equal size.");
    }

    // Rcpp::Rcout << "Fold ids: " << _foldid.t() << std::endl;
    // Rcpp::Rcout << "Beginning path search for lambdas " << _lambda.t() << std::endl;
    arma::mat loglik(l.size(), arma::max(_foldid), arma::fill::zeros);
    arma::uvec index(2, arma::fill::zeros);
    // arma::vec nParams(l.size(), arma::fill::zeros);
    std::vector<EdgeListGraph> cvGraphs = mgm.searchPathCV(l, _foldid, loglik, index);

    RcppThread::checkUserInterrupt();

    // double elapsedTime = mgm.getElapsedTime();

    // if (verbose) {
    // 	double factor = (elapsedTime < 10) ? std::pow(10, 2 - std::ceil(std::log10(std::abs(elapsedTime)))) : 1.0;
    // 	elapsedTime = std::round(elapsedTime * factor) / factor;
    //     Rcpp::Rcout << "  MGM Cross-Validation Elapsed Time =  " << elapsedTime << " s" << std::endl;
    // }


    // auto elapsedTime = mgm.getElapsedTime();

    // if (v) {
    // 	if (elapsedTime < 100*1000) {
    // 	    Rcpp::Rcout.precision(2);
    // 	} else {
    // 	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
    // 	}
    //     Rcpp::Rcout << "MGM Cross-Validation Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    // }

    // Rcpp::List graphList;
    
    // for (int i = 0; i < l.size(); i++) {
    //     graphList.push_back(mgmGraphs[i].toList());

    // }

    // uint minIdx = arma::find(lambdas == cvGraphs[0].getHyperParam("lambda")[1])[0];
    // uint seIdx = arma::find(lambdas == cvGraphs[1].getHyperParam("lambda")[1])[0];

    // Rcpp::List result;
    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.min"]=cvGraphs[0].toList(),
					   Rcpp::_["graph.1se"]=cvGraphs[1].toList(),
    					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
					   Rcpp::_["alphas"]=R_NilValue,
					   Rcpp::_["orientRules"]=R_NilValue,
					   Rcpp::_["idx.min"]=index(0) + 1,
					   Rcpp::_["idx.1se"]=index(1) + 1,
					   Rcpp::_["mean"] = arma::mean(loglik, 1),
					   Rcpp::_["se"] = arma::stddev(loglik, 0, 1),
					   Rcpp::_["size"] = R_NilValue,
					   Rcpp::_["foldid"]=_foldid);

    result.attr("class") = "graphCV";

    return result;
}



//' Calculates the optimal lambda values for the MGM algorithm using StEPS and run the algorithm using those values. Optimal values are printed
//'
//' @param data The dataframe
//' @param lambdas A range of lambda values assessed for stability by the StEPS algorithm. If NULL, lambdas is set to nLambda logarithmically spaced values from 10*sqrt(log10(p)/n) to sqrt(log10(p)/n). Defaults to NULL.
//' @param nLambda The number of lambda values to fit an MGM for when lambdas is NULL
//' @param g The gamma parameter for STEPS. Defaults to 0.05
//' @param numSub The number of subsets to split the data into. Defaults to 20
//' @param subSize The size of the subsamples used for STEPS. If the value is -1, the size of the subsamples is set to floor(10*sqrt(n)). If the value is in the range (0,1), the size of the subsamples is set to floor(subSize * n). Otherwise, if subSize is in the range [1,n), the size of the subsamples is set to subSize. Defaults to -1.
//' @param leaveOneOut If TRUE, performs leave-one-out subsampling. Defaults to FALSE.
//' @param computeStabs If TRUE, stability values are calculated. Defaults to FALSE.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::steps(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List steps(
    const Rcpp::DataFrame &data, 
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 30,
    const double gamma = 0.05,
    const int numSub = 20,
    const int subSize = -1,
    const bool leaveOneOut = false,
    const bool computeStabs = false,
    const int threads = -1,
    const bool rank = false,
    const bool verbose = false
) {

    std::vector<double> l;
    arma::vec _lambda;
    
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    MGM mgm(ds);

    double logLambdaMax = std::log10(mgm.calcLambdaMax());

    logLambdaMax = std::min(logLambdaMax,
    			    std::log10(10 * std::sqrt(std::log10(ds.getNumColumns()) /
						      ((double) ds.getNumRows()))));
    
    if (lambdas.isNotNull()) {
        _lambda = arma::vec(Rcpp::as<arma::vec>(lambdas)); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
	if (ds.getNumRows() > ds.getNumColumns()) {
	    _lambda = arma::logspace(logLambdaMax+std::log10(0.05), logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	} else {
	    _lambda = arma::logspace(logLambdaMax-1, logLambdaMax, nLambda); 
	    l = std::vector<double>(_lambda.begin(), _lambda.end());
	}
    }

    STEPS steps;
    if (subSize < 0)
	steps = STEPS(ds, l, gamma, numSub, leaveOneOut);
    else
	steps = STEPS(ds, l, gamma, numSub, subSize, leaveOneOut);
      
    if (threads > 0) steps.setThreads(threads);
    steps.setComputeStabs(computeStabs);
    steps.setVerbose(verbose);

    arma::mat instabs;
    if (ds.isCensored()) {
	instabs = arma::mat(l.size(), 6);
    } else {
	instabs = arma::mat(l.size(), 4);
    }
    instabs.fill(arma::datum::nan);

    arma::umat samps;

    std::vector<EdgeListGraph> graphs = steps.runStepsPath(instabs, samps);
    Rcpp::List graphSteps = graphs.at(0).toList();

    if (computeStabs) {
        graphSteps["stabilities"] = steps.getStabs();
        std::vector<std::string> names = ds.getVariableNames();
        Rcpp::rownames(graphSteps["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
        Rcpp::colnames(graphSteps["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.steps"]=graphSteps,
					   Rcpp::_["graph.stars"]=graphs.at(1).toList(),
    					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
					   Rcpp::_["gamma"]=gamma,
    					   Rcpp::_["instability"] = instabs,
					   Rcpp::_["subsamples"] = samps);

    result.attr("class") = "graphSTEPS";

    return result;

}

//' Runs the causal algorithm PC-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List pcStars(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority"),
    Rcpp::Nullable<Rcpp::NumericVector> alphas = R_NilValue,
    const double gamma = 0.01,
    const int numSub = 20,
    const int subSize = -1,
    const bool leaveOneOut = false,
    const bool computeStabs = false,
    const int threads = -1,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    arma::vec _alphas;

    if (alphas.isNotNull()) {
        _alphas = arma::vec(Rcpp::as<arma::vec>(alphas)); 
    } else {
      _alphas = {0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2 };
    }

    std::vector<OrientRule> rules;

    for (std::string rule : _orientRule) {
	rules.push_back(SepsetProducer::str2rule(rule));
    }

    STARS stars;

    if (subSize < 0)
	stars = STARS(ds, "pc", _alphas, gamma, numSub, leaveOneOut);
    else
	stars = STARS(ds, "pc", _alphas, gamma, numSub, subSize, leaveOneOut);

    stars.setVerbose(verbose);
    stars.setOrientRule(rules.at(0));
    // stars.setFDR(fdr);

    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        stars.setInitialGraph(&ig);
    }
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	stars.setKnowledge(k);
    }

    arma::mat instabs = arma::mat(_alphas.n_elem, 1);
    instabs.fill(arma::datum::nan);

    arma::umat samps;

    Rcpp::List graphStars = stars.runStarsPar(instabs, samps).toList();

    if (computeStabs) {
        graphStars["stabilities"] = stars.getStabs();
        std::vector<std::string> names = ds.getVariableNames();
        Rcpp::rownames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
        Rcpp::colnames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph"]=graphStars,
    					   Rcpp::_["alphas"]=arma::sort(_alphas, "ascend"),
					   Rcpp::_["gamma"]=gamma,
    					   Rcpp::_["instability"] = instabs,
					   Rcpp::_["subsamples"] = samps);

    result.attr("class") = "graphSTARS";

    return result;
}


//' Runs the causal algorithm FCI-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::fciStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List fciStars(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority"),
    Rcpp::Nullable<Rcpp::NumericVector> alphas = R_NilValue,
    const double gamma = 0.01,
    const int numSub = 20,
    const int subSize = -1,
    const bool leaveOneOut = false,
    const bool computeStabs = false,
    const int threads = -1,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    arma::vec _alphas;

    if (alphas.isNotNull()) {
        _alphas = arma::vec(Rcpp::as<arma::vec>(alphas)); 
    } else {
      _alphas = {0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2 };
    }

    std::vector<OrientRule> rules;

    for (std::string rule : _orientRule) {
	rules.push_back(SepsetProducer::str2rule(rule));
    }

    STARS stars;

    if (subSize < 0)
	stars = STARS(ds, "fci", _alphas, gamma, numSub, leaveOneOut);
    else
	stars = STARS(ds, "fci", _alphas, gamma, numSub, subSize, leaveOneOut);

    stars.setVerbose(verbose);
    stars.setOrientRule(rules.at(0));
    // stars.setFDR(fdr);

    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        stars.setInitialGraph(&ig);
    }
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	stars.setKnowledge(k);
    }

    arma::mat instabs = arma::mat(_alphas.n_elem, 1);
    instabs.fill(arma::datum::nan);

    arma::umat samps;

    Rcpp::List graphStars = stars.runStarsPar(instabs, samps).toList();

    if (computeStabs) {
        graphStars["stabilities"] = stars.getStabs();
        std::vector<std::string> names = ds.getVariableNames();
        Rcpp::rownames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
        Rcpp::colnames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph"]=graphStars,
    					   Rcpp::_["alphas"]=arma::sort(_alphas, "ascend"),
					   Rcpp::_["gamma"]=gamma,
    					   Rcpp::_["instability"] = instabs,
					   Rcpp::_["subsamples"] = samps);

    result.attr("class") = "graphSTARS";

    return result;
}


// //' Runs the causal algorithm PC-Stable on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List pcStable(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
//     const double alpha = 0.05,
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }

//     // bool v = Rcpp::is_true(Rcpp::all(verbose));
//     // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));

//     IndTestMultiCox itm(ds, alpha);

//     PcStable pcs((IndependenceTest*) &itm);
//     if (threads > 0) pcs.setThreads(threads);
//     pcs.setVerbose(verbose);
//     pcs.setFDR(fdr);
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         pcs.setInitialGraph(&ig);
//     }
//     Knowledge k;
//     if (!knowledge.isNull()) {
// 	Rcpp::List _knowledge(knowledge);
// 	k = Knowledge(ds.getVariables(), _knowledge);
// 	pcs.setKnowledge(k);
//     }

//     Rcpp::List result = pcs.search().toList();

//     // ds.deleteVariables();

//     return result;
// }

//' Runs the causal algorithm PC-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List pcStable(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority"),
    const double alpha = 0.05,
    const int threads = -1,
    const bool fdr = false,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    
    IndTestMultiCox itm(ds, alpha);
    
    PcStable pc((IndependenceTest*) &itm);
    if (threads > 0) pc.setThreads(threads);
    pc.setVerbose(verbose);
    pc.setFDR(fdr);
    pc.setOrientRule(SepsetProducer::str2rule(_orientRule.at(0)));
    
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        pc.setInitialGraph(&ig);
    }
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	pc.setKnowledge(k);
    }

    EdgeListGraph g = pc.search();
    Rcpp::List result;
    if (_orientRule.size() == 1) {
	result = g.toList();
    } else {
	bool first = true;
	for (std::string rule : _orientRule) {
	    if (first) {
		first = false;
	    } else {
		g = pc.reorientWithRule(SepsetProducer::str2rule(rule));
	    }
	    result.push_back(g.toList(), rule);
	}
    }
    
    return result;
}


//' Runs the causal algorithm PC-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List pcCV(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp", "conservative"),
    Rcpp::Nullable<Rcpp::NumericVector> alphas = R_NilValue,
    const int nfolds = 5,
    Rcpp::Nullable<Rcpp::NumericVector> foldid = R_NilValue,
    const int threads = -1,
    const bool fdr = false,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    arma::vec _alphas;

    if (alphas.isNotNull()) {
        _alphas = arma::vec(Rcpp::as<arma::vec>(alphas)); 
    } else {
      _alphas = {0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2 };
    }

    std::vector<OrientRule> rules;

    for (std::string rule : _orientRule) {
	rules.push_back(SepsetProducer::str2rule(rule));
    }

    SearchCV cv;

    arma::uvec _foldid;
    
    if (foldid.isNull()) {
	cv = SearchCV(ds, "pc", nfolds);
    } else {
	_foldid = Rcpp::as<arma::uvec>(foldid);
	cv = SearchCV(ds, "pc", _foldid);
    }

    cv.setVerbose(verbose);
    cv.setAlphas(_alphas);
    cv.setOrientRules(rules);
    cv.setFDR(fdr);

    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        cv.setInitialGraph(&ig);
    }
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	cv.setKnowledge(k);
    }

    std::vector<EdgeListGraph> cvGraphs = cv.causalCV();

    std::set<CvResult> results = cv.getCVResults();

    int nValues = results.size();
    arma::vec mbSize(nValues);
    arma::vec mean(nValues);
    arma::vec se(nValues);
    arma::vec alphasOut(nValues);
    std::vector<std::string> rulesOut(nValues);

    uint idx = 0;
    uint minIdx = 0;
    CvResult minResult(1e20, 1e20, 0);
    for (auto res : results) {
	mbSize(idx) = res.mbSize;
	mean(idx) = res.mean;
	se(idx) = res.se;
	alphasOut(idx) = res.alpha;
	rulesOut.at(idx) = SepsetProducer::rule2str(res.rule);
	if (res.mean < minResult.mean) {
	    minIdx = idx;
	    minResult = res;
	}
	idx++;
    }

    uint seIdx = 0;
    for (auto res : results) {
	if (res.mean < minResult.mean + minResult.se) {
	    // seResult = res;
	    break;
	}
	seIdx++;
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.min"]=cvGraphs[0].toList(),
					   Rcpp::_["graph.1se"]=cvGraphs[1].toList(),
    					   Rcpp::_["lambdas"]=R_NilValue,
					   Rcpp::_["alphas"]=alphasOut,
					   Rcpp::_["orientRules"]=rulesOut,
					   Rcpp::_["idx.min"]=minIdx+1,
					   Rcpp::_["idx.1se"]=seIdx+1,
					   Rcpp::_["mean"] = mean,
					   Rcpp::_["se"] = se,
					   Rcpp::_["size"] = mbSize,
					   Rcpp::_["foldid"]=cv.getFoldID());

    result.attr("class") = "graphCV";
        
    return result;
}


//' Runs the causal algorithm FCI-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::fciStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List fciCV(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp", "conservative"),
    Rcpp::Nullable<Rcpp::NumericVector> alphas = R_NilValue,
    const int nfolds = 5,
    Rcpp::Nullable<Rcpp::NumericVector> foldid = R_NilValue,
    const int threads = -1,
    const bool fdr = false,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    arma::vec _alphas;

    if (alphas.isNotNull()) {
        _alphas = arma::vec(Rcpp::as<arma::vec>(alphas)); 
    } else {
      _alphas = {0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2 };
    }
    
    std::vector<OrientRule> rules;

    for (std::string rule : _orientRule) {
	rules.push_back(SepsetProducer::str2rule(rule));
    }

    SearchCV cv;

    arma::uvec _foldid;
    
    if (foldid.isNull()) {
	cv = SearchCV(ds, "fci", nfolds);
    } else {
	_foldid = Rcpp::as<arma::uvec>(foldid);
	cv = SearchCV(ds, "fci", _foldid);
    }

    cv.setVerbose(verbose);
    cv.setAlphas(_alphas);
    cv.setOrientRules(rules);
    cv.setFDR(fdr);

    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        cv.setInitialGraph(&ig);
    }
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	cv.setKnowledge(k);
    }

    std::vector<EdgeListGraph> cvGraphs = cv.causalCV();

    std::set<CvResult> results = cv.getCVResults();

    int nValues = results.size();
    arma::vec mbSize(nValues);
    arma::vec mean(nValues);
    arma::vec se(nValues);
    arma::vec alphasOut(nValues);
    std::vector<std::string> rulesOut(nValues);

    uint idx = 0;
    uint minIdx = 0;
    CvResult minResult(1e20, 1e20, 0);
    for (auto res : results) {
	mbSize(idx) = res.mbSize;
	mean(idx) = res.mean;
	se(idx) = res.se;
	alphasOut(idx) = res.alpha;
	rulesOut.at(idx) = SepsetProducer::rule2str(res.rule);
	if (res.mean < minResult.mean) {
	    minIdx = idx;
	    minResult = res;
	}
	idx++;
    }

    uint seIdx = 0;
    for (auto res : results) {
	if (res.mean < minResult.mean + minResult.se) {
	    // seResult = res;
	    break;
	}
	seIdx++;
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.min"]=cvGraphs[0].toList(),
					   Rcpp::_["graph.1se"]=cvGraphs[1].toList(),
    					   Rcpp::_["lambdas"]=R_NilValue,
					   Rcpp::_["alphas"]=alphasOut,
					   Rcpp::_["orientRules"]=rulesOut,
					   Rcpp::_["idx.min"]=minIdx+1,
					   Rcpp::_["idx.1se"]=seIdx+1,
					   Rcpp::_["mean"] = mean,
					   Rcpp::_["se"] = se,
					   Rcpp::_["size"] = mbSize,
					   Rcpp::_["foldid"]=cv.getFoldID());

    result.attr("class") = "graphCV";
        
    return result;
}


//' Runs the causal algorithm PC-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List mgmpcCV(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const std::string cvType = "random",
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp", "conservative"),
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 20,
    Rcpp::Nullable<Rcpp::NumericVector> alphas = R_NilValue,
    const int numPoints = 60,
    const int nfolds = 5,
    Rcpp::Nullable<Rcpp::NumericVector> foldid = R_NilValue,
    const int threads = -1,
    const bool fdr = false,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    arma::vec _alphas;

    if (alphas.isNotNull()) {
        _alphas = arma::vec(Rcpp::as<arma::vec>(alphas)); 
    } else {
      _alphas = {0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2 };
    }

    std::vector<OrientRule> rules;

    for (std::string rule : _orientRule) {
	rules.push_back(SepsetProducer::str2rule(rule));
    }

    // std::vector<double> l;
    arma::vec _lambdas;

    MGM mgm(ds);

    double lambdaMax = mgm.calcLambdaMax();

    lambdaMax = std::min(lambdaMax,
			 10 * std::sqrt(std::log10(ds.getNumColumns()) /
					((double) ds.getNumRows())));
    
    if (lambdas.isNotNull()) {
        _lambdas = arma::vec(Rcpp::as<arma::vec>(lambdas)); 
	// l = std::vector<double>(_lambdas.begin(), _lambdas.end());
    } else {
	if (ds.getNumRows() > ds.getNumColumns()) {
	    _lambdas = arma::linspace(lambdaMax*0.05, lambdaMax, nLambda); 
	    // l = std::vector<double>(_lambdas.begin(), _lambdas.end());
	} else {
	    _lambdas = arma::linspace(lambdaMax*0.1, lambdaMax, nLambda); 
	    // l = std::vector<double>(_lambdas.begin(), _lambdas.end());
	}
    }


    SearchCV cv;

    arma::uvec _foldid;
    
    if (foldid.isNull()) {
	cv = SearchCV(ds, "pc", nfolds);
    } else {
	_foldid = Rcpp::as<arma::uvec>(foldid);
	cv = SearchCV(ds, "pc", _foldid);
    }

    cv.setVerbose(verbose);
    cv.setOrientRules(rules);
    cv.setFDR(fdr);
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	cv.setKnowledge(k);
    }

    std::vector<EdgeListGraph> cvGraphs;

    if (cvType=="grid") {
	
	cv.setAlphas(_alphas);
	cv.setLambdas(_lambdas);
	cvGraphs = cv.causalMGMGridCV();
	
    } else if (cvType=="random") {
	
	double aMax = _alphas.max();
	double aMin = _alphas.min();
	double lMax = _lambdas.max();
	double lMin = _lambdas.min();

	_lambdas = (lMax - lMin) * arma::randu(numPoints) + lMin;
	_alphas = (aMax - aMin) * arma::randu(numPoints) + aMin;
	
	cv.setAlphas(_alphas);
	cv.setLambdas(_lambdas);
	cvGraphs = cv.causalMGMRandCV();
	
    } else {
	throw std::invalid_argument("CV Type must be one of \"grid\" or \"random\"");
    }

    std::set<CvResult> results = cv.getCVResults();

    int nValues = results.size();
    arma::vec mbSize(nValues);
    arma::vec mean(nValues);
    arma::vec se(nValues);
    arma::vec alphasOut(nValues);
    arma::vec lambdasOut(nValues);
    std::vector<std::string> rulesOut(nValues);

    uint idx = 0;
    uint minIdx = 0;
    CvResult minResult(1e20, 1e20, 0);
    for (auto res : results) {
	mbSize(idx) = res.mbSize;
	mean(idx) = res.mean;
	se(idx) = res.se;
	alphasOut(idx) = res.alpha;
	lambdasOut(idx) = res.lambda.at(0);
	rulesOut.at(idx) = SepsetProducer::rule2str(res.rule);
	if (res.mean < minResult.mean) {
	    minIdx = idx;
	    minResult = res;
	}
	idx++;
    }

    uint seIdx = 0;
    for (auto res : results) {
	if (res.mean < minResult.mean + minResult.se) {
	    // seResult = res;
	    break;
	}
	seIdx++;
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.min"]=cvGraphs[0].toList(),
					   Rcpp::_["graph.1se"]=cvGraphs[1].toList(),
    					   Rcpp::_["lambdas"]=lambdasOut,
					   Rcpp::_["alphas"]=alphasOut,
					   Rcpp::_["orientRules"]=rulesOut,
					   Rcpp::_["idx.min"]=minIdx+1,
					   Rcpp::_["idx.1se"]=seIdx+1,
					   Rcpp::_["mean"] = mean,
					   Rcpp::_["se"] = se,
					   Rcpp::_["size"] = mbSize,
					   Rcpp::_["foldid"]=cv.getFoldID());

    result.attr("class") = "graphCV";
        
    return result;
}


//' Runs the causal algorithm FCI-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::fciStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List mgmfciCV(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const std::string cvType = "grid",
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp", "conservative"),
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 20,
    Rcpp::Nullable<Rcpp::NumericVector> alphas = R_NilValue,
    const int numPoints = 60,
    const int nfolds = 5,
    Rcpp::Nullable<Rcpp::NumericVector> foldid = R_NilValue,
    const int threads = -1,
    const bool fdr = false,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    arma::vec _alphas;

    if (alphas.isNotNull()) {
        _alphas = arma::vec(Rcpp::as<arma::vec>(alphas)); 
    } else {
      _alphas = {0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2 };
    }
    
    std::vector<OrientRule> rules;

    for (std::string rule : _orientRule) {
	rules.push_back(SepsetProducer::str2rule(rule));
    }

    // std::vector<double> l;
    arma::vec _lambdas;

    MGM mgm(ds);

    double lambdaMax = mgm.calcLambdaMax();

    lambdaMax = std::min(lambdaMax,
			 10 * std::sqrt(std::log10(ds.getNumColumns()) /
					((double) ds.getNumRows())));
    
    if (lambdas.isNotNull()) {
        _lambdas = arma::vec(Rcpp::as<arma::vec>(lambdas)); 
	// l = std::vector<double>(_lambdas.begin(), _lambdas.end());
    } else {
	if (ds.getNumRows() > ds.getNumColumns()) {
	    _lambdas = arma::linspace(lambdaMax*0.05, lambdaMax, nLambda); 
	    // l = std::vector<double>(_lambdas.begin(), _lambdas.end());
	} else {
	    _lambdas = arma::linspace(lambdaMax*0.1, lambdaMax, nLambda); 
	    // l = std::vector<double>(_lambdas.begin(), _lambdas.end());
	}
    }


    SearchCV cv;

    arma::uvec _foldid;
    
    if (foldid.isNull()) {
	cv = SearchCV(ds, "fci", nfolds);
    } else {
	_foldid = Rcpp::as<arma::uvec>(foldid);
	cv = SearchCV(ds, "fci", _foldid);
    }

    cv.setVerbose(verbose);
    cv.setOrientRules(rules);
    cv.setFDR(fdr);
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	cv.setKnowledge(k);
    }

    std::vector<EdgeListGraph> cvGraphs;

    if (cvType=="grid") {
	
	cv.setAlphas(_alphas);
	cv.setLambdas(_lambdas);
	cvGraphs = cv.causalMGMGridCV();
	
    } else if (cvType=="random") {
	
	double aMax = _alphas.max();
	double aMin = _alphas.min();
	double lMax = _lambdas.max();
	double lMin = _lambdas.min();

	_lambdas = (lMax - lMin) * arma::randu(numPoints) + lMin;
	_alphas = (aMax - aMin) * arma::randu(numPoints) + aMin;
	
	cv.setAlphas(_alphas);
	cv.setLambdas(_lambdas);
	cvGraphs = cv.causalMGMRandCV();
	
    } else {
	throw std::invalid_argument("CV Type must be one of \"grid\" or \"random\"");
    }

    std::set<CvResult> results = cv.getCVResults();

    int nValues = results.size();
    arma::vec mbSize(nValues);
    arma::vec mean(nValues);
    arma::vec se(nValues);
    arma::vec alphasOut(nValues);
    arma::vec lambdasOut(nValues);
    std::vector<std::string> rulesOut(nValues);

    uint idx = 0;
    uint minIdx = 0;
    CvResult minResult(1e20, 1e20, 0);
    for (auto res : results) {
	mbSize(idx) = res.mbSize;
	mean(idx) = res.mean;
	se(idx) = res.se;
	alphasOut(idx) = res.alpha;
	lambdasOut(idx) = res.lambda.at(0);
	rulesOut.at(idx) = SepsetProducer::rule2str(res.rule);
	if (res.mean < minResult.mean) {
	    minIdx = idx;
	    minResult = res;
	}
	idx++;
    }

    uint seIdx = 0;
    for (auto res : results) {
	if (res.mean < minResult.mean + minResult.se) {
	    // seResult = res;
	    break;
	}
	seIdx++;
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.min"]=cvGraphs[0].toList(),
					   Rcpp::_["graph.1se"]=cvGraphs[1].toList(),
    					   Rcpp::_["lambdas"]=lambdasOut,
					   Rcpp::_["alphas"]=alphasOut,
					   Rcpp::_["orientRules"]=rulesOut,
					   Rcpp::_["idx.min"]=minIdx+1,
					   Rcpp::_["idx.1se"]=seIdx+1,
					   Rcpp::_["mean"] = mean,
					   Rcpp::_["se"] = se,
					   Rcpp::_["size"] = mbSize,
					   Rcpp::_["foldid"]=cv.getFoldID());

    result.attr("class") = "graphCV";
        
    return result;
}



// //' Runs the causal algorithm PC-Stable on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::fciStable(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List pcStableBSC(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
//     const int N = 250,
//     const int threads = -1,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }

//     std::string rule;
//     // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
//     // bool v = Rcpp::is_true(Rcpp::all(verbose));

//     // std::vector<std::string> validNames = {"majority", "maxp", "conservative", "sepsets"};

//     // rule = "sepsets"; // orientRule[0];

//     // std::transform(rule.begin(), rule.end(), rule.begin(),
//     // 		   [](unsigned char c){ return std::tolower(c); });

//     // if (std::find(validNames.begin(), validNames.end(), rule) == validNames.end())
//     // 	throw std::invalid_argument("Orientation rule must be one of {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    
//     BayesIndTestMultiCox itm(ds, 0.05);
    
//     PcStable pc((IndependenceTest*) &itm);
//     if (threads > 0) pc.setThreads(threads);
//     pc.setVerbose(verbose);
//     pc.setOrientRule(ORIENT_SEPSETS);
//     // pc.setFDR(fdr);
//     // if (rule == "majority")
//     // 	pc.setOrientRule(ORIENT_MAJORITY);
//     // if (rule == "maxp")
//     // 	pc.setOrientRule(ORIENT_MAXP);
//     // if (rule == "conservative")
//     // 	pc.setOrientRule(ORIENT_CONSERVATIVE);
//     // if (rule == "sepsets")
//     // 	pc.setOrientRule(ORIENT_SEPSETS);
    
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         pc.setInitialGraph(&ig);
//     }
//     Knowledge k;
//     if (!knowledge.isNull()) {
// 	Rcpp::List _knowledge(knowledge);
// 	k = Knowledge(ds.getVariables(), _knowledge);
// 	pc.setKnowledge(k);
//     }

//     std::vector<EdgeListGraph> sampledGraphs;
//     arma::vec scores(N);

//     for (int i = 0; i < N; i++) {
// 	sampledGraphs.push_back(pc.search());
// 	scores(i) = pc.getScore();
//     }

//     Rcpp::List graphList;
    
//     for (int i = 0; i < N; i++) {
//         graphList.push_back(sampledGraphs[i].toList());
//     }

//     arma::vec probs = arma::exp(scores - arma::max(scores));

//     probs /= arma::accu(probs);

//     arma::uword mapIdx = arma::index_max(scores);
    
//     Rcpp::List result = Rcpp::List::create(
// 	Rcpp::_["graph.map"]=sampledGraphs[mapIdx].toList(),
// 	Rcpp::_["graphs"]=graphList,
// 	Rcpp::_["scores"]=scores,
// 	Rcpp::_["graphPosterior"]=probs);
    
//     // ds.deleteVariables();
    
//     return result;
// }




// //' Calculate the solution path for an PC graph on a dataset
// //'
// //' @param data The dataframe
// //' @param lambdas A range of lambda values used to calculate a solution path for MGM. If NULL, lambdas is set to nLambda logarithmically spaced values from 10*sqrt(log10(p)/n) to sqrt(log10(p)/n). Defaults to NULL.
// //' @param nLambda The number of lambda values to fit an MGM for when lambdas is NULL
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated MGM graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' g <- rCausalMGM::pcPath(data.n100.p25)
// // [[Rcpp::export]]
// Rcpp::List pcPath(
//     const Rcpp::DataFrame& data,
//     Rcpp::NumericVector alphas = Rcpp::NumericVector::create(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1),
//     Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("max", "conservative", "majority", "none"),
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }

//     int n = ds.getNumRows();
//     int p = ds.getNumColumns();

//     std::vector<double> a(alphas.begin(), alphas.end());

//     arma::vec _alphas(a);

//     _alphas = arma::sort(_alphas, "ascend");

//     arma::vec loglik(_alphas.size(), arma::fill::zeros);
//     arma::vec nParams(_alphas.size(), arma::fill::zeros);
//     std::vector<EdgeListGraph> pcGraphs;
//     std::vector<CausalMGMParams> pcParams;
//     std::vector<double> l = { 0.5, 0.5, 0.5 };

//     for (arma::uword i = 0; i < _alphas.n_elem; i++) {

// 	IndTestMultiCox itm(ds, _alphas(i));

// 	if (orientRule[0]=="none") {

// 	    PcStable pcs((IndependenceTest*) &itm);
// 	    if (threads > 0) pcs.setThreads(threads);
// 	    pcs.setVerbose(verbose);
// 	    pcs.setFDR(fdr);
// 	    pcGraphs.push_back(pcs.search());

// 	} else if (orientRule[0]=="max") {

// 	    PcMax pcm((IndependenceTest*) &itm);
// 	    if (threads > 0) pcm.setThreads(threads);
// 	    pcm.setVerbose(verbose);
// 	    pcm.setFDR(fdr);
// 	    pcGraphs.push_back(pcm.search());

// 	} else if (orientRule[0]=="conservative") {

// 	    CpcStable cpc((IndependenceTest*) &itm);
// 	    if (threads > 0) cpc.setThreads(threads);
// 	    cpc.setVerbose(verbose);
// 	    cpc.setFDR(fdr);
// 	    pcGraphs.push_back(cpc.search());

// 	} else if (orientRule[0]=="majority") {

// 	    Pc50 pc50((IndependenceTest*) &itm);
// 	    if (threads > 0) pc50.setThreads(threads);
// 	    pc50.setVerbose(verbose);
// 	    pc50.setFDR(fdr);
// 	    pcGraphs.push_back(pc50.search());

// 	}

// 	CausalMGM causalMGM(ds, pcGraphs.at(i));
// 	causalMGM.setVerbose(verbose);
// 	causalMGM.setLambda(l);
// 	pcParams.push_back(causalMGM.search());

// 	arma::vec par(pcParams.at(i).toMatrix1D());
// 	loglik(i) = -n * causalMGM.smoothValue(par);
// 	nParams(i) = causalMGM.getNParams();
	
// 	RcppThread::checkUserInterrupt();
//     }

//     // auto elapsedTime = mgm.getElapsedTime();

//     // if (v) {
//     // 	if (elapsedTime < 100*1000) {
//     // 	    Rcpp::Rcout.precision(2);
//     // 	} else {
//     // 	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
//     // 	}
//     //     Rcpp::Rcout << "MGM Path Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
//     // }

//     Rcpp::List graphList;
    
//     for (int i = 0; i < a.size(); i++) {
// 	Rcpp::List pcGraph = pcGraphs[i].toList();
// 	pcGraph["parameters"] = pcParams[i].toList();
//         graphList.push_back(pcGraph);
//     }

//     arma::vec aic = 2*nParams - 2*loglik;
//     arma::vec bic = std::log(n)*nParams - 2*loglik;

//     arma::uword aicIdx = arma::index_min(aic);
//     arma::uword bicIdx = arma::index_min(bic);
    
//     Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.bic"]=graphList[bicIdx],
// 					   Rcpp::_["graph.aic"]=graphList[aicIdx],
// 					   Rcpp::_["graphs"]=graphList,
// 					   Rcpp::_["lambdas"]=R_NilValue,
// 					   Rcpp::_["alphas"]=arma::sort(_alphas, "ascend"),
// 					   Rcpp::_["AIC"] = aic,
// 					   Rcpp::_["BIC"] = bic,
// 					   Rcpp::_["loglik"] = loglik,
// 					   Rcpp::_["nParams"] = nParams,
// 					   Rcpp::_["n"] = n);
    
//     result.attr("class") = "graphPath";

//     return result;
// }


// //' Calculate the solution path for an PC graph on a dataset
// //'
// //' @param data The dataframe
// //' @param lambdas A range of lambda values used to calculate a solution path for MGM. If NULL, lambdas is set to nLambda logarithmically spaced values from 10*sqrt(log10(p)/n) to sqrt(log10(p)/n). Defaults to NULL.
// //' @param nLambda The number of lambda values to fit an MGM for when lambdas is NULL
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated MGM graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' g <- rCausalMGM::pcPath(data.n100.p25)
// // [[Rcpp::export]]
// Rcpp::List pcCV(
//     const Rcpp::DataFrame& data,
//     Rcpp::NumericVector alphas = Rcpp::NumericVector::create(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1),
//     Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("max", "conservative", "majority", "none"),
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }

//     int n = ds.getNumRows();
//     int p = ds.getNumColumns();

//     std::vector<double> a(alphas.begin(), alphas.end());

//     arma::vec _alphas(a);

//     _alphas = arma::sort(_alphas, "ascend");

//     arma::vec loglik(_alphas.size(), arma::fill::zeros);
//     arma::vec nParams(_alphas.size(), arma::fill::zeros);
//     std::vector<EdgeListGraph> pcGraphs;
//     std::vector<CausalMGMParams> pcParams;
//     std::vector<double> l = { 0.5, 0.5, 0.5 };

//     for (arma::uword i = 0; i < _alphas.n_elem; i++) {

// 	IndTestMultiCox itm(ds, _alphas(i));

// 	if (orientRule[0]=="none") {

// 	    PcStable pcs((IndependenceTest*) &itm);
// 	    if (threads > 0) pcs.setThreads(threads);
// 	    pcs.setVerbose(verbose);
// 	    pcs.setFDR(fdr);
// 	    pcGraphs.push_back(pcs.search());

// 	} else if (orientRule[0]=="max") {

// 	    PcMax pcm((IndependenceTest*) &itm);
// 	    if (threads > 0) pcm.setThreads(threads);
// 	    pcm.setVerbose(verbose);
// 	    pcm.setFDR(fdr);
// 	    pcGraphs.push_back(pcm.search());

// 	} else if (orientRule[0]=="conservative") {

// 	    CpcStable cpc((IndependenceTest*) &itm);
// 	    if (threads > 0) cpc.setThreads(threads);
// 	    cpc.setVerbose(verbose);
// 	    cpc.setFDR(fdr);
// 	    pcGraphs.push_back(cpc.search());

// 	} else if (orientRule[0]=="majority") {

// 	    Pc50 pc50((IndependenceTest*) &itm);
// 	    if (threads > 0) pc50.setThreads(threads);
// 	    pc50.setVerbose(verbose);
// 	    pc50.setFDR(fdr);
// 	    pcGraphs.push_back(pc50.search());

// 	}

// 	CausalMGM causalMGM(ds, pcGraphs.at(i));
// 	causalMGM.setVerbose(verbose);
// 	causalMGM.setLambda(l);
// 	pcParams.push_back(causalMGM.search());

// 	arma::vec par(pcParams.at(i).toMatrix1D());
// 	loglik(i) = -n * causalMGM.smoothValue(par);
// 	nParams(i) = causalMGM.getNParams();
	
// 	RcppThread::checkUserInterrupt();
//     }

//     // auto elapsedTime = mgm.getElapsedTime();

//     // if (v) {
//     // 	if (elapsedTime < 100*1000) {
//     // 	    Rcpp::Rcout.precision(2);
//     // 	} else {
//     // 	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
//     // 	}
//     //     Rcpp::Rcout << "MGM Path Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
//     // }

//     Rcpp::List graphList;
    
//     for (int i = 0; i < a.size(); i++) {
// 	Rcpp::List pcGraph = pcGraphs[i].toList();
// 	pcGraph["parameters"] = pcParams[i].toList();
//         graphList.push_back(pcGraph);
//     }

//     arma::vec aic = 2*nParams - 2*loglik;
//     arma::vec bic = std::log(n)*nParams - 2*loglik;

//     arma::uword aicIdx = arma::index_min(aic);
//     arma::uword bicIdx = arma::index_min(bic);
    
//     // Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.bic"]=graphList[bicIdx],
//     // 					   Rcpp::_["graph.aic"]=graphList[aicIdx],
//     // 					   Rcpp::_["graphs"]=graphList,
//     // 					   Rcpp::_["lambdas"]=R_NilValue,
//     // 					   Rcpp::_["alphas"]=arma::sort(_alphas, "ascend"),
//     // 					   Rcpp::_["AIC"] = aic,
//     // 					   Rcpp::_["BIC"] = bic,
//     // 					   Rcpp::_["loglik"] = loglik,
//     // 					   Rcpp::_["nParams"] = nParams,
//     // 					   Rcpp::_["n"] = n);

//     Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.min"]=cvGraphs[1].toList(),
// 					   Rcpp::_["graph.1se"]=cvGraphs[0].toList(),
//     					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
// 					   Rcpp::_["lambda.min"]=cvGraphs[1].getHyperParam("lambda"),
// 					   Rcpp::_["lambda.1se"]=cvGraphs[0].getHyperParam("lambda"),
// 					   Rcpp::_["alphas"]=R_NilValue,
// 					   Rcpp::_["alpha.min"]=R_NilValue,
// 					   Rcpp::_["alpha.1se"]=R_NilValue,
// 					   Rcpp::_["foldid"]=_foldid,
//     					   Rcpp::_["loglik"] = loglik);
    
//     result.attr("class") = "graphCV";

//     return result;
// }


// //' Runs the causal algorithm CPC-Stable on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::cpcStable(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List cpcStable(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     const double alpha = 0.05,
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }
//     // bool v = Rcpp::is_true(Rcpp::all(verbose));
//     // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));

//     IndTestMultiCox itm(ds, alpha);

//     CpcStable cpc((IndependenceTest*) &itm);
//     if (threads > 0) cpc.setThreads(threads);
//     cpc.setVerbose(verbose);
//     cpc.setFDR(fdr);
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         cpc.setInitialGraph(&ig);
//     }
//     Rcpp::List result = cpc.search().toList();

//     // ds.deleteVariables();

//     return result;
// }

// //' Runs the causal algorithm PC-Max on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::pcMax(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List pcMax(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     const double alpha = 0.05,
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }

//     IndTestMultiCox itm(ds, alpha);

//     PcMax pcm((IndependenceTest*) &itm);
//     if (threads > 0) pcm.setThreads(threads);
//     pcm.setVerbose(verbose);
//     pcm.setFDR(fdr);
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         pcm.setInitialGraph(&ig);
//     }
//     Rcpp::List result = pcm.search().toList();

//     return result;
// }


// //' Runs the causal algorithm PC50 on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::pc50(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List pc50(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     const double alpha = 0.05,
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }

//     // bool v = Rcpp::is_true(Rcpp::all(verbose));
//     // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));

//     IndTestMultiCox itm(ds, alpha);

//     Pc50 pc50((IndependenceTest*) &itm);
//     if (threads > 0) pc50.setThreads(threads);
//     pc50.setVerbose(verbose);
//     pc50.setFDR(fdr);
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         pc50.setInitialGraph(&ig);
//     }
//     Rcpp::List result = pc50.search().toList();

//     // ds.deleteVariables();

//     return result;
// }



//' Runs the causal algorithm FCI-Stable on a dataset
//'
//' @param data The dataframe
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::fciStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List fciStable(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority"),
    const double alpha = 0.05,
    const int threads = -1,
    const bool fdr = false,
    const bool rank = false,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    
    IndTestMultiCox itm(ds, alpha);
    
    Fci fci((IndependenceTest*) &itm);
    if (threads > 0) fci.setThreads(threads);
    fci.setVerbose(verbose);
    fci.setFDR(fdr);
    fci.setOrientRule(SepsetProducer::str2rule(_orientRule.at(0)));
    
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        fci.setInitialGraph(&ig);
    }
    
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	fci.setKnowledge(k);
    }

    EdgeListGraph g = fci.search();
    Rcpp::List result;
    if (_orientRule.size() == 1) {
	result = g.toList();
    } else {
	bool first = true;
	for (std::string rule : _orientRule) {
	    if (first) {
		first = false;
	    } else {
		g = fci.reorientWithRule(SepsetProducer::str2rule(rule));
	    }
	    result.push_back(g.toList(), rule);
	}
    }
    
    return result;
}


// //' Runs the causal algorithm CFCI-Stable on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::cfci(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List cfci(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     const double alpha = 0.05,
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }
    
//     // bool v = Rcpp::is_true(Rcpp::all(verbose));
//     // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    
//     IndTestMultiCox itm(ds, alpha);
    
//     Cfci cfci((IndependenceTest*) &itm);
//     if (threads > 0) cfci.setThreads(threads);
//     cfci.setVerbose(verbose);
//     cfci.setFDR(fdr);
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         cfci.setInitialGraph(&ig);
//     }
    
//     Rcpp::List result = cfci.search().toList();
    
//     // ds.deleteVariables();
    
//     return result;
// }


// //' Runs the causal algorithm FCI-Max on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::fciMax(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List fciMax(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     const double alpha = 0.05,
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }
    
//     // bool v = Rcpp::is_true(Rcpp::all(verbose));
//     // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    
//     IndTestMultiCox itm(ds, alpha);
    
//     FciMax fcimax((IndependenceTest*) &itm);
//     if (threads > 0) fcimax.setThreads(threads);
//     fcimax.setVerbose(verbose);
//     fcimax.setFDR(fdr);
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         fcimax.setInitialGraph(&ig);
//     }

//     Rcpp::List result;
//     try {
// 	result = fcimax.search().toList();
//     } catch(std::exception& e) {
// 	Rcpp::Rcout << e.what() << std::endl;
//     }
        
//     return result;
// }


// //' Runs the causal algorithm FCI50 Stable on a dataset
// //'
// //' @param data The dataframe
// //' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
// //' @param alpha The p value below which results are considered significant. Defaults to 0.05.
// //' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
// //' @param fdr Whether or not to run with FDR correction for the adjacencies.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::fci50(data.n100.p25, initialGraph = ig)
// // [[Rcpp::export]]
// Rcpp::List fci50(
//     const Rcpp::DataFrame &data,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
//     const double alpha = 0.05,
//     const int threads = -1,
//     const bool fdr = false,
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }
//     // bool v = Rcpp::is_true(Rcpp::all(verbose));
//     // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    
//     IndTestMultiCox itm(ds, alpha);
    
//     Fci50 fci50((IndependenceTest*) &itm);
//     if (threads > 0) fci50.setThreads(threads);
//     fci50.setVerbose(verbose);
//     fci50.setFDR(fdr);
//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         fci50.setInitialGraph(&ig);
//     }
//     Knowledge k;
//     if (!knowledge.isNull()) {
// 	Rcpp::List _knowledge(knowledge);
// 	k = Knowledge(ds.getVariables(), _knowledge);
// 	fci50.setKnowledge(k);
//     }
    
//     Rcpp::List result = fci50.search().toList();
    
//     // ds.deleteVariables();
    
//     return result;
// }


// // no export // [[Rcpp::export]]
// Rcpp::List stars(
//     const Rcpp::DataFrame& data,
//     const std::string method, 
//     Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue,
//     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
//     const int maxDiscrete = 5,
//     const double g = 0.05,
//     const int numSub = 20,
//     const bool adjacency = true,
//     const bool leaveOneOut = false,
//     const int threads = -1,
//     const bool verbose = false
//     ) {

//     // Rcpp::Rcout << "running stars...\n";

//     std::string alg, _method;
//     std::vector<double> par;
//     bool adj = Rcpp::is_true(Rcpp::all(adjacency));
//     bool loo = Rcpp::is_true(Rcpp::all(leaveOneOut));
//     // bool cs = Rcpp::is_true(Rcpp::all(computeStabs));
//     bool v = Rcpp::is_true(Rcpp::all(verbose));

//     _method = method;

//     Rcpp::Rcout << _method << std::endl;
//     std::transform(_method.begin(), _method.end(), _method.begin(),
// 		   [](unsigned char c){ return std::tolower(c); });
//     Rcpp::Rcout << _method << std::endl;

//     if (_method == "mgm") {
// 	alg = "mgm";
//     } else if (_method == "pc" || _method == "pcs" || _method == "pcstable") {
// 	alg = "pc";
//     } else if (_method == "cpc" || _method == "cpcstable") {
// 	alg = "cpc";
//     } else if (_method == "pcm" || _method == "pcmax") {
// 	alg = "pcm";
//     } else if (_method == "fci" || _method == "fcistable") {
// 	alg = "fci";
//     } else if (_method == "cfci" || _method == "cfcistable") {
// 	alg = "cfci";
//     } else if (_method == "fcim" || _method == "fcimax") {
// 	alg = "fcim";
//     } else {
// 	throw std::invalid_argument("Invalid algorithm: " + _method
// 				    + "\n   Algorithm must be in the list: "
// 				    + "{ mgm, pc, cpc, pcm, fci, cfci, fcim }");
//     }
    
//     DataSet ds(data, maxDiscrete);

//     if (params.isNotNull()) {
//         Rcpp::NumericVector _params(params); 
// 	par = std::vector<double>(_params.begin(), _params.end());
//     } else {
// 	if (alg == "mgm") {
// 	    if (ds.getNumRows() > ds.getNumColumns()) {
// 		arma::vec _params = arma::logspace(std::log10(0.05), std::log10(0.8), 20); 
// 		par = std::vector<double>(_params.begin(), _params.end());
// 	    } else {
// 		arma::vec _params = arma::logspace(std::log10(0.1), std::log10(0.8), 20); 
// 		par = std::vector<double>(_params.begin(), _params.end());
// 	    }
// 	} else {
// 	    par = { 0.001, 0.005, 0.01, 0.05, 0.1 };
// 	}
//     }

//     // Rcpp::Rcout << "params vector filled\n";
//     STARS stars(ds, alg, par, g, numSub, adj, loo);
//     // Rcpp::Rcout << "stars object created\n";
//     if (threads > 0) stars.setThreads(threads);
//     // stars.setComputeStabs(cs);
//     stars.setVerbose(v);

//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         stars.setInitialGraph(&ig);
//     }

//     Rcpp::List result = stars.runStarsPar().toList();

//     // if (cs) {
//     //     result["stabilities"] = steps.getStabs();
//     //     std::vector<std::string> names = ds.getVariableNames();
//     //     Rcpp::rownames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
//     //     Rcpp::colnames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
//     // } 

//     // ds.deleteVariables();

//     return result;

// }


// Rcpp::List bootstrap(
//     const Rcpp::DataFrame& data,
//     Rcpp::StringVector algorithm = Rcpp::CharacterVector::create("mgm", "pc", "fci", "mgm-pc", "mgm-fci"),
//     Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
//     const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp",
// 									"conservative", "sepsets"),
//     Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
//     const double alpha = 0.05,
//     const int numBoots = 20,
//     const int threads = -1,
//     const bool replace = true,
//     const bool rank = false,
//     const bool verbose = false
//     )

//' Runs bootstrapping for a selected causal discovery algorithm on the dataset.
//'
//' @param data The dataframe
//' @param algorithm string indicating the name of the causal discovery algorithm to bootstrap. Causal discovery algorithms can be run alone or with mgm to learn an initial graph. Options include "mgm", "pc", "cpc", "pcm", "pc50", "fci", "cfci", "fcim", "mgm-pc", "mgm-cpc", "mgm-pcm", "mgm-pc50", "mgm-fci", "mgm-cfci", "mgm-fcim", "mgm-fci50." The default value is set to "mgm-pc50."
//' @param ensemble Method for construncting an ensemble graph from bootstrapped graphs. Options include "highest", which orients edges according to the orientation with the highest bootstrap probability, or "majority", which only orients edges if they have an orientation with a bootstrap probability > 0.5. Otherwise, the adjacency is included but the edge is left unoreineted. The default value is "highest."
//' @param lambda A vector of three lambda values - the first for continuous-continuous interaction, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param numBoots The number of bootstrap samples to run. Defaults to 20.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph with a table of edge stabilities
//' @export
//' @examples
//' data("data.n100.p25")
//' g.boot <- rCausalMGM::bootstrap(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List bootstrap(
    const Rcpp::DataFrame& data,
    Rcpp::List graph,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const int numBoots = 20,
    const int threads = -1,
    const bool replace = false,
    const bool rank = false,
    const bool verbose = false
    ) {

    // Rcpp::Rcout << "running bootstrap...\n";

    EdgeListGraph g(graph);

    std::string alg, _method, rule;

    _method = g.getAlgorithm();

    // Rcpp::Rcout << _method << std::endl;
    std::transform(_method.begin(), _method.end(), _method.begin(),
		   [](unsigned char c){ return std::tolower(c); });
    // Rcpp::Rcout << _method << std::endl;
    _method.erase(std::remove(_method.begin(), _method.end(), '-'), _method.end());

    std::string tempStr = "mgm";
    std::string::size_type found = _method.find(tempStr);
    
    if (found != std::string::npos) {
	if (_method[found + tempStr.length()] == 'm') {
	    rule = "majority";
	    _method.erase(found + tempStr.length(), 1);
	} else if (_method[found + tempStr.length()] == 'c') {
	    rule = "conservative";
	    _method.erase(found + tempStr.length(), 1);
	} else {
	    rule = "sepsets";
	}
    } else {
	if (_method[0] == 'm') {
	    rule = "majority";
	    _method.erase(0, 1);
	} else if (_method[0] == 'c') {
	    rule = "conservative";
	    _method.erase(0, 1);
	} else {
	    rule = "sepsets";
	}
    }

    tempStr = "max";
    found = _method.find(tempStr);
    if (found != std::string::npos) {
	rule = "maxp";
	_method.erase(found, tempStr.length());
    }
    
    tempStr = "stable";
    found = _method.find(tempStr);
    if (found != std::string::npos)
	_method.erase(found, tempStr.length());
    
    // Rcpp::Rcout << _method << std::endl;

    // Rcpp::Rcout << rule << std::endl;

    if (_method == "mgm" || _method == "coxmgm") {
	alg = "mgm";
    } else if (_method.substr(0,2) == "pc") {
	alg = "pc";
    } else if (_method.substr(0,3) == "fci") {
	alg = "fci";
    } else if (_method.substr(0,5) == "mgmpc" || _method.substr(0,8) == "coxmgmpc") {
	alg = "mgmpc";
    } else if (_method.substr(0,6) == "mgmfci" || _method.substr(0,9) == "coxmgmfci") {
	alg = "mgmfci";
    } else {
	throw std::invalid_argument("Invalid algorithm: " + _method
				    + "\n   Algorithm must be in the list: "
				    + "{ mgm, coxmgm, pc, fci, mgmpc, mgmfci, coxmgmpc, coxmgmfci }");
    }
    
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
    	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
    	ds.npnTransform();
    	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    arma::vec _lambda = g.getHyperParam("lambda");
    std::vector<double> l;

    if (_lambda.n_elem>0) {
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    }

    arma::vec _alpha = g.getHyperParam("alpha");
    double alpha = 0.05;
    
    if (_alpha.n_elem > 0) {
	alpha = _alpha(0);
    }

    Bootstrap boot(ds, alg, numBoots, replace);
    if (threads > 0) boot.setThreads(threads);

    Knowledge k;
    if (!knowledge.isNull()) {
    	Rcpp::List _knowledge(knowledge);
    	k = Knowledge(ds.getVariables(), _knowledge);
    	boot.setKnowledge(k);
    }

    boot.setVerbose(verbose);
    boot.setAlpha(alpha);
    boot.setLambda(l);
    boot.setOrientRule(SepsetProducer::str2rule(rule));
    
    // boot.setFdr(fdr);

    Rcpp::List result = boot.runBootstrap().toList();

    result["stabilities"] = boot.getStabs();

    result["subsamples"] = boot.getSubSamps();

    // result.push_back(boot.getSubSamps(), "subsamples");

    result.attr("class") = "graph";

    return result;

}


//' Runs the Grow-Shrink algorithm to find the Markov blanket of a feature in a dataset
//'
//' @param data The dataframe
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The list of features in the Markov Blanket and the BIC score
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::growShrinkMB(data.n100.p25)
// [[Rcpp::export]]
Rcpp::StringVector growShrinkMB(
    const Rcpp::DataFrame& data,
    const std::string& target,
    const double penalty = 1,
    const bool rank = false,
    const bool verbose = false
) {
    // Rcpp::Nullable<Rcpp::List> graph = R_NilValue,
   
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    // EdgeListGraph g;
    // if (!graph.isNull()) {
    //     Rcpp::List _graph(graph);
    //     g = EdgeListGraph(_graph, ds);
    //     gs.setGraph(g);
    // }

    Node targetNode = ds.getVariable(target);

    std::vector<Node> regressors(ds.getVariables());
    auto it = std::remove(regressors.begin(), regressors.end(), targetNode);

    regressors.erase(it, regressors.end());

    double score = 0;
    
    std::vector<Node> mb;

    if (!ds.isCensored()) {
	DegenerateGaussianScore scorer(ds, penalty);
	GrowShrink gs(&scorer);
	gs.setVerbose(verbose);
	
	mb = gs.search(targetNode, regressors, &score);
    } else {
        RegressionBicScore scorer(ds, penalty);
	GrowShrink gs(&scorer);
	gs.setVerbose(verbose);

	mb = gs.search(targetNode, regressors, &score);
    }
    
    RcppThread::checkUserInterrupt();

    Rcpp::StringVector _mb;

    for (Node n : mb) {
	_mb.push_back(n.getName());
    }

    // Rcpp::List result = Rcpp::List::create(Rcpp::_["markov.blanket"]=_mb,
    // 					   Rcpp::_["SCORE"]=score);
    
    _mb.attr("Score") = Rcpp::wrap(score);
    
    return _mb;
}


//' Runs the GRaSP causal discovery algorithm on the dataset 
//'
//' @param data The dataframe
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The list of features in the Markov Blanket and the BIC score
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::markovBlanket(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List grasp(
    const Rcpp::DataFrame& data,
    const int depth = 2,
    const int numStarts = 1,
    const double penalty = 2,
    const bool rank = false,
    const int threads = -1,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    Score *scorer;

    if (!ds.isCensored()) {
	scorer = new DegenerateGaussianScore(ds, penalty);
    } else {
	throw std::runtime_error("BOSS is not able to handle censored variables.");
        // scorer = new RegressionBicScore(ds, penalty);
    }	

    Grasp grasp(scorer, threads);
    grasp.setVerbose(verbose);
    grasp.setDepth(depth);
    grasp.setNumStarts(numStarts);
    // grasp.setPenalty(penalty);

    // EdgeListGraph ig;
    // if (!initialGraph.isNull()) {
    //     Rcpp::List _initialGraph(initialGraph);
    //     ig = EdgeListGraph(_initialGraph, ds);
    //     grasp.setInitialGraph(&ig);
    // }
    // Knowledge k;
    // if (!knowledge.isNull()) {
    // 	Rcpp::List _knowledge(knowledge);
    // 	k = Knowledge(ds.getVariables(), _knowledge);
    // 	grasp.setKnowledge(k);
    // }

    RcppThread::checkUserInterrupt();

    // std::map<EdgeListGraph, std::pair<int, double>> cpdagMap = grasp.search();

    // Rcpp::List graphList;
    // std::vector<int> graphCounts;
    // std::vector<double> graphBICs;

    // for(auto it = cpdagMap.begin(); it != cpdagMap.end(); ++it) {
    // 	graphList.push_back(it->first.toList());
    // 	graphCounts.push_back(it->second.first);
    // 	graphBICs.push_back(it->second.second);
    // }
    
    // Rcpp::List result = Rcpp::List::create(Rcpp::_["graphs"]=graphList,
    // 					   Rcpp::_["count"]=graphCounts,
    // 					   Rcpp::_["BIC"]=graphBICs
    // 	);

    EdgeListGraph g = grasp.search();

    Rcpp::List result = g.toList();

    double score = g.getScore();

    result.attr("Score") = Rcpp::wrap(score);

    delete scorer;

    return result;
}

//' Runs the BOSS causal discovery algorithm on the dataset 
//'
//' @param data The dataframe
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The list of features in the Markov Blanket and the BIC score
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::markovBlanket(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List boss(
    const Rcpp::DataFrame& data,
    const int numStarts = 1,
    const double penalty = 2,
    const bool rank = false,
    const int threads = -1,
    const bool verbose = false
) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int varIdx = StabilityUtils::checkForVariance(ds);
    if (varIdx >= 0) {
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category, Censored: <10 events).");
    }

    Score *scorer;

    if (!ds.isCensored()) {
	scorer = new DegenerateGaussianScore(ds, penalty);
    } else {
	throw std::runtime_error("BOSS is not able to handle censored variables.");
        // scorer = new RegressionBicScore(ds, penalty);
    }	

    Boss boss(scorer, threads);
    boss.setVerbose(verbose);
    boss.setNumStarts(numStarts);
    // boss.setPenalty(penalty);

    // EdgeListGraph ig;
    // if (!initialGraph.isNull()) {
    //     Rcpp::List _initialGraph(initialGraph);
    //     ig = EdgeListGraph(_initialGraph, ds);
    //     boss.setInitialGraph(&ig);
    // }
    // Knowledge k;
    // if (!knowledge.isNull()) {
    // 	Rcpp::List _knowledge(knowledge);
    // 	k = Knowledge(ds.getVariables(), _knowledge);
    // 	boss.setKnowledge(k);
    // }

    RcppThread::checkUserInterrupt();

    EdgeListGraph g = boss.search();

    Rcpp::List result = g.toList();

    double score = g.getScore();

    result.attr("Score") = Rcpp::wrap(score);

    delete scorer;

    return result;
}



// //' Parameterize a graph using the MGM framework
// //'
// //' @param data The dataframe
// //' @param graph An graph to parameterize.
// //' @param rank Whether or not to use rank-based associations as opposed to linear
// //' @param verbose Whether or not to output additional information. Defaults to FALSE.
// //' @return The calculated search graph
// //' @export
// //' @examples
// //' data("data.n100.p25")
// //' ig <- rCausalMGM::mgm(data.n100.p25)
// //' g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
// //' g <- parameterize(data.n100.p25, g)
// // [[Rcpp::export]]
// Rcpp::List parameterize(
//     const Rcpp::DataFrame& data,
//     Rcpp::List graph,
//     Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.5, 0.5, 0.5),
//     const bool rank = false,
//     const bool verbose = false
// ) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     if (rank) {
// 	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
// 	ds.npnTransform();
// 	if (verbose) Rcpp::Rcout << "done\n";
//     }

//     std::vector<double> l(lambda.begin(), lambda.end());

//     int lamLength = 3;

//     if (l.size() == 1) {
// 	for (int i = 1; i < lamLength; i++) {
// 	    l.push_back(l[0]);
// 	}
//     } else if (l.size() != lamLength) {
// 	throw std::runtime_error("The regularization parameter lambda should be either a vector of length " + std::to_string(lamLength) + " or a single value for this dataset.");
//     }

//     EdgeListGraph g(graph, ds);
//     CausalMGM causalMGM(ds, g);
//     causalMGM.setVerbose(verbose);
//     causalMGM.setLambda(l);

//     Rcpp::List result = g.toList();

//     result["parameters"] = causalMGM.search().toList();

//     // ds.deleteVariables();

//     return result;
// }
