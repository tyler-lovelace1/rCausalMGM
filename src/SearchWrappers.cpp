#include "MGM.hpp"
#include "CoxMGM.hpp"
#include "PcStable.hpp"
// #include "CpcStable.hpp"
// #include "PcMax.hpp"
// #include "Pc50.hpp"
#include "Fci.hpp"
// #include "Cfci.hpp"
// #include "FciMax.hpp"
// #include "Fci50.hpp"
#include "STEPS.hpp"
#include "StabilityUtils.hpp"
// #include "CausalMGMParams.hpp"
// #include "STARS.hpp"
#include "Bootstrap.hpp"
// #include "Tests.hpp"
#include "IndTestMultiCox.hpp"
// #include "BayesIndTestMultiCox.hpp"
#include "GrowShrink.hpp"
#include "DegenerateGaussianScore.hpp"
#include "RegressionBicScore.hpp"
#include "Grasp.hpp"
// #include "CausalMGM.hpp"
#include "Knowledge.hpp"

//' Calculate the MGM graph on a dataset
//'
//' @param data The dataframe
//' @param lambda A vector of three lambda values - the first for continuous-continuous interaction, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param rank Whether or not to use rank-based associations as opposed to linear
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::mgm(data.n100.p25)
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
    // mgm.calcLambdaMax();
    EdgeListGraph mgmGraph = mgm.search();
    
    // mgmGraph.setHyperParam("lambda", Rcpp::NumericVector(l.begin(), l.end()));

    RcppThread::checkUserInterrupt();

    auto elapsedTime = mgm.getElapsedTime();

    if (verbose) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "MGM Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

    Rcpp::List result = mgmGraph.toList();

    return result;
}


//' Calculate the CoxMGM graph on a dataset
//'
//' @param data The dataframe
//' @param lambda A vector of three lambda values - the first for continuous-continuous interaction, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::coxmgm(data.n100.p25)
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
    
    bool v = verbose; // Rcpp::is_true(Rcpp::all(verbose));

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
    
    // mgmGraph.setHyperParam("lambda", Rcpp::NumericVector(l.begin(), l.end()));

    RcppThread::checkUserInterrupt();

    auto elapsedTime = mgm.getElapsedTime();

    if (v) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "CoxMGM Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

    Rcpp::List result = mgmGraph.toList();

    return result;
}


//' Calculate the solution path for an MGM graph on a dataset
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
//' g <- rCausalMGM::mgmPath(data.n100.p25)
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

    bool v = verbose; // Rcpp::is_true(Rcpp::all(verbose));

    arma::vec _lambda;
    std::vector<double> l; // = {0.2, 0.2, 0.2};

    int n = ds.getNumRows();
    int p = ds.getNumColumns();

    MGM mgm(ds);
    mgm.setVerbose(v);

    double logLambdaMax = std::log10(mgm.calcLambdaMax());

    logLambdaMax = std::min(logLambdaMax,
    			    std::log10(10 * std::sqrt(std::log10(ds.getNumColumns()) /
						      ((double) ds.getNumRows()))));

    if (lambdas.isNotNull()) {
        _lambda = arma::vec(Rcpp::as<arma::vec>(lambdas)); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
	if (ds.getNumRows() > ds.getNumColumns()) {
	    _lambda = arma::logspace(logLambdaMax-2, logLambdaMax, nLambda); 
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

    auto elapsedTime = mgm.getElapsedTime();

    if (v) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "MGM Path Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

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
//' @param data The dataframe
//' @param lambdas A range of lambda values used to calculate a solution path for MGM. If NULL, lambdas is set to nLambda logarithmically spaced values from 10*sqrt(log10(p)/n) to sqrt(log10(p)/n). Defaults to NULL.
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::coxmgmPath(data.n100.p25)
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
	    _lambda = arma::logspace(logLambdaMax-2, logLambdaMax, nLambda); 
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

    auto elapsedTime = mgm.getElapsedTime();

    if (v) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "CoxMGM Path Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

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
    const int nfolds = 10,
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
	    _lambda = arma::logspace(logLambdaMax-2, logLambdaMax, nLambda); 
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
	arma::vec folds = arma::linspace(1, nfolds, nfolds);
	folds = arma::repmat(folds, 1, std::ceil(n/((double)nfolds)));
	_foldid = arma::conv_to<arma::uvec>::from(folds(Rcpp::as<arma::uvec>(Rcpp::sample(n,n))-1));
    } else {
	_foldid = Rcpp::as<arma::uvec>(foldid);
    }

    if (_foldid.n_elem != n) {
	throw std::invalid_argument("foldid has a length that does not equal number of samples in the dataset.");
    }

    // Rcpp::Rcout << "Fold ids: " << _foldid.t() << std::endl;
    // Rcpp::Rcout << "Beginning path search for lambdas " << _lambda.t() << std::endl;
    arma::mat loglik(l.size(), nfolds, arma::fill::zeros);
    // arma::vec nParams(l.size(), arma::fill::zeros);
    std::vector<EdgeListGraph> cvGraphs = mgm.searchPathCV(l, nfolds, _foldid, loglik);

    RcppThread::checkUserInterrupt();

    auto elapsedTime = mgm.getElapsedTime();

    if (v) {
	if (elapsedTime < 100*1000) {
	    Rcpp::Rcout.precision(2);
	} else {
	    elapsedTime = std::round(elapsedTime / 1000.0) * 1000;
	}
        Rcpp::Rcout << "MGM Cross-Validation Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

    // Rcpp::List graphList;
    
    // for (int i = 0; i < l.size(); i++) {
    //     graphList.push_back(mgmGraphs[i].toList());
    // }

    // Rcpp::List result;
    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.min"]=cvGraphs[1].toList(),
					   Rcpp::_["graph.1se"]=cvGraphs[0].toList(),
    					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
					   Rcpp::_["lambda.min"]=cvGraphs[1].getHyperParam("lambda"),
					   Rcpp::_["lambda.1se"]=cvGraphs[0].getHyperParam("lambda"),
					   Rcpp::_["alphas"]=R_NilValue,
					   Rcpp::_["alpha.min"]=R_NilValue,
					   Rcpp::_["alpha.1se"]=R_NilValue,
					   Rcpp::_["foldid"]=_foldid,
    					   Rcpp::_["loglik"] = loglik);

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
	    _lambda = arma::logspace(logLambdaMax-1, logLambdaMax, nLambda); 
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

    Rcpp::List graph = steps.runStepsPath(instabs, samps).toList();

    if (computeStabs) {
        graph["stabilities"] = steps.getStabs();
        std::vector<std::string> names = ds.getVariableNames();
        Rcpp::rownames(graph["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
        Rcpp::colnames(graph["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph"]=graph,
    					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
					   Rcpp::_["gamma"]=gamma,
    					   Rcpp::_["instability"] = instabs,
					   Rcpp::_["subsamples"] = samps);

    result.attr("class") = "graphSTEPS";

    // ds.deleteVariables();

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
    const Rcpp::DataFrame &data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp",
									"conservative", "sepsets"),
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

    std::string rule;
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    // bool v = Rcpp::is_true(Rcpp::all(verbose));

    std::vector<std::string> validNames = {"majority", "maxp", "conservative", "sepsets"};

    rule = orientRule[0];

    std::transform(rule.begin(), rule.end(), rule.begin(),
		   [](unsigned char c){ return std::tolower(c); });

    if (std::find(validNames.begin(), validNames.end(), rule) == validNames.end())
	throw std::invalid_argument("Orientation rule must be one of {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    
    IndTestMultiCox itm(ds, alpha);
    
    PcStable pc((IndependenceTest*) &itm);
    if (threads > 0) pc.setThreads(threads);
    pc.setVerbose(verbose);
    pc.setFDR(fdr);
    if (rule == "majority")
	pc.setOrientRule(ORIENT_MAJORITY);
    if (rule == "maxp")
	pc.setOrientRule(ORIENT_MAXP);
    if (rule == "conservative")
	pc.setOrientRule(ORIENT_CONSERVATIVE);
    if (rule == "sepsets")
	pc.setOrientRule(ORIENT_SEPSETS);
    
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

    Rcpp::List result = pc.search().toList();
    
    // ds.deleteVariables();
    
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
    const Rcpp::DataFrame &data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp",
									"conservative", "sepsets"),
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

    std::string rule;
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    // bool v = Rcpp::is_true(Rcpp::all(verbose));

    std::vector<std::string> validNames = {"majority", "maxp", "conservative", "sepsets"};

    rule = orientRule[0];

    std::transform(rule.begin(), rule.end(), rule.begin(),
		   [](unsigned char c){ return std::tolower(c); });

    if (std::find(validNames.begin(), validNames.end(), rule) == validNames.end())
	throw std::invalid_argument("Orientation rule must be one of {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    
    IndTestMultiCox itm(ds, alpha);
    
    Fci fci((IndependenceTest*) &itm);
    if (threads > 0) fci.setThreads(threads);
    fci.setVerbose(verbose);
    fci.setFDR(fdr);
    if (rule == "majority")
	fci.setOrientRule(ORIENT_MAJORITY);
    if (rule == "maxp")
	fci.setOrientRule(ORIENT_MAXP);
    if (rule == "conservative")
	fci.setOrientRule(ORIENT_CONSERVATIVE);
    if (rule == "sepsets")
	fci.setOrientRule(ORIENT_SEPSETS);
    
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

    Rcpp::List result = fci.search().toList();

    // result.push_back(, "majority");

    // result.push_back(fci.reorientWithRule(ORIENT_MAXP).toList(), "maxp");
    // result.push_back(fci.reorientWithRule(ORIENT_CONSERVATIVE).toList(), "conservative");
    // result.push_back(fci.reorientWithRule(ORIENT_SEPSETS).toList(), "sepsets");
    // fci.reorientWithRule(ORIENT_CONSERVATIVE);
    // fci.reorientWithRule(ORIENT_SEPSETS);
    
    // ds.deleteVariables();
    
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
    Rcpp::StringVector algorithm = Rcpp::CharacterVector::create("mgm", "pc", "fci", "mgm-pc", "mgm-fci"),
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority", "maxp",
									"conservative", "sepsets"),
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
    const double alpha = 0.05,
    const int numBoots = 20,
    const int threads = -1,
    const bool replace = true,
    const bool rank = false,
    const bool verbose = false
    ) {

    // Rcpp::Rcout << "running stars...\n";

    std::string alg, _method;
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    // bool v = Rcpp::is_true(Rcpp::all(verbose));

    _method = algorithm[0];

    // Rcpp::Rcout << _method << std::endl;
    std::transform(_method.begin(), _method.end(), _method.begin(),
		   [](unsigned char c){ return std::tolower(c); });
    // Rcpp::Rcout << _method << std::endl;
    _method.erase(std::remove(_method.begin(), _method.end(), '-'), _method.end());
    // Rcpp::Rcout << _method << std::endl;

    if (_method == "mgm") {
	alg = "mgm";
    } else if (_method.substr(0,2) == "pc") {
	alg = "pc";
    } else if (_method.substr(0,3) == "fci") {
	alg = "fci";
    } else if (_method.substr(0,5) == "mgmpc") {
	alg = "mgmpc";
    } else if (_method.substr(0,6) == "mgmfci") {
	alg = "mgmfci";
    } else {
	throw std::invalid_argument("Invalid algorithm: " + _method
				    + "\n   Algorithm must be in the list: "
				    + "{ mgm, pc, fci, mgm-pc, mgm-fci }");
    }

    std::string rule;

    std::vector<std::string> validNames = {"majority", "maxp", "conservative", "sepsets"};

    rule = orientRule[0];

    std::transform(rule.begin(), rule.end(), rule.begin(),
		   [](unsigned char c){ return std::tolower(c); });

    if (std::find(validNames.begin(), validNames.end(), rule) == validNames.end())
	throw std::invalid_argument("Orientation rule must be one of {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");

    std::vector<double> l(lambda.begin(), lambda.end());

    DataSet ds = DataSet(data);
    ds.dropMissing();

    if (rank) {
	if (verbose) Rcpp::Rcout << "Applying the nonparanormal transform to continuous variables...";
	ds.npnTransform();
	if (verbose) Rcpp::Rcout << "done\n";
    }

    int lamLength = 3;

    if (ds.isMixed()) {
    	if (ds.isCensored()) {
    	    lamLength = 5;
    	} else {
    	    lamLength = 3;
    	}
    	lamLength = 3;
    } // else {
    // 	throw std::runtime_error("MGM is not implemented for purely continuous or purely discrete datasets.");
    // }

    if (l.size() == 1) {
	for (int i = 1; i < lamLength; i++) {
	    l.push_back(l[0]);
	}
    } else if (l.size() != lamLength) {
	throw std::runtime_error("The regularization parameter lambda should be either a vector of length " + std::to_string(lamLength) + " or a single value for this dataset.");
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
    if (rule.substr(0,3) == "maj")
	boot.setOrientRule(ORIENT_MAJORITY);
    if (rule.substr(0,3) == "max")
	boot.setOrientRule(ORIENT_MAXP);
    if (rule.substr(0,1) == "c")
	boot.setOrientRule(ORIENT_CONSERVATIVE);
    if (rule.substr(0,1) == "s")
	boot.setOrientRule(ORIENT_SEPSETS);
    // boot.setFdr(fdr);

    Rcpp::List result = boot.runBootstrap().toList();

    result["stabilities"] = boot.getStabs();

    // result["subsamples"] = boot.getSubSamps();

    // result.push_back(boot.getSubSamps(), "subsamples");

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
    const int threads = -1,
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
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const int depth = 2,
    const int numStarts = 3,
    const double penalty = 1,
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

    Score *scorer;

    if (!ds.isCensored()) {
	scorer = new DegenerateGaussianScore(ds, penalty);
    } else {
        scorer = new RegressionBicScore(ds, penalty);
    }	

    Grasp grasp(scorer, threads);
    grasp.setVerbose(verbose);
    grasp.setDepth(depth);
    grasp.setNumStarts(numStarts);
    grasp.setPenalty(penalty);

    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        grasp.setInitialGraph(&ig);
    }
    Knowledge k;
    if (!knowledge.isNull()) {
	Rcpp::List _knowledge(knowledge);
	k = Knowledge(ds.getVariables(), _knowledge);
	grasp.setKnowledge(k);
    }

    RcppThread::checkUserInterrupt();

    std::map<EdgeListGraph, std::pair<int, double>> cpdagMap = grasp.search();

    Rcpp::List graphList;
    std::vector<int> graphCounts;
    std::vector<double> graphBICs;

    for(auto it = cpdagMap.begin(); it != cpdagMap.end(); ++it) {
	graphList.push_back(it->first.toList());
	graphCounts.push_back(it->second.first);
	graphBICs.push_back(it->second.second);
    }
    
    Rcpp::List result = Rcpp::List::create(Rcpp::_["graphs"]=graphList,
					   Rcpp::_["count"]=graphCounts,
					   Rcpp::_["BIC"]=graphBICs
	);

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
