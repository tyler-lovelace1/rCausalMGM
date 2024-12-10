#include "MGM.hpp"
#include "PcStable.hpp"
#include "Fci.hpp"
#include "STEPS.hpp"
#include "STARS.hpp"
#include "StabilityUtils.hpp"
#include "Bootstrap.hpp"
#include "IndTestMulti.hpp"
#include "SearchCV.hpp"
#include "GrowShrink.hpp"
#include "DegenerateGaussianScore.hpp"
#include "RegressionBicScore.hpp"
#include "Knowledge.hpp"

//' Calculate the MGM graph on a dataset
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param lambda A numeric vector of three values for the regularization parameter lambda: the first for continuous-continuous edges, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print updates on the progress of optimizing MGM. The default is FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- mgm(sim$data)
//' print(g)
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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


//' Estimates a solution path for MGM
//'
//' @description Calculate the solution path for an MGM graph on a dataset. It also returns the models selected by the BIC and AIC scores.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param lambdas A numeric vector containing the values of lambda to learn an MGM with. The default value is NULL, in which case a log-spaced vector of nLambda values for lambda will be supplied instead.
//' @param nLambda A numeric value indicating the number of lambda values to test when the lambdas vector is NULL. The default is 30.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphPath object that contains MGM graphs learned by the solution path, as well as the BIC and AIC selected models
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' ig.path <- mgmPath(sim$data)
//' print(ig.path)
//' }
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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

//' Implements k-fold cross-validation for MGM
//'
//' @description Calculate the solution path for an MGM graph on a dataset with k-fold cross-validation. This function returns the graph that minimizes negative log(pseudolikelihood) and the graph selected by the one standard error rule.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param lambdas A numeric vector containing the values of lambda to learn an MGM with. The default value is NULL, in which case a log-spaced vector of nLambda values for lambda will be supplied instead.
//' @param nLambda A numeric value indicating the number of lambda values to test when the lambdas vector is NULL. The default is 30.
//' @param nfolds An integer value defining the number of folds to be used for cross-validation if foldid is NULL. The default value is 5.
//' @param foldid An integer vector containing values in the range of 1 to K for each sample that identifies which test set that sample belongs to. This enables users to define their own cross-validation splits, for example in the case stratified cross-validation is needed. The default value is NULL.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphCV object that contains the minimum and one standard error rule selected graphs.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' ig.cv <- mgmCV(sim$data)
//' print(ig.cv)
//' }
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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

    arma::mat loglik(l.size(), arma::max(_foldid), arma::fill::zeros);
    arma::uvec index(2, arma::fill::zeros);
    std::vector<EdgeListGraph> cvGraphs = mgm.searchPathCV(l, _foldid, loglik, index);

    RcppThread::checkUserInterrupt();

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



//' Implements StEPS and StARS for MGM
//'
//' @description Calculates the optimal lambda values for the MGM algorithm using StEPS and StARS. Returns a graphSTEPS object that contains the MGMs selected by StEPS and StARS as well as the instability at each value of lambda.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param lambdas A numeric vector containing the values of lambda to learn an MGM with. The default value is NULL, in which case a log-spaced vector of nLambda values for lambda will be supplied instead.
//' @param nLambda A numeric value indicating the number of lambda values to test when the lambdas vector is NULL. The default is 30.
//' @param gamma The threshold for edge instability. The default value is 0.05, and it is not recommended to change this value.
//' @param numSub The number of subsamples of the dataset used to estimate edge instability. The default value is 20.
//' @param subSize The number of samples to be drawn without replacement for each subsample. The default value is -1. When subSize is -1, it is set to min(floor(0.75 * N), floor(10*sqrt(N))), where N is the number of samples.
//' @param leaveOneOut If TRUE, performs leave-one-out subsampling. Defaults to FALSE.
//' @param threads An integer value denoting the number of threads to use for parallelization of learning MGMs across subsamples. The default value is -1, which will all available CPUs.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphSTEPS object containing the MGMs selected by StEPS and StARS, as well as the instability of each edge type at each value of lambda.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' ig.steps <- steps(sim$data)
//' print(ig.steps)
// [[Rcpp::export]]
Rcpp::List steps(
    const Rcpp::DataFrame &data, 
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 30,
    const double gamma = 0.05,
    const int numSub = 20,
    const int subSize = -1,
    const bool leaveOneOut = false,
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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
    // steps.setComputeStabs(false);
    steps.setVerbose(verbose);

    arma::mat instabs;
    // if (ds.isCensored()) {
    // 	instabs = arma::mat(l.size(), 6);
    // } else {
    instabs = arma::mat(l.size(), 4);
    // }
    instabs.fill(arma::datum::nan);

    arma::umat samps;

    std::vector<EdgeListGraph> graphs = steps.runStepsPath(instabs, samps);
    Rcpp::List graphSteps = graphs.at(0).toList();

    // if (computeStabs) {
    //     graphSteps["stabilities"] = steps.getStabs();
    //     std::vector<std::string> names = ds.getVariableNames();
    //     Rcpp::rownames(graphSteps["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    //     Rcpp::colnames(graphSteps["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    // }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph.steps"]=graphSteps,
					   Rcpp::_["graph.stars"]=graphs.at(1).toList(),
    					   Rcpp::_["lambdas"]=arma::sort(_lambda, "descend"),
					   Rcpp::_["gamma"]=gamma,
    					   Rcpp::_["instability"] = instabs,
					   Rcpp::_["subsamples"] = samps);

    result.attr("class") = "graphSTEPS";

    return result;

}


//' Implements StARS for PC-Stable
//'
//' @description Runs StARS to select the value of alpha for PC-Stable based on adjacency stability. Returns a graphSTARS object containing the CPDAG selected by StARS and the adjacency instabilities for each alpha.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param initialGraph An undirected rCausalMGM graph object containing the initial skeleton of adjacencies used in the causal discovery algorithm. This graph can be learned by `mgm` or learned by another method and imported into an undirected rCausalMGM graph object from its adjacency matrix. The default is NULL, in which case a fully connected graph is used as the initial skeleton.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param orientRule Determines which of the four possible orientation rules will be utilized to orient colliders in the PC-Stable algorithm. Possible options are "majority", "maxp", "conservative", and "sepsets". The default value is "majority". 
//' @param alphas A numeric vector containing values of alpha to test in the cross-validation procedure. The default value is NULL, in which case we set alpha = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2).
//' @param gamma The threshold for edge instability. The default value is 0.01, and it is not recommended to change this value.
//' @param numSub The number of subsamples of the dataset used to estimate edge instability. The default value is 20.
//' @param subSize The number of samples to be drawn without replacement for each subsample. The default value is -1. When subSize is -1, it is set to min(floor(0.75 * N), floor(10*sqrt(N))), where N is the number of samples.
//' @param leaveOneOut If TRUE, performs leave-one-out subsampling. Defaults to FALSE.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphSTARS object containing the CPDAG selected by StARS and the instabilities at each value of alpha.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' g.stars <- pcStars(sim$data)
//' print(g.stars)
//' }
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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

    // if (computeStabs) {
    //     graphStars["stabilities"] = stars.getStabs();
    //     std::vector<std::string> names = ds.getVariableNames();
    //     Rcpp::rownames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    //     Rcpp::colnames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    // }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph"]=graphStars,
    					   Rcpp::_["alphas"]=arma::sort(_alphas, "ascend"),
					   Rcpp::_["gamma"]=gamma,
    					   Rcpp::_["instability"] = instabs,
					   Rcpp::_["subsamples"] = samps);

    result.attr("class") = "graphSTARS";

    return result;
}

//' Implements StARS for FCI-Stable
//'
//' @description Runs StARS to select the value of alpha for FCI-Stable based on adjacency stability. Returns a graphSTARS object containing the PAG selected by StARS and the adjacency instabilities for each alpha.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param initialGraph An undirected rCausalMGM graph object containing the initial skeleton of adjacencies used in the causal discovery algorithm. This graph can be learned by `mgm` or learned by another method and imported into an undirected rCausalMGM graph object from its adjacency matrix. The default is NULL, in which case a fully connected graph is used as the initial skeleton.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param orientRule Determines which of the four possible orientation rules will be utilized to orient colliders in the FCI-Stable algorithm. Possible options are "majority", "maxp", "conservative", and "sepsets". The default value is "majority". 
//' @param alphas A numeric vector containing values of alpha to test in the cross-validation procedure. The default value is NULL, in which case we set alpha = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2).
//' @param gamma The threshold for edge instability. The default value is 0.01, and it is not recommended to change this value.
//' @param numSub The number of subsamples of the dataset used to estimate edge instability. The default value is 20.
//' @param subSize The number of samples to be drawn without replacement for each subsample. The default value is -1. When subSize is -1, it is set to min(floor(0.75 * N), floor(10*sqrt(N))), where N is the number of samples.
//' @param leaveOneOut If TRUE, performs leave-one-out subsampling. Defaults to FALSE.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphSTARS object containing the PAG selected by StARS and the instabilities at each value of alpha.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' g.stars <- fciStars(sim$data)
//' print(g.stars)
//' }
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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

    // if (computeStabs) {
    //     graphStars["stabilities"] = stars.getStabs();
    //     std::vector<std::string> names = ds.getVariableNames();
    //     Rcpp::rownames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    //     Rcpp::colnames(graphStars["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    // }

    Rcpp::List result = Rcpp::List::create(Rcpp::_["graph"]=graphStars,
    					   Rcpp::_["alphas"]=arma::sort(_alphas, "ascend"),
					   Rcpp::_["gamma"]=gamma,
    					   Rcpp::_["instability"] = instabs,
					   Rcpp::_["subsamples"] = samps);

    result.attr("class") = "graphSTARS";

    return result;
}


//' Runs the causal discovery algorithm PC-Stable on a dataset.
//'
//' @description Runs the causal discovery algorithm PC-Stable on a dataset. The PC-Stable algorithm is designed to recover the Markov equivalence class of causal DAGs that could give rise to the observed conditional independence relationships under the assumption of causal sufficiency. A dataset is said to be causally sufficient if all variables relevant to the causal process are observed (i.e. there are no latent confounders). The resulting graph is a completed partially directed acyclic graph (CPDAG) containing directed edges where the causal orientation can be uniquely determined and an undirected edge where multiple orientations are possible.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param initialGraph An undirected rCausalMGM graph object containing the initial skeleton of adjacencies used in the causal discovery algorithm. This graph can be learned by `mgm` or learned by another method and imported into an undirected rCausalMGM graph object from its adjacency matrix. The default is NULL, in which case a fully connected graph is used as the initial skeleton.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param orientRule Determines which of the four possible orientation rules will be utilized to orient colliders in the PC-Stable algorithm. Possible options are "majority", "maxp", "conservative", and "sepsets". The default value is "majority". Additionally, a vector of valid orientation rules can be provided, and pcStable will return a list containing the graphs learned with each.
//' @param alpha A numeric value containing the significance threshold alpha for the conditional independence tests used during constraint-based causal discovery. This parameter directly controls graph sparsity, with low values of alpha yielding sparse graphs and high values yielding dense graphs. The default value is 0.05.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param fdr A logical value indicating whether to use false discovery rate control for the discovery of adjacencies in the causal graph. The default value is FALSE.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return The CPDAG learned by PC-Stable.
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' print(g)
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    
    IndTestMulti itm(ds, alpha);
    
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


//' Runs the causal discovery algorithm FCI-Stable on a dataset.
//'
//' @description Runs the causal discovery algorithm FCI-Stable on a dataset. The FCI-Stable algorithm is designed to recover the Markov equivalence class of causal MAGs that could give rise to the observed conditional independence relationships in the causally insufficient case. This means that FCI-Stable can still learn the Markov equivalence class of the true MAG even in the presence of latent confounders and/or selection bias. The resulting graph is a partial ancestral graph (PAG).
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param initialGraph An undirected rCausalMGM graph object containing the initial skeleton of adjacencies used in the causal discovery algorithm. This graph can be learned by `mgm` or learned by another method and imported into an undirected rCausalMGM graph object from its adjacency matrix. The default is NULL, in which case a fully connected graph is used as the initial skeleton.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param orientRule Determines which of the four possible orientation rules will be utilized to orient colliders in the FCI-Stable algorithm. Possible options are "majority", "maxp", "conservative", and "sepsets". The default value is "majority". Additionally, a vector of valid orientation rules can be provided, and fciStable will return a list containing the graphs learned with each.
//' @param alpha A numeric value containing the significance threshold alpha for the conditional independence tests used during constraint-based causal discovery. This parameter directly controls graph sparsity, with low values of alpha yielding sparse graphs and high values yielding dense graphs. The default value is 0.05.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param possDsep A logical value indicating whether to perform the possible-D-Sep search stage of the FCI algorithm. The possible-D-Sep search is necessaey fro correctness but can be computationally expensive in dense or high-dimensional or graphs. If set to FALSE, the RFCI rule R0 will be applied to remove some of the extraneous adjacencies that would have been removed by possible-D-Sep search. The default value is TRUE.
//' @param fdr A logical value indicating whether to use false discovery rate control for the discovery of adjacencies in the causal graph. The default value is FALSE.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return The PAG learned by FCI-Stable.
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' g <- fciStable(sim$data)
//' print(g)
// [[Rcpp::export]]
Rcpp::List fciStable(
    const Rcpp::DataFrame& data,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    Rcpp::Nullable<Rcpp::List> knowledge = R_NilValue,
    const Rcpp::StringVector orientRule = Rcpp::CharacterVector::create("majority"),
    const double alpha = 0.05,
    const int threads = -1,
    const bool possDsep = true,
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
    }

    if (orientRule.size()==0) {
	throw std::invalid_argument("At least one orientation rule must be provided. Options are {\"majority\", \"maxp\", \"conservative\", \"sepsets\"}");
    }
    
    std::vector<std::string> _orientRule = Rcpp::as<std::vector<std::string>>(orientRule);
    
    IndTestMulti itm(ds, alpha);
    
    Fci fci((IndependenceTest*) &itm);
    if (threads > 0) fci.setThreads(threads);
    fci.setVerbose(verbose);
    fci.setFDR(fdr);
    fci.setPossibleDsepSearchDone(possDsep);
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

//' Implements k-fold cross-validation for PC-Stable
//'
//' @description Runs k-fold cross-validation to select the value of alpha and orientation rule for PC-Stable. Returns a graphCV object containing the causal graphical models that minimize the negative log(pseudo-likelihood) and the sparsest model within one standard error of the minimum.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param initialGraph An undirected rCausalMGM graph object containing the initial skeleton of adjacencies used in the causal discovery algorithm. This graph can be learned by `mgm` or learned by another method and imported into an undirected rCausalMGM graph object from its adjacency matrix. The default is NULL, in which case a fully connected graph is used as the initial skeleton.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param orientRule A vector of strings to determine which of the orientation rules to test in the cross-validation procedure to select the optimal model. The default is a vector that contains the "majority", "maxp", and "conservative" orientation rules.
//' @param alphas A numeric vector containing values of alpha to test in the cross-validation procedure. The default value is NULL, in which case we set alpha = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2).
//' @param nfolds An integer value defining the number of folds to be used for cross-validation if foldid is NULL. The default value is 5.
//' @param foldid An integer vector containing values in the range of 1 to K for each sample that identifies which test set that sample belongs to. This enables users to define their own cross-validation splits, for example in the case stratified cross-validation is needed. The default value is NULL.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param fdr A logical value indicating whether to use false discovery rate control for the discovery of adjacencies in the causal graph. The default value is FALSE.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphCV object containing the CPDAGs selected by the minimum and one standard error rule.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' g.cv <- pcCV(sim$data)
//' print(g.cv)
//' }
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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


//' Implements k-fold cross-validation for FCI-Stable
//'
//' @description Runs k-fold cross-validation to select the value of alpha and orientation rule for FCI-Stable. Returns a graphCV object containing the causal graphical models that minimize the negative log(pseudo-likelihood) and the sparsest model within one standard error of the minimum.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param initialGraph An undirected rCausalMGM graph object containing the initial skeleton of adjacencies used in the causal discovery algorithm. This graph can be learned by `mgm` or learned by another method and imported into an undirected rCausalMGM graph object from its adjacency matrix. The default is NULL, in which case a fully connected graph is used as the initial skeleton.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param orientRule A vector of strings to determine which of the orientation rules to test in the cross-validation procedure to select the optimal model. The default is a vector that contains the "majority", "maxp", and "conservative" orientation rules.
//' @param alphas A numeric vector containing values of alpha to test in the cross-validation procedure. The default value is NULL, in which case we set alpha = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2).
//' @param nfolds An integer value defining the number of folds to be used for cross-validation if foldid is NULL. The default value is 5.
//' @param foldid An integer vector containing values in the range of 1 to K for each sample that identifies which test set that sample belongs to. This enables users to define their own cross-validation splits, for example in the case stratified cross-validation is needed. The default value is NULL.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param fdr A logical value indicating whether to use false discovery rate control for the discovery of adjacencies in the causal graph. The default value is FALSE.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphCV object containing the PAGs selected by the minimum and one standard error rule.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' g.cv <- fciCV(sim$data)
//' print(g.cv)
//' }
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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


//' Implements k-fold cross-validation for MGM-PC-Stable
//'
//' @description Runs k-fold cross-validation to select the value of lambda, alpha, and the orientation rule for MGM-PC-Stable. Returns a graphCV object containing the causal graphical models that minimize the negative log(pseudo-likelihood) and the sparsest model within one standard error of the minimum.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param cvType A string determining whether to perform random search or grid search cross-validation, indicated by "random" or "grid" respectively. The default value is "random".
//' @param orientRule A vector of strings to determine which of the orientation rules to test in the cross-validation procedure to select the optimal model. The default is a vector that contains the "majority", "maxp", and "conservative" orientation rules.
//' @param lambdas A numeric vector containing the values of lambda to learn an MGM with. The default value is NULL, in which case a log-spaced vector of nLambda values for lambda will be supplied instead.
//' @param nLambda A numeric value indicating the number of lambda values to test when the lambdas vector is NULL. The default is 20.
//' @param alphas A numeric vector containing values of alpha to test in the cross-validation procedure. The default value is NULL, in which case we set alpha = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2).
//' @param numPoints An integer value containing indicating the number of samples to draw uniformly from the search space if performing random search cross-validation. The default is 60, the number of points required to have a 5\% chance of sampling a model in the top 5\% of the search space.
//' @param nfolds An integer value defining the number of folds to be used for cross-validation if foldid is NULL. The default value is 5.
//' @param foldid An integer vector containing values in the range of 1 to K for each sample that identifies which test set that sample belongs to. This enables users to define their own cross-validation splits, for example in the case stratified cross-validation is needed. The default value is NULL.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param fdr A logical value indicating whether to use false discovery rate control for the discovery of adjacencies in the causal graph. The default value is FALSE.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphCV object containing the CPDAGs selected by the minimum and one standard error rule.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' g.cv <- mgmpcCV(sim$data)
//' print(g.cv)
//' }
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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


//' Implements k-fold cross-validation for MGM-FCI-Stable
//'
//' @description Runs k-fold cross-validation to select the value of lambda, alpha, and the orientation rule for MGM-FCI-Stable. Returns a graphCV object containing the causal graphical models that minimize the negative log(pseudo-likelihood) and the sparsest model within one standard error of the minimum.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param cvType A string determining whether to perform random search or grid search cross-validation, indicated by "random" or "grid" respectively. The default value is "random".
//' @param orientRule A vector of strings to determine which of the orientation rules to test in the cross-validation procedure to select the optimal model. The default is a vector that contains the "majority", "maxp", and "conservative" orientation rules.
//' @param lambdas A numeric vector containing the values of lambda to learn an MGM with. The default value is NULL, in which case a log-spaced vector of nLambda values for lambda will be supplied instead.
//' @param nLambda A numeric value indicating the number of lambda values to test when the lambdas vector is NULL. The default is 20.
//' @param alphas A numeric vector containing values of alpha to test in the cross-validation procedure. The default value is NULL, in which case we set alpha = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2).
//' @param numPoints An integer value containing indicating the number of samples to draw uniformly from the search space if performing random search cross-validation. The default is 60, the number of points required to have a 5\% chance of sampling a model in the top 5\% of the search space.
//' @param nfolds An integer value defining the number of folds to be used for cross-validation if foldid is NULL. The default value is 5.
//' @param foldid An integer vector containing values in the range of 1 to K for each sample that identifies which test set that sample belongs to. This enables users to define their own cross-validation splits, for example in the case stratified cross-validation is needed. The default value is NULL.
//' @param threads An integer value denoting the number of threads to use for parallelization of independence tests. The default value is -1, which will all available CPUs.
//' @param fdr A logical value indicating whether to use false discovery rate control for the discovery of adjacencies in the causal graph. The default value is FALSE.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return A graphCV object containing the PAGs selected by the minimum and one standard error rule.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' g.cv <- mgmfciCV(sim$data)
//' print(g.cv)
//' }
// [[Rcpp::export]]
Rcpp::List mgmfciCV(
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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



//' Runs bootstrapping for a causal graph on the dataset.
//'
//' @description Runs bootstrapping for a causal graph on the dataset. This function can be used to estimate the stability of edge adjacencies and orientations in the causal graph. It returns an ensemble graph which consists of the most common edges accross bootstrap samples. The ensemble graph is constructed based on edge-wise probabilities, so it is not guaranteed to be a valid CPDAG or PAG. The ensemble graph's stabilites entry contains information about the frequency of each possible orientation for each edge that appears at least once across bootstrap samples.
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param graph A graph object containing the graph to estimate the stability of through bootstrapping.
//' @param knowledge A knowledge object containing prior knowledge about the causal interactions in a dataset. This knowledge can be used to forbid or require certain edges in the causal graph, helping to inform causal discovery an prevent orientations known to be nonsensical. The default is NULL, in which case no prior knowledge is provided to the causal discovery algorithm.
//' @param numBoots The number of bootstrap samples to run. The default is 20.
//' @param threads An integer value denoting the number of threads to use for parallelization. The default value is -1, which will all available CPUs.
//' @param replace A logical value indicating whether to use sampling with replacement or to draw subsamples of size floor(0.632 * N). The default value is FALSE.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @export
//' @examples
//' \dontrun{
//' sim <- simRandomDAG(200, 25)
//' g <- pcStable(sim$data)
//' g.boot <- bootstrap(sim$data, g)
//' print(g.boot)
//' print(head(g.boot$stabilities))
//' }
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
				    + "{ mgm, pc, fci, mgmpc, mgmfci }");
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
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


//' Implements Grow-Shrink algorithm for Markov blanket identification
//'
//' @description Runs the Grow-Shrink algorithm to find the Markov blanket of a feature in a dataset
//'
//' @param data A data.frame containing the dataset to be used for estimating the MGM, with each row representing a sample and each column representing a variable. All continuous variables must be of the numeric type, while categorical variables must be factor or character. Any rows with missing values will be dropped.
//' @param target A string denoting the name of the target variable to identify the Markov blanket of.
//' @param penalty A numeric value that represents the strength of the penalty for model complexity. The default value is 1, which corresponds to the BIC score.
//' @param rank A logical value indicating whether to use the nonparanormal transform to learn rank-based associations. The default is FALSE.
//' @param verbose A logical value indicating whether to print progress updates. The default is FALSE.
//' @return The list of features in the Markov Blanket and the BIC score
//' @export
//' @examples
//' sim <- simRandomDAG(200, 25)
//' mb <- growShrinkMB(sim$data, "X1")
//' print(mb)
// [[Rcpp::export]]
Rcpp::StringVector growShrinkMB(
    const Rcpp::DataFrame& data,
    const std::string& target,
    const double penalty = 1,
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
	throw std::invalid_argument("The variable " + ds.getVariable(varIdx).getName() + " has an invalid variance (Continuous: all values are the same, Categorical: <5 samples in a category).");
    }

    Node targetNode = ds.getVariable(target);

    std::vector<Node> regressors(ds.getVariables());
    auto it = std::remove(regressors.begin(), regressors.end(), targetNode);

    regressors.erase(it, regressors.end());

    double score = 0;
    
    std::vector<Node> mb;

    // if (!ds.isCensored()) {
    DegenerateGaussianScore scorer(ds, penalty);
    GrowShrink gs(&scorer);
    gs.setVerbose(verbose);
	
    mb = gs.search(targetNode, regressors, &score);
    // } else {
    //     RegressionBicScore scorer(ds, penalty);
    // 	GrowShrink gs(&scorer);
    // 	gs.setVerbose(verbose);

    // 	mb = gs.search(targetNode, regressors, &score);
    // }
    
    RcppThread::checkUserInterrupt();

    Rcpp::StringVector _mb;

    for (Node n : mb) {
	_mb.push_back(n.getName());
    }

    _mb.attr("Score") = Rcpp::wrap(score);
    
    return _mb;
}

