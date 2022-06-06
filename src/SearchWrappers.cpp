#include "MGM.hpp"
#include "PcStable.hpp"
#include "CpcStable.hpp"
#include "PcMax.hpp"
#include "Pc50.hpp"
#include "Fci.hpp"
#include "Cfci.hpp"
#include "FciMax.hpp"
#include "Fci50.hpp"
#include "STEPS.hpp"
// #include "STARS.hpp"
#include "Bootstrap.hpp"
// #include "Tests.hpp"
// #include "IndTestMulti.hpp"


//' Calculate the MGM graph on a dataset
//'
//' @param df The dataframe
//' @param lambda A vector of three lambda values - the first for continuous-continuous interaction, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::mgm(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List mgm(
    const Rcpp::DataFrame &df, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
    const int maxDiscrete = 5,
    const bool verbose = false
) {
    DataSet ds = DataSet(df, maxDiscrete);

    bool v = verbose; // Rcpp::is_true(Rcpp::all(verbose));

    std::vector<double> l(lambda.begin(), lambda.end());

    int lamLength = 3;

    // if (ds.isMixed()) {
    // 	// if (ds.isCensored()) {
    // 	//     lamLength = 5;
    // 	// } else {
    // 	//     lamLength = 3;
    // 	// }
    // 	lamLength = 3;
    // } // else {
    // // 	throw std::runtime_error("MGM is not implemented for purely continuous or purely discrete datasets.");
    // // }

    if (l.size() == 1) {
	for (int i = 1; i < lamLength; i++) {
	    l.push_back(l[0]);
	}
    } else if (l.size() != lamLength) {
	throw std::runtime_error("The regularization parameter lambda should be either a vector of length " + std::to_string(lamLength) + " or a single value for this dataset.");
    }

    MGM mgm(ds, l);
    mgm.setVerbose(v);
    // mgm.calcLambdaMax();
    EdgeListGraph mgmGraph = mgm.search();
    
    mgmGraph.setHyperParam("lambda", Rcpp::NumericVector(l.begin(), l.end()));

    RcppThread::checkUserInterrupt();

    auto elapsedTime = mgm.getElapsedTime();

    if (v) {
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


//' Calculate the solution path for an MGM graph on a dataset
//'
//' @param df The dataframe
//' @param lambdas A range of lambda values used to calculate a solution path for MGM. If NULL, lambdas is set to nLambda logarithmically spaced values from 10*sqrt(log10(p)/n) to sqrt(log10(p)/n). Defaults to NULL.
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::mgm(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List mgmPath(
    const Rcpp::DataFrame &df,
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 30,
    const int maxDiscrete = 5,
    const bool verbose = false
) {
    DataSet ds = DataSet(df, maxDiscrete);

    bool v = verbose; // Rcpp::is_true(Rcpp::all(verbose));

    arma::vec _lambda;
    std::vector<double> l; // = {0.2, 0.2, 0.2};

    int n = ds.getNumRows();

    MGM mgm(ds);
    mgm.setVerbose(v);

    double logLambdaMax = std::log10(mgm.calcLambdaMax());

    logLambdaMax = std::min(logLambdaMax,
    			    std::log10(10 * std::sqrt(std::log10(ds.getNumColumns()) /
						      ((double) ds.getNumRows()))));

    if (lambdas.isNotNull()) {
        Rcpp::NumericVector _lambda(lambdas); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
	arma::vec _lambda = arma::logspace(logLambdaMax-1, logLambdaMax, nLambda); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
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
        Rcpp::Rcout << "MGM Elapsed time =  " << elapsedTime / 1000.0 << " s" << std::endl;
    }

    Rcpp::List graphList;
    
    for (int i = 0; i < nLambda; i++) {
        graphList.push_back(mgmGraphs[i].toList());
    }
    
    Rcpp::List result = Rcpp::List::create(Rcpp::_["graphs"]=graphList,
					   Rcpp::_["lambdas"]=arma::reverse(_lambda),
					   Rcpp::_["loglik"] = loglik,
					   Rcpp::_["AIC"] = 2*nParams - 2*loglik,
					   Rcpp::_["BIC"] = std::log(n)*nParams - 2*loglik);

    

    return result;
}



//' Calculates the optimal lambda values for the MGM algorithm using StEPS and run the algorithm using those values. Optimal values are printed
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param lambdas A range of lambda values assessed for stability by the StEPS algorithm. If NULL, lambdas is set to nLambda logarithmically spaced values from 10*sqrt(log10(p)/n) to sqrt(log10(p)/n). Defaults to NULL.
//' @param g The gamma parameter for STEPS. Defaults to 0.05
//' @param numSub The number of subsets to split the data into. Defaults to 20
//' @param subSize The size of the subsamples used for STEPS. If the value is -1, the size of the subsamples is set to floor(10*sqrt(n)). If the value is in the range (0,1), the size of the subsamples is set to floor(subSize * n). Otherwise, if subSize is in the range [1,n), the size of the subsamples is set to subSize. Defaults to -1.
//' @param leaveOneOut If TRUE, performs leave-one-out subsampling. Defaults to FALSE.
//' @param computeStabs If TRUE, stability values are calculated. Defaults to FALSE.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated MGM graph
//' @export
//' @examples
//' data("data.n100.p25")
//' g <- rCausalMGM::steps(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List steps(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::NumericVector> lambdas = R_NilValue,
    const int nLambda = 20,
    const double g = 0.05,
    const int numSub = 20,
    const int subSize = -1,
    const bool leaveOneOut = false,
    const bool computeStabs = false,
    const int threads = -1,
    const bool verbose = false
) {

    std::vector<double> l;
    
    DataSet ds(df, maxDiscrete);

    MGM mgm(ds);

    double logLambdaMax = std::log10(mgm.calcLambdaMax());

    logLambdaMax = std::min(logLambdaMax,
    			    std::log10(10 * std::sqrt(std::log10(ds.getNumColumns()) /
						      ((double) ds.getNumRows()))));

    if (lambdas.isNotNull()) {
        Rcpp::NumericVector _lambda(lambdas); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
	arma::vec _lambda = arma::logspace(logLambdaMax-1, logLambdaMax, nLambda); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    }

    STEPS steps;
    if (subSize < 0)
	steps = STEPS(ds, l, g, numSub, leaveOneOut);
    else
	steps = STEPS(ds, l, g, numSub, subSize, leaveOneOut);
      
    if (threads > 0) steps.setThreads(threads);
    steps.setComputeStabs(computeStabs);
    steps.setVerbose(verbose);

    Rcpp::List result = steps.runStepsPath().toList();

    if (computeStabs) {
        result["stabilities"] = steps.getStabs();
        std::vector<std::string> names = ds.getVariableNames();
        Rcpp::rownames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
        Rcpp::colnames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    } 

    // ds.deleteVariables();

    return result;

}

//' Runs the causal algorithm PC-Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::pcStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List pcStable(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue, 
    const double alpha = 0.1, 
    const int threads = -1,
    const bool fdr = true,
    const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);

    // bool v = Rcpp::is_true(Rcpp::all(verbose));
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));

    IndTestMulti itm(ds, alpha);

    PcStable pcs((IndependenceTest*) &itm);
    if (threads > 0) pcs.setThreads(threads);
    pcs.setVerbose(verbose);
    pcs.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        pcs.setInitialGraph(&ig);
    }

    Rcpp::List result = pcs.search().toList();

    // ds.deleteVariables();

    return result;
}

//' Runs the causal algorithm CPC-Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::cpcStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List cpcStable(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue, 
    const double alpha = 0.1, 
    const int threads = -1,
    const bool fdr = true,
    const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);

    // bool v = Rcpp::is_true(Rcpp::all(verbose));
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));

    IndTestMulti itm(ds, alpha);

    CpcStable cpc((IndependenceTest*) &itm);
    if (threads > 0) cpc.setThreads(threads);
    cpc.setVerbose(verbose);
    cpc.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        cpc.setInitialGraph(&ig);
    }
    Rcpp::List result = cpc.search().toList();

    // ds.deleteVariables();

    return result;
}

//' Runs the causal algorithm PC-Max on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::pcMax(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List pcMax(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue, 
    const double alpha = 0.1, 
    const int threads = -1,
    const bool fdr = true,
    const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);

    IndTestMulti itm(ds, alpha);

    PcMax pcm((IndependenceTest*) &itm);
    if (threads > 0) pcm.setThreads(threads);
    pcm.setVerbose(verbose);
    pcm.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        pcm.setInitialGraph(&ig);
    }
    Rcpp::List result = pcm.search().toList();

    return result;
}


//' Runs the causal algorithm PC50 on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::pc50(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List pc50(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue, 
    const double alpha = 0.1, 
    const int threads = -1,
    const bool fdr = true,
    const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);

    // bool v = Rcpp::is_true(Rcpp::all(verbose));
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));

    IndTestMulti itm(ds, alpha);

    Pc50 pc50((IndependenceTest*) &itm);
    if (threads > 0) pc50.setThreads(threads);
    pc50.setVerbose(verbose);
    pc50.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        pc50.setInitialGraph(&ig);
    }
    Rcpp::List result = pc50.search().toList();

    // ds.deleteVariables();

    return result;
}



//' Runs the causal algorithm FCI-Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::fciStable(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List fciStable(
        const Rcpp::DataFrame &df,
        const int maxDiscrete = 5,
        Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
        const double alpha = 0.1,
        const int threads = -1,
	const bool fdr = true,
        const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);
    
    // bool v = Rcpp::is_true(Rcpp::all(verbose));
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    
    IndTestMulti itm(ds, alpha);
    
    Fci fci((IndependenceTest*) &itm);
    if (threads > 0) fci.setThreads(threads);
    fci.setVerbose(verbose);
    fci.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        fci.setInitialGraph(&ig);
    }
    
    Rcpp::List result = fci.search().toList();
    
    // ds.deleteVariables();
    
    return result;
}


//' Runs the causal algorithm CFCI-Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::cfci(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List cfci(
        const Rcpp::DataFrame &df,
        const int maxDiscrete = 5,
        Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
        const double alpha = 0.1,
        const int threads = -1,
	const bool fdr = true,
        const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);
    
    // bool v = Rcpp::is_true(Rcpp::all(verbose));
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    
    IndTestMulti itm(ds, alpha);
    
    Cfci cfci((IndependenceTest*) &itm);
    if (threads > 0) cfci.setThreads(threads);
    cfci.setVerbose(verbose);
    cfci.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        cfci.setInitialGraph(&ig);
    }
    
    Rcpp::List result = cfci.search().toList();
    
    // ds.deleteVariables();
    
    return result;
}


//' Runs the causal algorithm FCI-Max on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::fciMax(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List fciMax(
        const Rcpp::DataFrame &df,
        const int maxDiscrete = 5,
        Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
        const double alpha = 0.1,
        const int threads = -1,
	const bool fdr = true,
        const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);
    
    // bool v = Rcpp::is_true(Rcpp::all(verbose));
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    
    IndTestMulti itm(ds, alpha);
    
    FciMax fcimax((IndependenceTest*) &itm);
    if (threads > 0) fcimax.setThreads(threads);
    fcimax.setVerbose(verbose);
    fcimax.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        fcimax.setInitialGraph(&ig);
    }

    Rcpp::List result;
    try {
	result = fcimax.search().toList();
    } catch(std::exception& e) {
	Rcpp::Rcout << e.what() << std::endl;
    }
        
    return result;
}


//' Runs the causal algorithm FCI50 Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph An initial undirected graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param fdr Whether or not to run with FDR correction for the adjacencies.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph
//' @export
//' @examples
//' data("data.n100.p25")
//' ig <- rCausalMGM::mgm(data.n100.p25)
//' g <- rCausalMGM::fci50(data.n100.p25, initialGraph = ig)
// [[Rcpp::export]]
Rcpp::List fci50(
        const Rcpp::DataFrame &df,
        const int maxDiscrete = 5,
        Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
        const double alpha = 0.1,
        const int threads = -1,
	const bool fdr = true,
        const bool verbose = false
) {
    DataSet ds(df, maxDiscrete);
    
    // bool v = Rcpp::is_true(Rcpp::all(verbose));
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    
    IndTestMulti itm(ds, alpha);
    
    Fci50 fci50((IndependenceTest*) &itm);
    if (threads > 0) fci50.setThreads(threads);
    fci50.setVerbose(verbose);
    fci50.setFDR(fdr);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        fci50.setInitialGraph(&ig);
    }
    
    Rcpp::List result = fci50.search().toList();
    
    // ds.deleteVariables();
    
    return result;
}


// // no export // [[Rcpp::export]]
// Rcpp::List stars(
//     const Rcpp::DataFrame& df,
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
    
//     DataSet ds(df, maxDiscrete);

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
//' @param df The dataframe
//' @param algorithm string indicating the name of the causal discovery algorithm to bootstrap. Causal discovery algorithms can be run alone or with mgm to learn an initial graph. Options include "mgm", "pc", "cpc", "pcm", "pc50", "fci", "cfci", "fcim", "mgm-pc", "mgm-cpc", "mgm-pcm", "mgm-pc50", "mgm-fci", "mgm-cfci", "mgm-fcim", "mgm-fci50." The default value is set to "mgm-pc50."
//' @param ensemble Method for construncting an ensemble graph from bootstrapped graphs. Options include "highest", which orients edges according to the orientation with the highest bootstrap probability, or "majority", which only orients edges if they have an orientation with a bootstrap probability > 0.5. Otherwise, the adjacency is included but the edge is left unoreineted. The default value is "highest."
//' @param lambda A vector of three lambda values - the first for continuous-continuous interaction, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2). If a single value is provided, all three values in the vector will be set to that value.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param numBoots The number of bootstrap samples to run. Defaults to 20.
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
//' @param verbose Whether or not to output additional information. Defaults to FALSE.
//' @return The calculated search graph with a table of edge stabilities
//' @export
//' @examples
//' data("data.n100.p25")
//' g.boot <- rCausalMGM::bootstrap(data.n100.p25)
// [[Rcpp::export]]
Rcpp::List bootstrap(
    const Rcpp::DataFrame& df,
    Rcpp::StringVector algorithm = Rcpp::CharacterVector::create("mgm-pc50", "mgm", "pc", "cpc", "pcm", "pc50", "fci", "cfci", "fcim", "mgm-pc", "mgm-cpc", "mgm-pcm", "mgm-fci", "mgm-cfci", "mgm-fcim", "mgm-fci50"),
    Rcpp::StringVector ensemble = Rcpp::CharacterVector::create("highest", "majority"),
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
    const double alpha = 0.05,
    const int numBoots = 20,
    const int maxDiscrete = 5,
    const int threads = -1,
    const bool verbose = false
    ) {

    // Rcpp::Rcout << "running stars...\n";

    std::string alg, _method, _ensemble;
    // bool _fdr = Rcpp::is_true(Rcpp::all(fdr));
    // bool v = Rcpp::is_true(Rcpp::all(verbose));

    _method = algorithm[0];

    _ensemble = ensemble[0];

    // Rcpp::Rcout << _method << std::endl;
    std::transform(_method.begin(), _method.end(), _method.begin(),
		   [](unsigned char c){ return std::tolower(c); });
    // Rcpp::Rcout << _method << std::endl;
    _method.erase(std::remove(_method.begin(), _method.end(), '-'), _method.end());
    // Rcpp::Rcout << _method << std::endl;

    std::transform(_ensemble.begin(), _ensemble.end(), _ensemble.begin(),
		   [](unsigned char c){ return std::tolower(c); });

    if (_method == "mgm") {
	alg = "mgm";
    } else if (_method == "pc" || _method == "pcs" || _method == "pcstable") {
	alg = "pc";
    } else if (_method == "cpc" || _method == "cpcstable") {
	alg = "cpc";
    } else if (_method == "pcm" || _method == "pcmax") {
	alg = "pcm";
    } else if (_method == "pc50") {
	alg = "pc50";
    }else if (_method == "fci" || _method == "fcistable") {
	alg = "fci";
    } else if (_method == "cfci" || _method == "cfcistable") {
	alg = "cfci";
    } else if (_method == "fcim" || _method == "fcimax") {
	alg = "fcim";
    } else if (_method == "fci50") {
	alg = "fci50";
    } else if (_method == "mgmpc" || _method == "mgmpcs" || _method == "mgmpcstable") {
	alg = "mgmpc";
    } else if (_method == "mgmcpc" || _method == "mgmcpcstable") {
	alg = "mgmcpc";
    } else if (_method == "mgmpcm" || _method == "mgmpcmax") {
	alg = "mgmpcm";
    } else if (_method == "mgmpc50") {
	alg = "mgmpc50";
    }else if (_method == "mgmfci" || _method == "mgmfcistable") {
	alg = "mgmfci";
    } else if (_method == "mgmcfci" || _method == "mgmcfcistable") {
	alg = "mgmcfci";
    } else if (_method == "mgmfcim" || _method == "mgmfcimax") {
	alg = "mgmfcim";
    } else if (_method == "mgmfci50") {
	alg = "mgmfci50";
    } else {
	throw std::invalid_argument("Invalid algorithm: " + _method
				    + "\n   Algorithm must be in the list: "
				    + "{ mgm, pc, cpc, pcm, pc50, fci, cfci, fcim, fci50, "
				    + "mgm-pc, mgm-cpc, mgm-pcm, mgm-pc50, mgm-fci, "
				    + "mgm-cfci, mgm-fcim, mgm-fci50 }");
    }

    std::vector<double> l(lambda.begin(), lambda.end());

    int lamLength = 3;

    // if (ds.isMixed()) {
    // 	// if (ds.isCensored()) {
    // 	//     lamLength = 5;
    // 	// } else {
    // 	//     lamLength = 3;
    // 	// }
    // 	lamLength = 3;
    // } // else {
    // // 	throw std::runtime_error("MGM is not implemented for purely continuous or purely discrete datasets.");
    // // }

    if (l.size() == 1) {
	for (int i = 1; i < lamLength; i++) {
	    l.push_back(l[0]);
	}
    } else if (l.size() != lamLength) {
	throw std::runtime_error("The regularization parameter lambda should be either a vector of length " + std::to_string(lamLength) + " or a single value for this dataset.");
    }
    
    DataSet ds(df, maxDiscrete);

    Bootstrap boot(ds, alg, _ensemble, numBoots);
    if (threads > 0) boot.setThreads(threads);
    boot.setVerbose(verbose);
    boot.setAlpha(alpha);
    boot.setLambda(l);
    // boot.setFdr(fdr);

    Rcpp::List result = boot.runBootstrap().toList();

    result["stabilities"] = boot.getStabs();

    return result;

}
