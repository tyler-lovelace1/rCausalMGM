#include "MGM.hpp"
#include "PcStable.hpp"
#include "CpcStable.hpp"
#include "PcMax.hpp"
#include "Pc50.hpp"
#include "Fci.hpp"
#include "Cfci.hpp"
#include "FciMax.hpp"
#include "STEPS.hpp"
#include "STARS.hpp"
#include "Tests.hpp"
#include "IndTestMultiCC.hpp"


//' Calculate the MGM graph on a dataset
//'
//' @param df The dataframe
//' @param lambda A vector of three lambda values - the first for continuous-continuous interaction, the second for continuous-discrete, and the third for discrete-discrete. Defaults to c(0.2, 0.2, 0.2)
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
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2), 
    const int maxDiscrete = 5,
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds = DataSet(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    std::vector<double> l(lambda.begin(), lambda.end());

    int lamLength;

    if (ds.isMixed()) {
	if (ds.isCensored()) {
	    lamLength = 5;
	} else {
	    lamLength = 3;
	}
    } else {
	throw std::runtime_error("MGM is not implemented for purely continuous or purely discrete datasets.");
    }

    if (l.size() == 1) {
	for (int i = 1; i < lamLength; i++) {
	    l.push_back(l[0]);
	}
    } else if (l.size() != lamLength) {
	throw std::runtime_error("The length of the vector of regularization parameters should be either " + std::to_string(lamLength) + " for this dataset.");
    }

    MGM mgm(ds, l);
    mgm.setVerbose(v);
    EdgeListGraph mgmGraph = mgm.search();

    if (v) {
        Rcpp::Rcout.precision(2);
        Rcpp::Rcout << "MGM Elapsed time =  " << (mgm.getElapsedTime() / 1000.0) << " s" << std::endl;
    }

    Rcpp::List result = mgmGraph.toList();

    // ds.deleteVariables();

    return result;
}
 
//' Calculate the optimal lambda values for the MGM algorithm and run the algorithm using those values. Optimal values are printed
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param lambda A vector of the lambda values to test. Defaults to a logspaced vector with 20 values ranging from 0.9 to 0.09 if n < p, or from 0.9 to 0.009 if n > p.
//' @param g The gamma parameter for STEPS. Defaults to 0.05
//' @param numSub The number of subsets to split the data into. Defaults to 20
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
    const Rcpp::DataFrame& df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::NumericVector> lambda = R_NilValue, 
    const double g = 0.05,
    const int numSub = 20,
    Rcpp::LogicalVector leaveOneOut = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector computeStabs = Rcpp::LogicalVector::create(FALSE),
    const int threads = -1,
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {

    Rcpp::Rcout << "running steps...\n";

    std::vector<double> l;
    bool loo = Rcpp::is_true(Rcpp::all(leaveOneOut));
    bool cs = Rcpp::is_true(Rcpp::all(computeStabs));
    bool v = Rcpp::is_true(Rcpp::all(verbose));

    Rcpp::Rcout << "creating dataset object...\n";
    DataSet ds(df, maxDiscrete);
    Rcpp::Rcout << "dataset object created\n";

    if (lambda.isNotNull()) {
        Rcpp::NumericVector _lambda(lambda); 
	l = std::vector<double>(_lambda.begin(), _lambda.end());
    } else {
        if (ds.getNumRows() > ds.getNumColumns()) {
            arma::vec _lambda = arma::logspace(std::log10(0.9)-2, std::log10(0.9), 20); 
            l = std::vector<double>(_lambda.begin(), _lambda.end());
        } else {
            arma::vec _lambda = arma::logspace(std::log10(0.9)-1, std::log10(0.9), 20); 
            l = std::vector<double>(_lambda.begin(), _lambda.end());
        }
    }

    Rcpp::Rcout << "lambda vector filled\n";
    STEPS steps(ds, l, g, numSub, loo);
    Rcpp::Rcout << "steps object created\n";
    if (threads > 0) steps.setThreads(threads);
    steps.setComputeStabs(cs);
    steps.setVerbose(v);

    Rcpp::List result = steps.runStepsPar().toList();

    if (cs) {
        result["stabilities"] = steps.getStabs();
        std::vector<std::string> names = ds.getVariableNames();
        Rcpp::rownames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
        Rcpp::colnames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    } 

    // ds.deleteVariables();

    return result;

}

//' Runs the causal algorithm PC Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
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
    const double alpha = 0.05, 
    const int threads = -1,
    Rcpp::LogicalVector censoredCauses = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    IndTestMultiCC itm(ds, alpha);

    PcStable pcs((IndependenceTest*) &itm);
    if (threads > 0) pcs.setThreads(threads);
    pcs.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
	ig.setCensoredCauses(Rcpp::is_true(Rcpp::all(censoredCauses)));
        pcs.setInitialGraph(&ig);
    }

    Rcpp::List result = pcs.search().toList();

    // ds.deleteVariables();

    return result;
}

//' Runs the causal algorithm CPC Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
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
    const double alpha = 0.05, 
    const int threads = -1,
    Rcpp::LogicalVector censoredCauses = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    IndTestMultiCC itm(ds, alpha);

    CpcStable cpc((IndependenceTest*) &itm);
    if (threads > 0) cpc.setThreads(threads);
    cpc.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
	ig.setCensoredCauses(Rcpp::is_true(Rcpp::all(censoredCauses)));
        cpc.setInitialGraph(&ig);
    }
    Rcpp::List result = cpc.search().toList();

    // ds.deleteVariables();

    return result;
}

//' Runs the causal algorithm PC Max on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
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
    const double alpha = 0.05, 
    const int threads = -1,
    Rcpp::LogicalVector censoredCauses = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    IndTestMultiCC itm(ds, alpha);

    PcMax pcm((IndependenceTest*) &itm);
    if (threads > 0) pcm.setThreads(threads);
    pcm.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
	ig.setCensoredCauses(Rcpp::is_true(Rcpp::all(censoredCauses)));
        pcm.setInitialGraph(&ig);
    }
    Rcpp::List result = pcm.search().toList();

    // ds.deleteVariables();

    return result;
}


//' Runs the causal algorithm PC50 on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
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
    const double alpha = 0.05, 
    const int threads = -1,
    Rcpp::LogicalVector censoredCauses = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    IndTestMultiCC itm(ds, alpha);

    Pc50 pc50((IndependenceTest*) &itm);
    if (threads > 0) pc50.setThreads(threads);
    pc50.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
	ig.setCensoredCauses(Rcpp::is_true(Rcpp::all(censoredCauses)));
        pc50.setInitialGraph(&ig);
    }
    Rcpp::List result = pc50.search().toList();

    // ds.deleteVariables();

    return result;
}



//' Runs the causal algorithm FCI Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
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
    const double alpha = 0.05,
    const int threads = -1,
    Rcpp::LogicalVector censoredCauses = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds(df, maxDiscrete);
    
    bool v = Rcpp::is_true(Rcpp::all(verbose));
    
    IndTestMultiCC itm(ds, alpha);
    
    Fci fci((IndependenceTest*) &itm);
    if (threads > 0) fci.setThreads(threads);
    fci.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
	ig.setCensoredCauses(Rcpp::is_true(Rcpp::all(censoredCauses)));
        fci.setInitialGraph(&ig);
    }
    
    Rcpp::List result = fci.search().toList();
    
    // ds.deleteVariables();
    
    return result;
}


//' Runs the causal algorithm CFCI Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
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
    const double alpha = 0.05,
    const int threads = -1,
    Rcpp::LogicalVector censoredCauses = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds(df, maxDiscrete);
    
    bool v = Rcpp::is_true(Rcpp::all(verbose));
    
    IndTestMultiCC itm(ds, alpha);
    
    Cfci cfci((IndependenceTest*) &itm);
    if (threads > 0) cfci.setThreads(threads);
    cfci.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
	ig.setCensoredCauses(Rcpp::is_true(Rcpp::all(censoredCauses)));
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
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
//' @param threads The number of consumer threads to create during multi-threaded steps. If -1, defaults to number of availible processors.
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
    const double alpha = 0.05,
    const int threads = -1,
    Rcpp::LogicalVector censoredCauses = Rcpp::LogicalVector::create(FALSE),
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {
    DataSet ds(df, maxDiscrete);
    
    bool v = Rcpp::is_true(Rcpp::all(verbose));
    
    IndTestMultiCC itm(ds, alpha);
    
    FciMax fcimax((IndependenceTest*) &itm);
    if (threads > 0) fcimax.setThreads(threads);
    fcimax.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
	ig.setCensoredCauses(Rcpp::is_true(Rcpp::all(censoredCauses)));
        fcimax.setInitialGraph(&ig);
    }

    Rcpp::List result;
    try {
	result = fcimax.search().toList();
    } catch(std::exception& e) {
	Rcpp::Rcout << e.what() << std::endl;
    }
    
    // ds.deleteVariables();
    
    return result;
}

// Rcpp::LogicalVector computeStabs = Rcpp::LogicalVector::create(FALSE),
// [[Rcpp::export]]
Rcpp::List stars(
    const Rcpp::DataFrame& df,
    const std::string method, 
    Rcpp::Nullable<Rcpp::NumericVector> params = R_NilValue,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue,
    const int maxDiscrete = 5,
    const double g = 0.05,
    const int numSub = 20,
    Rcpp::LogicalVector leaveOneOut = Rcpp::LogicalVector::create(FALSE),
    const int threads = -1,
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(FALSE)
    ) {

    Rcpp::Rcout << "running stars...\n";

    std::string alg, _method;
    std::vector<double> par;
    bool loo = Rcpp::is_true(Rcpp::all(leaveOneOut));
    // bool cs = Rcpp::is_true(Rcpp::all(computeStabs));
    bool v = Rcpp::is_true(Rcpp::all(verbose));

    _method = method;

    Rcpp::Rcout << _method << std::endl;
    std::transform(_method.begin(), _method.end(), _method.begin(),
		   [](unsigned char c){ return std::tolower(c); });
    Rcpp::Rcout << _method << std::endl;

    if (_method == "mgm") {
	alg = "mgm";
    } else if (_method == "pc" || _method == "pcs" || _method == "pcstable") {
	alg = "pc";
    } else if (_method == "cpc" || _method == "cpcstable") {
	alg = "cpc";
    } else if (_method == "pcm" || _method == "pcmax") {
	alg = "pcm";
    } else if (_method == "fci" || _method == "fcistable") {
	alg = "fci";
    } else if (_method == "cfci" || _method == "cfcistable") {
	alg = "cfci";
    } else if (_method == "fcim" || _method == "fcimax") {
	alg = "fcim";
    } else {
	throw std::invalid_argument("Invalid algorithm: " + _method
				    + "\n   Algorithm must be in the list: "
				    + "{ mgm, pc, cpc, pcm, fci, cfci, fcim }");
    }
    
    DataSet ds(df, maxDiscrete);

    if (params.isNotNull()) {
        Rcpp::NumericVector _params(params); 
	par = std::vector<double>(_params.begin(), _params.end());
    } else {
	if (alg == "mgm") {
	    if (ds.getNumRows() > ds.getNumColumns()) {
		arma::vec _params = arma::logspace(std::log10(0.9)-2, std::log10(0.9), 20); 
		par = std::vector<double>(_params.begin(), _params.end());
	    } else {
		arma::vec _params = arma::logspace(std::log10(0.9)-1, std::log10(0.9), 20); 
		par = std::vector<double>(_params.begin(), _params.end());
	    }
	} else {
	    par = { 0.001, 0.005, 0.01, 0.05, 0.1 };
	}
    }

    Rcpp::Rcout << "params vector filled\n";
    STARS stars(ds, alg, par, g, numSub, loo);
    Rcpp::Rcout << "stars object created\n";
    if (threads > 0) stars.setThreads(threads);
    // stars.setComputeStabs(cs);
    stars.setVerbose(v);

    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        stars.setInitialGraph(&ig);
    }

    Rcpp::List result = stars.runStarsPar().toList();

    // if (cs) {
    //     result["stabilities"] = steps.getStabs();
    //     std::vector<std::string> names = ds.getVariableNames();
    //     Rcpp::rownames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    //     Rcpp::colnames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    // } 

    // ds.deleteVariables();

    return result;

}
