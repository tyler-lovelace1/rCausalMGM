#include "MGM.hpp"
#include "PcStable.hpp"
#include "CpcStable.hpp"
#include "PcMax.hpp"
#include "STEPS.hpp"
#include "Tests.hpp"
#include "IndTestMulti.hpp"


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
//' g <- rCausalMGM::mgm("data.n100.p25")
// [[Rcpp::export]]
Rcpp::List mgm(
    const Rcpp::DataFrame &df, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2), 
    const int maxDiscrete = 5,
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    std::vector<double> l(lambda.begin(), lambda.end());

    MGM mgm(ds, l);
    mgm.setVerbose(v);
    EdgeListGraph mgmGraph = mgm.search();

    if (v) {
        Rcpp::Rcout.precision(2);
        Rcpp::Rcout << "MGM Elapsed time =  " << (mgm.getElapsedTime() / 1000.0) << " s" << std::endl;
    }

    Rcpp::List result = mgmGraph.toList();

    ds.deleteVariables();

    return result;
}

// Rcpp::NumericVector::create(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85),
 
//' Calculate the optimal lambda values for the MGM algorithm and run the algorithm using those values. Optimal values are printed
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param lambda A vector of the lambda values to test. Defaults to a logspaced vector with 20 values ranging from 0.9 to 0.09 if n < p, or from 0.9 to 0.009 if n > p.
//' @param g The gamma parameter for STEPS. Defaults to 0.05
//' @param numSub The number of subsets to split the data into. Defaults to 20
//' @param leaveOneOut If TRUE, performs leave-one-out subsampling. Defaults to FALSE.
//' @param computeStabs If TRUE, stability values are calculated. Defaults to FALSE.
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
    Rcpp::Nullable<Rcpp::NumericVector> lambda = R_NilValue, 
    const double g = 0.05,
    const int numSub = 20,
    Rcpp::LogicalVector leaveOneOut = Rcpp::LogicalVector::create(0), // FALSE
    Rcpp::LogicalVector computeStabs = Rcpp::LogicalVector::create(0), // FALSE
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {

    std::vector<double> l;
    bool loo = Rcpp::is_true(Rcpp::all(leaveOneOut));
    bool cs = Rcpp::is_true(Rcpp::all(computeStabs));
    bool v = Rcpp::is_true(Rcpp::all(verbose));

    DataSet ds(df, maxDiscrete);

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

    STEPS steps(ds, l, g, numSub, loo);
    steps.setComputeStabs(cs);
    steps.setVerbose(v);

    Rcpp::List result = steps.runStepsPar().toList();

    if (cs) {
        result["stabilities"] = steps.getStabs();
        std::vector<std::string> names = ds.getVariableNames();
        Rcpp::rownames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
        Rcpp::colnames(result["stabilities"]) = Rcpp::CharacterVector::import(names.begin(), names.end());
    } 

    ds.deleteVariables();

    return result;

}

//' Runs the causal algorithm PC Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
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
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    IndTestMulti itm(ds, alpha);

    PcStable pcs((IndependenceTest*) &itm);
    pcs.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        pcs.setInitialGraph(&ig);
    }

    Rcpp::List result = pcs.search().toList();

    ds.deleteVariables();

    return result;
}

//' Runs the causal algorithm CPC Stable on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
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
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    IndTestMulti itm(ds, alpha);

    CpcStable cpc((IndependenceTest*) &itm);
    cpc.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        cpc.setInitialGraph(&ig);
    }
    Rcpp::List result = cpc.search().toList();

    ds.deleteVariables();

    return result;
}

//' Runs the causal algorithm PC Max on a dataset
//'
//' @param df The dataframe
//' @param maxDiscrete The maximum number of unique values a variable can have before being considered continuous. Defaults to 5
//' @param initialGraph The MGM graph to use as a starting point. If NULL, a full graph will be used. Defaults to NULL.
//' @param alpha The p value below which results are considered significant. Defaults to 0.05.
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
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    IndTestMulti itm(ds, alpha);

    PcMax pcm((IndependenceTest*) &itm);
    pcm.setVerbose(v);
    EdgeListGraph ig;
    if (!initialGraph.isNull()) {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
        pcm.setInitialGraph(&ig);
    }
    Rcpp::List result = pcm.search().toList();

    ds.deleteVariables();

    return result;
}
