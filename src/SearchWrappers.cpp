#include "MGM.hpp"
#include "PcStable.hpp"
#include "CpcStable.hpp"
#include "PcMax.hpp"
#include "STEPS.hpp"
#include "Tests.hpp"
#include "IndTestMulti.hpp"

// [[Rcpp::export]]
Rcpp::List mgm(
    const Rcpp::DataFrame &df, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2), 
    const int maxDiscrete = 5,
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = verbose[0];

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

// [[Rcpp::export]]
Rcpp::List steps(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85), 
    const double g = 0.05,
    const int numSub = 20,
    Rcpp::LogicalVector leaveOneOut = Rcpp::LogicalVector::create(0), // FALSE
    Rcpp::LogicalVector computeStabs = Rcpp::LogicalVector::create(0), // FALSE
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {

    std::vector<double> l(lambda.begin(), lambda.end());
    bool loo = Rcpp::is_true(Rcpp::all(leaveOneOut));
    bool cs = Rcpp::is_true(Rcpp::all(computeStabs));
    bool v = Rcpp::is_true(Rcpp::all(verbose));

    DataSet ds(df, maxDiscrete);

    STEPS steps(ds, l, g, numSub, loo);
    steps.setComputeStabs(cs);
    steps.setVerbose(v);

    Rcpp::List result = steps.runStepsPar().toList();

    ds.deleteVariables();

    return result;

}

// [[Rcpp::export]]
Rcpp::List pcStable(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
    const double alpha = 0.05, 
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    EdgeListGraph ig;
    if (initialGraph.isNull()) {
        std::vector<double> l(lambda.begin(), lambda.end());
        MGM mgm(ds, l);
        mgm.setVerbose(v);
        ig = mgm.search();
        if (v)  {
            Rcpp::Rcout << "MGM graph: \n" << ig << std::endl;
            Rcpp::Rcout.precision(2);
            Rcpp::Rcout << "MGM Elapsed time =  " << (mgm.getElapsedTime() / 1000.0) << " s" << std::endl;
        }
    } else {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
    }

    IndTestMulti itm(ds, alpha);

    PcStable pcs((IndependenceTest*) &itm);
    pcs.setVerbose(v);
    pcs.setInitialGraph(&ig);
    Rcpp::List result = pcs.search().toList();

    ds.deleteVariables();

    return result;
}

// [[Rcpp::export]]
Rcpp::List cpcStable(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
    const double alpha = 0.05, 
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    EdgeListGraph ig;
    if (initialGraph.isNull()) {
        std::vector<double> l(lambda.begin(), lambda.end());
        MGM mgm(ds, l);
        mgm.setVerbose(v);
        ig = mgm.search();
        if (v)  {
            Rcpp::Rcout << "MGM graph: \n" << ig << std::endl;
            Rcpp::Rcout.precision(2);
            Rcpp::Rcout << "MGM Elapsed time =  " << (mgm.getElapsedTime() / 1000.0) << " s" << std::endl;
        }
    } else {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
    }

    IndTestMulti itm(ds, alpha);

    CpcStable cpc((IndependenceTest*) &itm);
    cpc.setVerbose(v);
    cpc.setInitialGraph(&ig);
    Rcpp::List result = cpc.search().toList();

    ds.deleteVariables();

    return result;
}

// [[Rcpp::export]]
Rcpp::List pcMax(
    const Rcpp::DataFrame &df, 
    const int maxDiscrete = 5,
    Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue, 
    Rcpp::NumericVector lambda = Rcpp::NumericVector::create(0.2, 0.2, 0.2),
    const double alpha = 0.05, 
    Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(0) // FALSE
) {
    DataSet ds(df, maxDiscrete);

    bool v = Rcpp::is_true(Rcpp::all(verbose));

    EdgeListGraph ig;
    if (initialGraph.isNull()) {
        std::vector<double> l(lambda.begin(), lambda.end());
        MGM mgm(ds, l);
        mgm.setVerbose(v);
        ig = mgm.search();
        if (v)  {
            Rcpp::Rcout << "MGM graph: \n" << ig << std::endl;
            Rcpp::Rcout.precision(2);
            Rcpp::Rcout << "MGM Elapsed time =  " << (mgm.getElapsedTime() / 1000.0) << " s" << std::endl;
        }
    } else {
        Rcpp::List _initialGraph(initialGraph);
        ig = EdgeListGraph(_initialGraph, ds);
    }

    IndTestMulti itm(ds, alpha);

    PcMax pcm((IndependenceTest*) &itm);
    pcm.setVerbose(v);
    pcm.setInitialGraph(&ig);
    Rcpp::List result = pcm.search().toList();

    ds.deleteVariables();

    return result;
}

// [[Rcpp::export]]
void runTests(const Rcpp::DataFrame &df, const int maxDiscrete = 5) {

    // Tests::testMGMFunctions(df, maxDiscrete);

    // Tests::testConcurrentQueue();

    // Tests::testPcStable(df, maxDiscrete);

    // Tests::testCpcStable(df, maxDiscrete);

    // Tests::testPcMax(df, maxDiscrete);

    // Tests::testMGMTiming(df, maxDiscrete);

    Tests::testSTEPS(df, maxDiscrete);

    // Tests::testGraphFromFile(df, "data/graph/graph5.txt", maxDiscrete);

}
