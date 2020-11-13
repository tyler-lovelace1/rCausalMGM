// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// DataSetTest
void DataSetTest(const Rcpp::DataFrame& df, const int maxDiscrete);
RcppExport SEXP _rCausalMGM_DataSetTest(SEXP dfSEXP, SEXP maxDiscreteSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    DataSetTest(df, maxDiscrete);
    return R_NilValue;
END_RCPP
}
// saveGraph
void saveGraph(const Rcpp::List& list, const std::string& filename);
RcppExport SEXP _rCausalMGM_saveGraph(SEXP listSEXP, SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type list(listSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    saveGraph(list, filename);
    return R_NilValue;
END_RCPP
}
// loadGraph
Rcpp::List loadGraph(const std::string& filename);
RcppExport SEXP _rCausalMGM_loadGraph(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(loadGraph(filename));
    return rcpp_result_gen;
END_RCPP
}
// printGraph
void printGraph(const Rcpp::List& graph, const Rcpp::DataFrame& df);
RcppExport SEXP _rCausalMGM_printGraph(SEXP graphSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    printGraph(graph, df);
    return R_NilValue;
END_RCPP
}
// indTestMultiTest
void indTestMultiTest(const Rcpp::DataFrame& df);
RcppExport SEXP _rCausalMGM_indTestMultiTest(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    indTestMultiTest(df);
    return R_NilValue;
END_RCPP
}
// LinearRegressionTest
void LinearRegressionTest(const Rcpp::DataFrame& df);
RcppExport SEXP _rCausalMGM_LinearRegressionTest(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    LinearRegressionTest(df);
    return R_NilValue;
END_RCPP
}
// LogisticRegressionTest
void LogisticRegressionTest(const Rcpp::DataFrame& df);
RcppExport SEXP _rCausalMGM_LogisticRegressionTest(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    LogisticRegressionTest(df);
    return R_NilValue;
END_RCPP
}
// mgm
Rcpp::List mgm(const Rcpp::DataFrame& df, Rcpp::NumericVector lambda, const int maxDiscrete, Rcpp::LogicalVector verbose);
RcppExport SEXP _rCausalMGM_mgm(SEXP dfSEXP, SEXP lambdaSEXP, SEXP maxDiscreteSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mgm(df, lambda, maxDiscrete, verbose));
    return rcpp_result_gen;
END_RCPP
}
// steps
Rcpp::List steps(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::NumericVector lambda, const double g, const int numSub, Rcpp::LogicalVector leaveOneOut, Rcpp::LogicalVector computeStabs, Rcpp::LogicalVector verbose);
RcppExport SEXP _rCausalMGM_steps(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP lambdaSEXP, SEXP gSEXP, SEXP numSubSEXP, SEXP leaveOneOutSEXP, SEXP computeStabsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type leaveOneOut(leaveOneOutSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type computeStabs(computeStabsSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(steps(df, maxDiscrete, lambda, g, numSub, leaveOneOut, computeStabs, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pcStable
Rcpp::List pcStable(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::NumericVector lambda, const double alpha, Rcpp::LogicalVector verbose);
RcppExport SEXP _rCausalMGM_pcStable(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcStable(df, maxDiscrete, initialGraph, lambda, alpha, verbose));
    return rcpp_result_gen;
END_RCPP
}
// cpcStable
Rcpp::List cpcStable(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::NumericVector lambda, const double alpha, Rcpp::LogicalVector verbose);
RcppExport SEXP _rCausalMGM_cpcStable(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cpcStable(df, maxDiscrete, initialGraph, lambda, alpha, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pcMax
Rcpp::List pcMax(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::NumericVector lambda, const double alpha, Rcpp::LogicalVector verbose);
RcppExport SEXP _rCausalMGM_pcMax(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcMax(df, maxDiscrete, initialGraph, lambda, alpha, verbose));
    return rcpp_result_gen;
END_RCPP
}
// runTests
void runTests(const Rcpp::DataFrame& df, const int maxDiscrete);
RcppExport SEXP _rCausalMGM_runTests(SEXP dfSEXP, SEXP maxDiscreteSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    runTests(df, maxDiscrete);
    return R_NilValue;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _rCausalMGM_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _rCausalMGM_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _rCausalMGM_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _rCausalMGM_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rCausalMGM_DataSetTest", (DL_FUNC) &_rCausalMGM_DataSetTest, 2},
    {"_rCausalMGM_saveGraph", (DL_FUNC) &_rCausalMGM_saveGraph, 2},
    {"_rCausalMGM_loadGraph", (DL_FUNC) &_rCausalMGM_loadGraph, 1},
    {"_rCausalMGM_printGraph", (DL_FUNC) &_rCausalMGM_printGraph, 2},
    {"_rCausalMGM_indTestMultiTest", (DL_FUNC) &_rCausalMGM_indTestMultiTest, 1},
    {"_rCausalMGM_LinearRegressionTest", (DL_FUNC) &_rCausalMGM_LinearRegressionTest, 1},
    {"_rCausalMGM_LogisticRegressionTest", (DL_FUNC) &_rCausalMGM_LogisticRegressionTest, 1},
    {"_rCausalMGM_mgm", (DL_FUNC) &_rCausalMGM_mgm, 4},
    {"_rCausalMGM_steps", (DL_FUNC) &_rCausalMGM_steps, 8},
    {"_rCausalMGM_pcStable", (DL_FUNC) &_rCausalMGM_pcStable, 6},
    {"_rCausalMGM_cpcStable", (DL_FUNC) &_rCausalMGM_cpcStable, 6},
    {"_rCausalMGM_pcMax", (DL_FUNC) &_rCausalMGM_pcMax, 6},
    {"_rCausalMGM_runTests", (DL_FUNC) &_rCausalMGM_runTests, 2},
    {"_rCausalMGM_rcpparma_hello_world", (DL_FUNC) &_rCausalMGM_rcpparma_hello_world, 0},
    {"_rCausalMGM_rcpparma_outerproduct", (DL_FUNC) &_rCausalMGM_rcpparma_outerproduct, 1},
    {"_rCausalMGM_rcpparma_innerproduct", (DL_FUNC) &_rCausalMGM_rcpparma_innerproduct, 1},
    {"_rCausalMGM_rcpparma_bothproducts", (DL_FUNC) &_rCausalMGM_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rCausalMGM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
