// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
// adjMat2Graph
Rcpp::List adjMat2Graph(arma::mat adj, Rcpp::StringVector nodes, bool directed);
RcppExport SEXP _rCausalMGM_adjMat2Graph(SEXP adjSEXP, SEXP nodesSEXP, SEXP directedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type adj(adjSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< bool >::type directed(directedSEXP);
    rcpp_result_gen = Rcpp::wrap(adjMat2Graph(adj, nodes, directed));
    return rcpp_result_gen;
END_RCPP
}
// printGraph
void printGraph(const Rcpp::List& graph);
RcppExport SEXP _rCausalMGM_printGraph(SEXP graphSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph(graphSEXP);
    printGraph(graph);
    return R_NilValue;
END_RCPP
}
// mgm
Rcpp::List mgm(const Rcpp::DataFrame& df, Rcpp::NumericVector lambda, const int maxDiscrete, const bool verbose);
RcppExport SEXP _rCausalMGM_mgm(SEXP dfSEXP, SEXP lambdaSEXP, SEXP maxDiscreteSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mgm(df, lambda, maxDiscrete, verbose));
    return rcpp_result_gen;
END_RCPP
}
// steps
Rcpp::List steps(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::NumericVector> lambda, const double g, const int numSub, const int subSize, const bool leaveOneOut, const bool computeStabs, const int threads, const bool verbose);
RcppExport SEXP _rCausalMGM_steps(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP lambdaSEXP, SEXP gSEXP, SEXP numSubSEXP, SEXP subSizeSEXP, SEXP leaveOneOutSEXP, SEXP computeStabsSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type g(gSEXP);
    Rcpp::traits::input_parameter< const int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< const int >::type subSize(subSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type leaveOneOut(leaveOneOutSEXP);
    Rcpp::traits::input_parameter< const bool >::type computeStabs(computeStabsSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(steps(df, maxDiscrete, lambda, g, numSub, subSize, leaveOneOut, computeStabs, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pcStable
Rcpp::List pcStable(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, const double alpha, const int threads, const bool fdr, const bool verbose);
RcppExport SEXP _rCausalMGM_pcStable(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcStable(df, maxDiscrete, initialGraph, alpha, threads, fdr, verbose));
    return rcpp_result_gen;
END_RCPP
}
// cpcStable
Rcpp::List cpcStable(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, const double alpha, const int threads, const bool fdr, const bool verbose);
RcppExport SEXP _rCausalMGM_cpcStable(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cpcStable(df, maxDiscrete, initialGraph, alpha, threads, fdr, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pcMax
Rcpp::List pcMax(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, const double alpha, const int threads, const bool fdr, const bool verbose);
RcppExport SEXP _rCausalMGM_pcMax(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcMax(df, maxDiscrete, initialGraph, alpha, threads, fdr, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pc50
Rcpp::List pc50(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, const double alpha, const int threads, const bool fdr, const bool verbose);
RcppExport SEXP _rCausalMGM_pc50(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pc50(df, maxDiscrete, initialGraph, alpha, threads, fdr, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fciStable
Rcpp::List fciStable(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, const double alpha, const int threads, const bool fdr, const bool verbose);
RcppExport SEXP _rCausalMGM_fciStable(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(fciStable(df, maxDiscrete, initialGraph, alpha, threads, fdr, verbose));
    return rcpp_result_gen;
END_RCPP
}
// cfci
Rcpp::List cfci(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, const double alpha, const int threads, const bool fdr, const bool verbose);
RcppExport SEXP _rCausalMGM_cfci(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cfci(df, maxDiscrete, initialGraph, alpha, threads, fdr, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fciMax
Rcpp::List fciMax(const Rcpp::DataFrame& df, const int maxDiscrete, Rcpp::Nullable<Rcpp::List> initialGraph, const double alpha, const int threads, const bool fdr, const bool verbose);
RcppExport SEXP _rCausalMGM_fciMax(SEXP dfSEXP, SEXP maxDiscreteSEXP, SEXP initialGraphSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(fciMax(df, maxDiscrete, initialGraph, alpha, threads, fdr, verbose));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap
Rcpp::List bootstrap(const Rcpp::DataFrame& df, Rcpp::StringVector method, Rcpp::StringVector ensembleMethod, Rcpp::NumericVector lambda, const double alpha, const int numBoots, const int maxDiscrete, const int threads, const bool verbose);
RcppExport SEXP _rCausalMGM_bootstrap(SEXP dfSEXP, SEXP methodSEXP, SEXP ensembleMethodSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP numBootsSEXP, SEXP maxDiscreteSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type method(methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type ensembleMethod(ensembleMethodSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type numBoots(numBootsSEXP);
    Rcpp::traits::input_parameter< const int >::type maxDiscrete(maxDiscreteSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap(df, method, ensembleMethod, lambda, alpha, numBoots, maxDiscrete, threads, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rCausalMGM_saveGraph", (DL_FUNC) &_rCausalMGM_saveGraph, 2},
    {"_rCausalMGM_loadGraph", (DL_FUNC) &_rCausalMGM_loadGraph, 1},
    {"_rCausalMGM_adjMat2Graph", (DL_FUNC) &_rCausalMGM_adjMat2Graph, 3},
    {"_rCausalMGM_printGraph", (DL_FUNC) &_rCausalMGM_printGraph, 1},
    {"_rCausalMGM_mgm", (DL_FUNC) &_rCausalMGM_mgm, 4},
    {"_rCausalMGM_steps", (DL_FUNC) &_rCausalMGM_steps, 10},
    {"_rCausalMGM_pcStable", (DL_FUNC) &_rCausalMGM_pcStable, 7},
    {"_rCausalMGM_cpcStable", (DL_FUNC) &_rCausalMGM_cpcStable, 7},
    {"_rCausalMGM_pcMax", (DL_FUNC) &_rCausalMGM_pcMax, 7},
    {"_rCausalMGM_pc50", (DL_FUNC) &_rCausalMGM_pc50, 7},
    {"_rCausalMGM_fciStable", (DL_FUNC) &_rCausalMGM_fciStable, 7},
    {"_rCausalMGM_cfci", (DL_FUNC) &_rCausalMGM_cfci, 7},
    {"_rCausalMGM_fciMax", (DL_FUNC) &_rCausalMGM_fciMax, 7},
    {"_rCausalMGM_bootstrap", (DL_FUNC) &_rCausalMGM_bootstrap, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_rCausalMGM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
