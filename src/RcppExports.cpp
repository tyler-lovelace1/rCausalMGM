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
void saveGraph(const Rcpp::List& graph, const std::string& filename);
RcppExport SEXP _rCausalMGM_saveGraph(SEXP graphSEXP, SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    saveGraph(graph, filename);
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
// cpdag
Rcpp::List cpdag(const Rcpp::List& graph);
RcppExport SEXP _rCausalMGM_cpdag(SEXP graphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph(graphSEXP);
    rcpp_result_gen = Rcpp::wrap(cpdag(graph));
    return rcpp_result_gen;
END_RCPP
}
// moral
Rcpp::List moral(const Rcpp::List& graph);
RcppExport SEXP _rCausalMGM_moral(SEXP graphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph(graphSEXP);
    rcpp_result_gen = Rcpp::wrap(moral(graph));
    return rcpp_result_gen;
END_RCPP
}
// skeleton
Rcpp::List skeleton(const Rcpp::List& graph);
RcppExport SEXP _rCausalMGM_skeleton(SEXP graphSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph(graphSEXP);
    rcpp_result_gen = Rcpp::wrap(skeleton(graph));
    return rcpp_result_gen;
END_RCPP
}
// pag
Rcpp::List pag(const Rcpp::List& graph, Rcpp::Nullable<Rcpp::StringVector> latent);
RcppExport SEXP _rCausalMGM_pag(SEXP graphSEXP, SEXP latentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::StringVector> >::type latent(latentSEXP);
    rcpp_result_gen = Rcpp::wrap(pag(graph, latent));
    return rcpp_result_gen;
END_RCPP
}
// SHD
Rcpp::NumericVector SHD(const Rcpp::List& graph1, const Rcpp::List& graph2);
RcppExport SEXP _rCausalMGM_SHD(SEXP graph1SEXP, SEXP graph2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph1(graph1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graph2(graph2SEXP);
    rcpp_result_gen = Rcpp::wrap(SHD(graph1, graph2));
    return rcpp_result_gen;
END_RCPP
}
// prMetricsAdjacency
Rcpp::NumericVector prMetricsAdjacency(const Rcpp::List& estimate, const Rcpp::List& groundTruth);
RcppExport SEXP _rCausalMGM_prMetricsAdjacency(SEXP estimateSEXP, SEXP groundTruthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type estimate(estimateSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type groundTruth(groundTruthSEXP);
    rcpp_result_gen = Rcpp::wrap(prMetricsAdjacency(estimate, groundTruth));
    return rcpp_result_gen;
END_RCPP
}
// prMetricsOrientation
Rcpp::NumericVector prMetricsOrientation(const Rcpp::List& estimate, const Rcpp::List& groundTruth, const Rcpp::Nullable<Rcpp::List>& groundTruthDAG);
RcppExport SEXP _rCausalMGM_prMetricsOrientation(SEXP estimateSEXP, SEXP groundTruthSEXP, SEXP groundTruthDAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type estimate(estimateSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type groundTruth(groundTruthSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::List>& >::type groundTruthDAG(groundTruthDAGSEXP);
    rcpp_result_gen = Rcpp::wrap(prMetricsOrientation(estimate, groundTruth, groundTruthDAG));
    return rcpp_result_gen;
END_RCPP
}
// prMetricsCausal
Rcpp::NumericVector prMetricsCausal(const Rcpp::List& estimate, const Rcpp::List& groundTruthDAG);
RcppExport SEXP _rCausalMGM_prMetricsCausal(SEXP estimateSEXP, SEXP groundTruthDAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type estimate(estimateSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type groundTruthDAG(groundTruthDAGSEXP);
    rcpp_result_gen = Rcpp::wrap(prMetricsCausal(estimate, groundTruthDAG));
    return rcpp_result_gen;
END_RCPP
}
// prMetrics
Rcpp::NumericVector prMetrics(const Rcpp::List& estimate, const Rcpp::List& groundTruth, const Rcpp::Nullable<Rcpp::List>& groundTruthDAG);
RcppExport SEXP _rCausalMGM_prMetrics(SEXP estimateSEXP, SEXP groundTruthSEXP, SEXP groundTruthDAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type estimate(estimateSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type groundTruth(groundTruthSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::List>& >::type groundTruthDAG(groundTruthDAGSEXP);
    rcpp_result_gen = Rcpp::wrap(prMetrics(estimate, groundTruth, groundTruthDAG));
    return rcpp_result_gen;
END_RCPP
}
// allMetrics
Rcpp::NumericVector allMetrics(const Rcpp::List& estimate, const Rcpp::List& groundTruth, const Rcpp::Nullable<Rcpp::List>& groundTruthDAG);
RcppExport SEXP _rCausalMGM_allMetrics(SEXP estimateSEXP, SEXP groundTruthSEXP, SEXP groundTruthDAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type estimate(estimateSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type groundTruth(groundTruthSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Rcpp::List>& >::type groundTruthDAG(groundTruthDAGSEXP);
    rcpp_result_gen = Rcpp::wrap(allMetrics(estimate, groundTruth, groundTruthDAG));
    return rcpp_result_gen;
END_RCPP
}
// mgm
Rcpp::List mgm(const Rcpp::DataFrame& data, Rcpp::NumericVector lambda, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_mgm(SEXP dataSEXP, SEXP lambdaSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mgm(data, lambda, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// coxmgm
Rcpp::List coxmgm(const Rcpp::DataFrame& data, Rcpp::NumericVector lambda, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_coxmgm(SEXP dataSEXP, SEXP lambdaSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(coxmgm(data, lambda, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// mgmPath
Rcpp::List mgmPath(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::NumericVector> lambdas, const int nLambda, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_mgmPath(SEXP dataSEXP, SEXP lambdasSEXP, SEXP nLambdaSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const int >::type nLambda(nLambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mgmPath(data, lambdas, nLambda, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// coxmgmPath
Rcpp::List coxmgmPath(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::NumericVector> lambdas, const int nLambda, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_coxmgmPath(SEXP dataSEXP, SEXP lambdasSEXP, SEXP nLambdaSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const int >::type nLambda(nLambdaSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(coxmgmPath(data, lambdas, nLambda, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// mgmCV
Rcpp::List mgmCV(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::NumericVector> lambdas, const int nLambda, const int nfolds, Rcpp::Nullable<Rcpp::NumericVector> foldid, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_mgmCV(SEXP dataSEXP, SEXP lambdasSEXP, SEXP nLambdaSEXP, SEXP nfoldsSEXP, SEXP foldidSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const int >::type nLambda(nLambdaSEXP);
    Rcpp::traits::input_parameter< const int >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type foldid(foldidSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mgmCV(data, lambdas, nLambda, nfolds, foldid, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// steps
Rcpp::List steps(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::NumericVector> lambdas, const int nLambda, const double gamma, const int numSub, const int subSize, const bool leaveOneOut, const int threads, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_steps(SEXP dataSEXP, SEXP lambdasSEXP, SEXP nLambdaSEXP, SEXP gammaSEXP, SEXP numSubSEXP, SEXP subSizeSEXP, SEXP leaveOneOutSEXP, SEXP threadsSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const int >::type nLambda(nLambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< const int >::type subSize(subSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type leaveOneOut(leaveOneOutSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(steps(data, lambdas, nLambda, gamma, numSub, subSize, leaveOneOut, threads, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pcStars
Rcpp::List pcStars(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::Nullable<Rcpp::List> knowledge, const Rcpp::StringVector orientRule, Rcpp::Nullable<Rcpp::NumericVector> alphas, const double gamma, const int numSub, const int subSize, const bool leaveOneOut, const int threads, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_pcStars(SEXP dataSEXP, SEXP initialGraphSEXP, SEXP knowledgeSEXP, SEXP orientRuleSEXP, SEXP alphasSEXP, SEXP gammaSEXP, SEXP numSubSEXP, SEXP subSizeSEXP, SEXP leaveOneOutSEXP, SEXP threadsSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< const int >::type subSize(subSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type leaveOneOut(leaveOneOutSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcStars(data, initialGraph, knowledge, orientRule, alphas, gamma, numSub, subSize, leaveOneOut, threads, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fciStars
Rcpp::List fciStars(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::Nullable<Rcpp::List> knowledge, const Rcpp::StringVector orientRule, Rcpp::Nullable<Rcpp::NumericVector> alphas, const double gamma, const int numSub, const int subSize, const bool leaveOneOut, const int threads, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_fciStars(SEXP dataSEXP, SEXP initialGraphSEXP, SEXP knowledgeSEXP, SEXP orientRuleSEXP, SEXP alphasSEXP, SEXP gammaSEXP, SEXP numSubSEXP, SEXP subSizeSEXP, SEXP leaveOneOutSEXP, SEXP threadsSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< const int >::type subSize(subSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type leaveOneOut(leaveOneOutSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(fciStars(data, initialGraph, knowledge, orientRule, alphas, gamma, numSub, subSize, leaveOneOut, threads, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pcStable
Rcpp::List pcStable(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::Nullable<Rcpp::List> knowledge, const Rcpp::StringVector orientRule, const double alpha, const int threads, const bool fdr, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_pcStable(SEXP dataSEXP, SEXP initialGraphSEXP, SEXP knowledgeSEXP, SEXP orientRuleSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcStable(data, initialGraph, knowledge, orientRule, alpha, threads, fdr, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fciStable
Rcpp::List fciStable(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::Nullable<Rcpp::List> knowledge, const Rcpp::StringVector orientRule, const double alpha, const int threads, const bool possDsep, const bool fdr, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_fciStable(SEXP dataSEXP, SEXP initialGraphSEXP, SEXP knowledgeSEXP, SEXP orientRuleSEXP, SEXP alphaSEXP, SEXP threadsSEXP, SEXP possDsepSEXP, SEXP fdrSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type possDsep(possDsepSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(fciStable(data, initialGraph, knowledge, orientRule, alpha, threads, possDsep, fdr, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pcCV
Rcpp::List pcCV(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::Nullable<Rcpp::List> knowledge, const Rcpp::StringVector orientRule, Rcpp::Nullable<Rcpp::NumericVector> alphas, const int nfolds, Rcpp::Nullable<Rcpp::NumericVector> foldid, const int threads, const bool fdr, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_pcCV(SEXP dataSEXP, SEXP initialGraphSEXP, SEXP knowledgeSEXP, SEXP orientRuleSEXP, SEXP alphasSEXP, SEXP nfoldsSEXP, SEXP foldidSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const int >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type foldid(foldidSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcCV(data, initialGraph, knowledge, orientRule, alphas, nfolds, foldid, threads, fdr, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fciCV
Rcpp::List fciCV(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> initialGraph, Rcpp::Nullable<Rcpp::List> knowledge, const Rcpp::StringVector orientRule, Rcpp::Nullable<Rcpp::NumericVector> alphas, const int nfolds, Rcpp::Nullable<Rcpp::NumericVector> foldid, const int threads, const bool fdr, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_fciCV(SEXP dataSEXP, SEXP initialGraphSEXP, SEXP knowledgeSEXP, SEXP orientRuleSEXP, SEXP alphasSEXP, SEXP nfoldsSEXP, SEXP foldidSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type initialGraph(initialGraphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const int >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type foldid(foldidSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(fciCV(data, initialGraph, knowledge, orientRule, alphas, nfolds, foldid, threads, fdr, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// mgmpcCV
Rcpp::List mgmpcCV(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> knowledge, const std::string cvType, const Rcpp::StringVector orientRule, Rcpp::Nullable<Rcpp::NumericVector> lambdas, const int nLambda, Rcpp::Nullable<Rcpp::NumericVector> alphas, const int numPoints, const int nfolds, Rcpp::Nullable<Rcpp::NumericVector> foldid, const int threads, const bool fdr, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_mgmpcCV(SEXP dataSEXP, SEXP knowledgeSEXP, SEXP cvTypeSEXP, SEXP orientRuleSEXP, SEXP lambdasSEXP, SEXP nLambdaSEXP, SEXP alphasSEXP, SEXP numPointsSEXP, SEXP nfoldsSEXP, SEXP foldidSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type cvType(cvTypeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const int >::type nLambda(nLambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const int >::type numPoints(numPointsSEXP);
    Rcpp::traits::input_parameter< const int >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type foldid(foldidSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mgmpcCV(data, knowledge, cvType, orientRule, lambdas, nLambda, alphas, numPoints, nfolds, foldid, threads, fdr, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// mgmfciCV
Rcpp::List mgmfciCV(const Rcpp::DataFrame& data, Rcpp::Nullable<Rcpp::List> knowledge, const std::string cvType, const Rcpp::StringVector orientRule, Rcpp::Nullable<Rcpp::NumericVector> lambdas, const int nLambda, Rcpp::Nullable<Rcpp::NumericVector> alphas, const int numPoints, const int nfolds, Rcpp::Nullable<Rcpp::NumericVector> foldid, const int threads, const bool fdr, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_mgmfciCV(SEXP dataSEXP, SEXP knowledgeSEXP, SEXP cvTypeSEXP, SEXP orientRuleSEXP, SEXP lambdasSEXP, SEXP nLambdaSEXP, SEXP alphasSEXP, SEXP numPointsSEXP, SEXP nfoldsSEXP, SEXP foldidSEXP, SEXP threadsSEXP, SEXP fdrSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type cvType(cvTypeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type orientRule(orientRuleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const int >::type nLambda(nLambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const int >::type numPoints(numPointsSEXP);
    Rcpp::traits::input_parameter< const int >::type nfolds(nfoldsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type foldid(foldidSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type fdr(fdrSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(mgmfciCV(data, knowledge, cvType, orientRule, lambdas, nLambda, alphas, numPoints, nfolds, foldid, threads, fdr, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap
Rcpp::List bootstrap(const Rcpp::DataFrame& data, Rcpp::List graph, Rcpp::Nullable<Rcpp::List> knowledge, const int numBoots, const int threads, const bool replace, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_bootstrap(SEXP dataSEXP, SEXP graphSEXP, SEXP knowledgeSEXP, SEXP numBootsSEXP, SEXP threadsSEXP, SEXP replaceSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type knowledge(knowledgeSEXP);
    Rcpp::traits::input_parameter< const int >::type numBoots(numBootsSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type replace(replaceSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap(data, graph, knowledge, numBoots, threads, replace, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// growShrinkMB
Rcpp::StringVector growShrinkMB(const Rcpp::DataFrame& data, const std::string& target, const double penalty, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_growShrinkMB(SEXP dataSEXP, SEXP targetSEXP, SEXP penaltySEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(growShrinkMB(data, target, penalty, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// grasp
Rcpp::List grasp(const Rcpp::DataFrame& data, const int depth, const int numStarts, const double penalty, const bool bossInit, const int threads, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_grasp(SEXP dataSEXP, SEXP depthSEXP, SEXP numStartsSEXP, SEXP penaltySEXP, SEXP bossInitSEXP, SEXP threadsSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type depth(depthSEXP);
    Rcpp::traits::input_parameter< const int >::type numStarts(numStartsSEXP);
    Rcpp::traits::input_parameter< const double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< const bool >::type bossInit(bossInitSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(grasp(data, depth, numStarts, penalty, bossInit, threads, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}
// boss
Rcpp::List boss(const Rcpp::DataFrame& data, const int numStarts, const double penalty, const int threads, const bool rank, const bool verbose);
RcppExport SEXP _rCausalMGM_boss(SEXP dataSEXP, SEXP numStartsSEXP, SEXP penaltySEXP, SEXP threadsSEXP, SEXP rankSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type numStarts(numStartsSEXP);
    Rcpp::traits::input_parameter< const double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(boss(data, numStarts, penalty, threads, rank, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rCausalMGM_saveGraph", (DL_FUNC) &_rCausalMGM_saveGraph, 2},
    {"_rCausalMGM_loadGraph", (DL_FUNC) &_rCausalMGM_loadGraph, 1},
    {"_rCausalMGM_adjMat2Graph", (DL_FUNC) &_rCausalMGM_adjMat2Graph, 3},
    {"_rCausalMGM_printGraph", (DL_FUNC) &_rCausalMGM_printGraph, 1},
    {"_rCausalMGM_cpdag", (DL_FUNC) &_rCausalMGM_cpdag, 1},
    {"_rCausalMGM_moral", (DL_FUNC) &_rCausalMGM_moral, 1},
    {"_rCausalMGM_skeleton", (DL_FUNC) &_rCausalMGM_skeleton, 1},
    {"_rCausalMGM_pag", (DL_FUNC) &_rCausalMGM_pag, 2},
    {"_rCausalMGM_SHD", (DL_FUNC) &_rCausalMGM_SHD, 2},
    {"_rCausalMGM_prMetricsAdjacency", (DL_FUNC) &_rCausalMGM_prMetricsAdjacency, 2},
    {"_rCausalMGM_prMetricsOrientation", (DL_FUNC) &_rCausalMGM_prMetricsOrientation, 3},
    {"_rCausalMGM_prMetricsCausal", (DL_FUNC) &_rCausalMGM_prMetricsCausal, 2},
    {"_rCausalMGM_prMetrics", (DL_FUNC) &_rCausalMGM_prMetrics, 3},
    {"_rCausalMGM_allMetrics", (DL_FUNC) &_rCausalMGM_allMetrics, 3},
    {"_rCausalMGM_mgm", (DL_FUNC) &_rCausalMGM_mgm, 4},
    {"_rCausalMGM_coxmgm", (DL_FUNC) &_rCausalMGM_coxmgm, 4},
    {"_rCausalMGM_mgmPath", (DL_FUNC) &_rCausalMGM_mgmPath, 5},
    {"_rCausalMGM_coxmgmPath", (DL_FUNC) &_rCausalMGM_coxmgmPath, 5},
    {"_rCausalMGM_mgmCV", (DL_FUNC) &_rCausalMGM_mgmCV, 7},
    {"_rCausalMGM_steps", (DL_FUNC) &_rCausalMGM_steps, 10},
    {"_rCausalMGM_pcStars", (DL_FUNC) &_rCausalMGM_pcStars, 12},
    {"_rCausalMGM_fciStars", (DL_FUNC) &_rCausalMGM_fciStars, 12},
    {"_rCausalMGM_pcStable", (DL_FUNC) &_rCausalMGM_pcStable, 9},
    {"_rCausalMGM_fciStable", (DL_FUNC) &_rCausalMGM_fciStable, 10},
    {"_rCausalMGM_pcCV", (DL_FUNC) &_rCausalMGM_pcCV, 11},
    {"_rCausalMGM_fciCV", (DL_FUNC) &_rCausalMGM_fciCV, 11},
    {"_rCausalMGM_mgmpcCV", (DL_FUNC) &_rCausalMGM_mgmpcCV, 14},
    {"_rCausalMGM_mgmfciCV", (DL_FUNC) &_rCausalMGM_mgmfciCV, 14},
    {"_rCausalMGM_bootstrap", (DL_FUNC) &_rCausalMGM_bootstrap, 8},
    {"_rCausalMGM_growShrinkMB", (DL_FUNC) &_rCausalMGM_growShrinkMB, 5},
    {"_rCausalMGM_grasp", (DL_FUNC) &_rCausalMGM_grasp, 8},
    {"_rCausalMGM_boss", (DL_FUNC) &_rCausalMGM_boss, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_rCausalMGM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
