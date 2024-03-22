#ifndef STARS_HPP_
#define STARS_HPP_

#include "armaLapack.hpp"
#include "StabilityUtils.hpp"
#include "PcStable.hpp"
#include "Fci.hpp"
#include "MGM.hpp"
#include "CoxMGM.hpp"
#include "IndTestMultiCox.hpp"

class STARS {

private:

    DataSet d;
    int N;
    int b;

    std::vector<double> lambda;
    arma::vec alphas;
    EdgeListGraph *initialGraph = NULL;
    OrientRule orientRule = ORIENT_MAJORITY;
    std::string alg = "pc";
    Knowledge knowledge;
    bool fdr = false;

    double gamma;
    int iterLimit = 500;
    bool leaveOneOut = false;
    arma::mat stabilities;
    bool computeStabs = false;
    bool verbose = false;
    int threads = -1;

public:
    STARS() {}
    
    STARS(DataSet& dat, std::string alg, arma::vec& alphas, double g, int numSub, bool loo = false) :
        d(dat),
        leaveOneOut(loo),
        N(numSub),
        gamma(g),
	alg(alg),
        alphas(alphas),
        b(StabilityUtils::getSubSize(dat.getNumRows())) {}

    STARS(DataSet& dat, std::string alg, arma::vec& alphas, double g, int numSub, int subSize, bool loo = false) :
        d(dat),
        leaveOneOut(loo),
        N(numSub),
        gamma(g),
	alg(alg),
        alphas(alphas),
        b(subSize) {}

    EdgeListGraph runStarsPar(arma::mat& instabs, arma::umat& samps);

    void setVerbose(bool v) { verbose = v; }
    void setThreads(int threads) { this->threads = threads; }

    void setInitialGraph(EdgeListGraph *initialGraph) { this->initialGraph = initialGraph; }
    void setAlphas(arma::vec alphas) { this->alphas = alphas; }
    
    void setLambda(std::vector<double> l) { lambda = l; }
    void setAlgorithm(std::string alg) { this->alg = alg; }
    void setOrientRule(OrientRule rule) {this->orientRule = rule; }
    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }

    void setComputeStabs(bool cs) { computeStabs = cs; }
    bool getComputeStabs() { return computeStabs; }

    arma::mat getStabs() { return stabilities; }
};

#endif /* STARS_HPP_ */
