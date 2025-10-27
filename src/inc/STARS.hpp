#ifndef STARS_HPP_
#define STARS_HPP_

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
    arma::vec alphas;
    double gamma;
    std::string alg = "pc";
    bool leaveOneOut = false;

    
    EdgeListGraph *initialGraph = NULL;
    OrientRule orientRule = ORIENT_MAJORITY;
    Knowledge knowledge;
    bool fdr = false;
    std::vector<double> lambda;
    // int iterLimit = 500;
    arma::mat stabilities;
    bool computeStabs = false;
    bool verbose = false;
    int threads = -1;

    static const std::size_t MAX_QUEUE_SIZE = 100;
    BlockingQueue<ParallelTask> taskQueue;

    void parallelTaskConsumer() {
	while (true) {
	    ParallelTask t = taskQueue.pop();

	    if (t.is_poison())
		break;

	    if (RcppThread::isInterrupted())
		break;	

	    t();
	}
    }

public:
    STARS() : taskQueue(MAX_QUEUE_SIZE) {}
    
    STARS(DataSet& dat, std::string alg, arma::vec& alphas, double g, int numSub, bool loo = false) :
        d(dat),
        N(numSub),
	b(StabilityUtils::getSubSize(dat.getNumRows())),
        alphas(alphas),
        gamma(g),
	alg(alg),
        leaveOneOut(loo),
	taskQueue(MAX_QUEUE_SIZE) {}

    STARS(DataSet& dat, std::string alg, arma::vec& alphas, double g, int numSub, int subSize, bool loo = false) :
	d(dat),
        N(numSub),
	b(subSize),
        alphas(alphas),
        gamma(g),
	alg(alg),
        leaveOneOut(loo),
	taskQueue(MAX_QUEUE_SIZE) {}

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
