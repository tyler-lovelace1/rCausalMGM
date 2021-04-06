#include "Tests.hpp"

#include "MGM.hpp"
#include "IndTestMulti.hpp"
#include "SepsetMap.hpp"
#include "ChoiceGenerator.hpp"
#include "PcStable.hpp"
#include "CpcStable.hpp"
#include "PcMax.hpp"
#include "BlockingQueue.hpp"
#include "STEPS.hpp"
#include <thread>
#include <atomic>
#include <cstdlib>


void Tests::testConcurrentQueue() {
    BlockingQueue<int> q(100);

    std::mutex cout_lock;

    int NTHREADS = 32;

    int COUNT = NTHREADS * 10000;

    std::atomic<int> countDownLatch(NTHREADS);

    std::vector<bool> present(COUNT);
    for (int i = 0; i < COUNT; i++)
	present[i] = false;

    std::vector<int> outputs[NTHREADS];

    std::thread producers[NTHREADS];

    std::thread consumers[NTHREADS];

    std::atomic<int> prodCount(0);

    std::atomic<int> consCount(0);

    for (int i = 0; i < NTHREADS; i++) {
    	producers[i] = std::thread([&]() {
    		int id = prodCount.fetch_add(1, std::memory_order_relaxed);
    		for(int j = id * COUNT / NTHREADS + 1; j <= (id+1) * COUNT / NTHREADS; ++j) {
    		    q.push(j);
    		}
    		countDownLatch--;
    		if (countDownLatch == 0) {
    		    q.push(-1);
    		}
    	    });
    }

    for (int i = 0; i < NTHREADS; i++) {
	consumers[i] = std::thread([&]() {
		int id = consCount.fetch_add(1, std::memory_order_relaxed);
		int v = q.pop();
		while (v != -1) {
		    outputs[id].push_back(v);
		    std::this_thread::sleep_for(std::chrono::nanoseconds(100));
		    v = q.pop();
		}
		q.push(-1);
		q.push(-1);
	    });
    }
    

    for (int i = 0; i < NTHREADS; i++)
    	producers[i].join();
    for (int i = 0; i < NTHREADS; i++)
	consumers[i].join();

    for (int i = 0; i < NTHREADS; i++) {
	Rcpp::Rcout << "\nconsumer " << i+1 << " size: " << outputs[i].size() << "\n";
	for (int j = 0; j < outputs[i].size(); j++) {
	    // Rcpp::Rcout << outputs[i][j] << " ";
	    present[outputs[i][j]-1] = true;
	}
	// Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";

    bool missing;
    for (int i = 0; i < COUNT; i++) {
	missing = !present[i];
	if (missing)
	    break;
    }

    if (missing)
	Rcpp::Rcout << "Value missing from output\n";
    else
	Rcpp::Rcout << "All values present in output\n";
}

void Tests::testMGMFunctions(const Rcpp::DataFrame &df, const int maxDiscrete) {
    DataSet ds(df, maxDiscrete);

    std::vector<double> lambda = {0.2, 0.2, 0.2};

    MGM mgm(ds, lambda);

    arma::vec params = mgm.params.toMatrix1D();

    Rcpp::Rcout << "Smooth value: \n" << mgm.smoothValue(params) << std::endl;

    arma::vec smoothGrad = mgm.smoothGradient(params);
    Rcpp::Rcout << "Calculated smooth gradient" << std::endl;
    MGMParams params2(smoothGrad, 5, 20);
    Rcpp::Rcout << "Smooth gradient: \n" << params2;

    // Test smooth()
    arma::vec smoothOut(arma::size(smoothGrad));
    Rcpp::Rcout << "Calling smooth()..." << std::endl;
    Rcpp::Rcout << "Smooth() value: \n" << mgm.smooth(params, smoothOut) << std::endl;
    MGMParams params3(smoothOut, 5, 20);
    Rcpp::Rcout << "Smooth() gradient: \n" << params3;

    arma::vec proxOperator = mgm.proximalOperator(1, smoothGrad);
    MGMParams params4(proxOperator, 5, 20);
    Rcpp::Rcout << "proxOperator: \n" << params4;

    Rcpp::Rcout << "nonSmoothValue: \n" << mgm.nonSmoothValue(smoothGrad) << std::endl;

    // Test nonsmooth()
    arma::vec nonSmoothOut(arma::size(smoothGrad));
    Rcpp::Rcout << "Calling nonSmooth()..." << std::endl;
    Rcpp::Rcout << "nonSmooth() value: \n" << mgm.nonSmooth(1, smoothGrad, nonSmoothOut) << std::endl;
    MGMParams params5(nonSmoothOut, 5, 20);
    Rcpp::Rcout << "nonSmooth() gradient: \n" << params5;

    mgm.learnEdges(500);

    Rcpp::Rcout << "mgm.params.alpha1: \n" << mgm.params.getAlpha1() << std::endl;
    Rcpp::Rcout << "mgm.params.alpha2: \n" << mgm.params.getAlpha2() << std::endl;
    Rcpp::Rcout << "mgm.params.betad: \n" << mgm.params.getBetad() << std::endl;
    Rcpp::Rcout << "mgm.params.beta: \n" << mgm.params.getBeta() << std::endl;
    Rcpp::Rcout << "mgm.params.theta: \n" << mgm.params.getTheta() << std::endl;
    Rcpp::Rcout << "mgm.params.phi: \n" << mgm.params.getPhi() << std::endl;

    EdgeListGraph mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time = " << mgm.getElapsedTime() << " ms" << std::endl;

    Rcpp::Rcout << "MGM GRAPH\n" << mgmGraph << std::endl;
}

void Tests::testMGMTiming(const Rcpp::DataFrame &df, const int maxDiscrete) {
    DataSet ds(df, maxDiscrete);

    std::vector<double> lambda = {0.05, 0.05, 0.05};
    MGM mgm(ds, lambda);
    EdgeListGraph mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time (l: " << lambda[0] << ") = " << mgm.getElapsedTime() << " ms" << std::endl;

    lambda = {0.1, 0.1, 0.1};
    mgm = MGM(ds, lambda);
    mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time (l: " << lambda[0] << ") = " << mgm.getElapsedTime() << " ms" << std::endl;

    lambda = {0.15, 0.15, 0.15};
    mgm = MGM(ds, lambda);
    mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time (l: " << lambda[0] << ") = " << mgm.getElapsedTime() << " ms" << std::endl;

    lambda = {0.2, 0.2, 0.2};
    mgm = MGM(ds, lambda);
    mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time (l: " << lambda[0] << ") = " << mgm.getElapsedTime() << " ms" << std::endl;

    lambda = {0.25, 0.25, 0.25};
    mgm = MGM(ds, lambda);
    mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time (l: " << lambda[0] << ") = " << mgm.getElapsedTime() << " ms" << std::endl;

    lambda = {0.3, 0.3, 0.3};
    mgm = MGM(ds, lambda);
    mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time (l: " << lambda[0] << ") = " << mgm.getElapsedTime() << " ms" << std::endl;
}

void Tests::testPcStable(const Rcpp::DataFrame &df, const int maxDiscrete) {
    DataSet ds(df, maxDiscrete);
    std::vector<double> lambda = {0.2, 0.2, 0.2};
    MGM mgm(ds, lambda);
    EdgeListGraph mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time = " << mgm.getElapsedTime() << " ms" << std::endl;
    Rcpp::Rcout << "MGM GRAPH\n" << mgmGraph << std::endl;
    IndTestMulti itm(ds, 0.05);

    PcStable pcs((IndependenceTest*) &itm);
    pcs.setInitialGraph(&mgmGraph);
    EdgeListGraph pcGraph = pcs.search();
    Rcpp::Rcout << "PC GRAPH\n" << pcGraph << std::endl;
}

void Tests::testCpcStable(const Rcpp::DataFrame &df, const int maxDiscrete) {
    DataSet ds(df, maxDiscrete);
    std::vector<double> lambda = {0.2, 0.2, 0.2};
    MGM mgm(ds, lambda);
    EdgeListGraph mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time = " << mgm.getElapsedTime() << " ms" << std::endl;
    Rcpp::Rcout << "MGM GRAPH\n" << mgmGraph << std::endl;
    IndTestMulti itm(ds, 0.05);

    CpcStable cpcs((IndependenceTest*) &itm);
    cpcs.setInitialGraph(&mgmGraph);
    EdgeListGraph cpcGraph = cpcs.search();
    Rcpp::Rcout << "CPC GRAPH\n" << cpcGraph << std::endl;
    Rcpp::Rcout << "Ambiguous triples (i.e. list of triples for which there is ambiguous data about whether they are colliders or not):" << std::endl;
    for (Triple t : cpcGraph.getAmbiguousTriples()) Rcpp::Rcout << t << std::endl;

}

void Tests::testPcMax(const Rcpp::DataFrame &df, const int maxDiscrete) {
    DataSet ds(df, maxDiscrete);
    std::vector<double> lambda = {0.2, 0.2, 0.2};
    MGM mgm(ds, lambda);
    EdgeListGraph mgmGraph = mgm.search();
    Rcpp::Rcout << "MGM elapsed time = " << mgm.getElapsedTime() << " ms" << std::endl;
    Rcpp::Rcout << "MGM GRAPH\n" << mgmGraph << std::endl;
    IndTestMulti itm(ds, 0.05);

    PcMax pcm((IndependenceTest*) &itm);
    pcm.setInitialGraph(&mgmGraph);
    EdgeListGraph pcmGraph = pcm.search();
    Rcpp::Rcout << "PCM GRAPH\n" << pcmGraph << std::endl;
}

void Tests::testSTEPS(const Rcpp::DataFrame &df, const int maxDiscrete) {
    DataSet ds(df, maxDiscrete);
    std::vector<double> lambdas = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
                                   0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85};

    STEPS steps(ds, lambdas, 0.05, 50);
    steps.setVerbose(true);

    EdgeListGraph g = steps.runStepsPar();

    Rcpp::Rcout << "STEPS graph:\n" << g << std::endl;
}