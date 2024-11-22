#include "STARS.hpp"

EdgeListGraph STARS::runStarsPar(arma::mat& instabs, arma::umat& samps) {
    // Sort in ascending order
    alphas = arma::sort(alphas, "ascend");

    int parallelism = 0;

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }

    int currIndex = 0;
    double oneAlph = -1;
    double allMax = 0;
    int allMaxI = -1;
    // bool censFlag = d.isCensored();
    bool mgmFlag = initialGraph != NULL;

    arma::vec lambdaVec;
    if (mgmFlag) {
	std::string alg = initialGraph->getAlgorithm();
	if (alg == "MGM") {
	    lambdaVec = initialGraph->getHyperParam("lambda");
	// } else if (alg == "CoxMGM") {
	//     lambdaVec = initialGraph->getHyperParam("lambda");
	} else {
	    throw std::invalid_argument("Unsupported initial graph type for cross-validation. If a method other than MGM is being used, cross-validation must be done externally in R.");
	}
    }

    if (verbose) Rcpp::Rcout << "Running STARS for " << alphas.n_elem << " alphas from "
			     << alphas[0] << " to " << alphas[alphas.n_elem-1] << "..."
			     << std::endl;
    
    if (leaveOneOut) {
	samps = StabilityUtils::subSampleLOO(d);
    } else {
	samps = StabilityUtils::subSampleNoReplacement(d, b, N);
    }
    
    int numVars = d.getNumColumns();
    arma::mat thetaMat(numVars, numVars, arma::fill::zeros);
    arma::mat thetaMatOld(numVars, numVars, arma::fill::zeros);
    std::vector<DataSet> dataSubVec;
    std::vector<IndTestMulti> itSubVec;
    
    for (int i = 0; i < samps.n_rows; i++) {
	dataSubVec.push_back(DataSet(d, samps.row(i)));
	itSubVec.push_back(IndTestMulti(dataSubVec.at(i), alphas[0]));
    }
    
    std::vector<EdgeListGraph> igVec(samps.n_rows);
    std::vector<std::future<EdgeListGraph>> futures(samps.n_rows);

    if (mgmFlag) {
      RcppThread::ThreadPool pool(std::max(1, std::min(parallelism, (int) samps.n_rows)));
        lambda = std::vector<double>(lambdaVec.begin(), lambdaVec.end());

	auto initialGraphTask = [&] (int i) {
				    EdgeListGraph ig;
				    
				    if (RcppThread::isInterrupted())
					return ig;
				    
				    // if (!censFlag) {
				    MGM mgm(dataSubVec.at(i), lambda);
				    ig = mgm.search();
				    // } else {
				    // 	CoxMGM coxmgm(dataSubVec.at(i), lambda);
				    // 	ig = coxmgm.search();
				    // }
				    return ig;
				};

	for (int i = 0; i < samps.n_rows; i++) {
	    futures[i] = pool.pushReturn(initialGraphTask, i);
	}

	for (int i = 0; i < samps.n_rows; i++) {
	    igVec[i] = futures[i].get();
	}

	RcppThread::checkUserInterrupt();

	pool.join();
    }


    // go until we break by having instability better than threshold
    for (currIndex = 0; currIndex < alphas.n_elem; currIndex++) {
	
        if (verbose) Rcpp::Rcout << "  Testing alpha = " << alphas[currIndex] << std::endl;

	thetaMatOld = thetaMat;

	thetaMat.fill(0);

	for (int i = 0; i < samps.n_rows; i++) {
	    
	    itSubVec[i].setAlpha(alphas[currIndex]);
	    if (alg == "pc") {
		PcStable causalAlg((IndependenceTest*) &itSubVec[i]);

		if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		causalAlg.setKnowledge(knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(fdr);
		causalAlg.setOrientRule(ORIENT_SEPSETS);
	    
		if (mgmFlag) {
		    causalAlg.setInitialGraph(&igVec[i]);
		}
	    
		EdgeListGraph g = causalAlg.search();

		thetaMat += StabilityUtils::skeletonToMatrix(g, d);
		
	    } else if (alg == "fci") {
		Fci causalAlg((IndependenceTest*) &itSubVec[i]);

		if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		causalAlg.setKnowledge(knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(fdr);
		causalAlg.setOrientRule(orientRule);
	    
		if (mgmFlag) {
		    causalAlg.setInitialGraph(&igVec[i]);
		}
	    
		EdgeListGraph g = causalAlg.search();

		thetaMat += StabilityUtils::skeletonToMatrix(g, d);
	    }
	}

        arma::mat adjMat = thetaMat / samps.n_rows;

	// Rcpp::Rcout << "adjMat:\n" << adjMat << std::endl;

	double numAll = numVars * (numVars-1);
        double allDestable = arma::accu(2 * adjMat % (1-adjMat)) / numAll;

	// Rcpp::Rcout << "numAll:  " << numAll << std::endl;

	instabs(currIndex, 0) = std::max(allMax, allDestable);
        
	if (verbose) {
	    Rcpp::Rcout << "  Overall instability for alpha = " << alphas[currIndex]
			<< ":  " << allDestable << std::endl;
	}

	if (allDestable >= gamma && oneAlph == -1 && currIndex > 0) {
	    oneAlph = alphas[currIndex-1];
        }
	
        if (allDestable >= allMax) {
            allMax = allDestable;
            allMaxI = currIndex;
        }

	if (oneAlph != -1) {
	    break;
	}

    }
    // EXIT_LOOP:

    if (oneAlph == -1) {
        oneAlph = alphas[allMaxI];
	stabilities = thetaMat / samps.n_rows;
    } else {
	stabilities = thetaMatOld / samps.n_rows;
    }

    if (verbose) Rcpp::Rcout << "StARS alpha: " << oneAlph << std::endl;
    
    // EdgeListGraph ig, g;
    // if (mgmFlag) {
    // 	if (censFlag) {
    // 	    CoxMGM m(d, lambda);
    // 	    ig = m.search();
    // 	} else {
    // 	    MGM m(d, lambda);
    // 	    ig = m.search();
    // 	}
    // }

    EdgeListGraph g;

    IndTestMulti itm(d, oneAlph);
    
    if (alg == "pc") {
	PcStable causalAlg((IndependenceTest*) &itm);

	if (parallelism > 0) causalAlg.setThreads(parallelism);
	
	causalAlg.setKnowledge(knowledge);
	causalAlg.setVerbose(false);
	causalAlg.setFDR(fdr);
	causalAlg.setOrientRule(orientRule);
	    
	if (mgmFlag) {
	    causalAlg.setInitialGraph(initialGraph);
	}
	    
	g = causalAlg.search();
		
    } else if (alg == "fci") {
	Fci causalAlg((IndependenceTest*) &itm);

	if (parallelism > 0) causalAlg.setThreads(parallelism);
	
	causalAlg.setKnowledge(knowledge);
	causalAlg.setVerbose(false);
	causalAlg.setFDR(fdr);
	causalAlg.setOrientRule(orientRule);
	
	if (mgmFlag) {
	    causalAlg.setInitialGraph(initialGraph);
	}
	    
	g = causalAlg.search();
    }

    
    // if (computeStabs) {
    // 	if (allMaxI
    // }
    // lastLambda = lambda;
    
    // EdgeListGraph g = m.graphFromMGM();
    // gSteps.setHyperParam("lambda", arma::vec(lambda));
    
    return g;    
}
