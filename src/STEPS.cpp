#include "STEPS.hpp"

EdgeListGraph STEPS::runStepsPar() {

    // Sort in descending order
    std::sort(lambda.begin(), lambda.end(), std::greater<double>());

    int currIndex = 0;
    double CC = -1;
    double CD = -1;
    double DD = -1;
    double CCMax = 0;
    double CCMaxI = -1;
    double CDMax = 0;
    double CDMaxI = -1;
    double DDMax = 0;
    double DDMaxI = -1;
    double oneLamb = -1;
    double allMax = 0;
    double allMaxI = -1;
    int p = 0;
    int q = 0;

    if (verbose) Rcpp::Rcout << "Running STEPS for " << lambda.size() << " lambdas from "
			     << lambda[lambda.size()-1] << " to " << lambda[0] << "..."
			     << std::endl;

    for (Node n : d.getVariables()) {
        if (n.isDiscrete()) {
            q++;
        } else if (n.isContinuous()) {
	    p++;
	} else {
	    throw std::runtime_error("Invalid variable type for node " + n.getName());
	}
    }

    // go until we break by having instability better than threshold
    while(true) {
        if (verbose) Rcpp::Rcout << "  Testing lambda = " << lambda[currIndex] << std::endl;

        std::vector<double> lambdaCurr = { lambda[currIndex], lambda[currIndex], lambda[currIndex] };

        arma::mat adjMat;

        if (leaveOneOut) {
            adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr, threads);
        } else {
            adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr, threads, N, b);
        }

        // Rcpp::Rcout << "adjMat = " << adjMat << std::endl;

        double ccDestable = 0;
        double cdDestable = 0;
        double ddDestable = 0;
        ////////////////// TODO Decide if this is too harsh
        double numCC = 0;
        double numCD = 0;
        double numDD = 0;

        // We assume here that the subsamples have the variables in the same order
        for (int j = 0; j < d.getNumColumns(); j++) {
            for (int k = j + 1; k < d.getNumColumns(); k++) {
                Node one = d.getVariable(j);
                Node two = d.getVariable(k);

                if (one.isDiscrete() && two.isDiscrete()) {
                    numDD++;
                    ddDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                } else if (one.isDiscrete() || two.isDiscrete()) {
                    numCD++;
                    cdDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                } else {
                    numCC++;
                    ccDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                }
            }
        }

        double allDestable = ccDestable + cdDestable + ddDestable;
        allDestable = allDestable / (numCC + numCD + numDD);
	
        ccDestable = ccDestable / numCC;
        cdDestable = cdDestable / numCD;
        ddDestable = ddDestable / numDD;

	if (numCC == 0)
	    ccDestable = 0.5;
	if (numCD == 0)
	    cdDestable = 0.5;
	if (numDD == 0)
	    ddDestable = 0.5;

	if (verbose) Rcpp::Rcout << "    Instabilities for lambda = " << lambda[currIndex]
				 << ":  {" << ccDestable << ", " << cdDestable << ", "
				 << ddDestable << "}" << std::endl;

	if (allDestable >= gamma && oneLamb == -1 && currIndex > 0) {
	    oneLamb = lambda[currIndex-1];
        }
	
        if (allDestable >= allMax) {
            allMax = allDestable;
            allMaxI = lambda[currIndex];
        }

        if (ccDestable >= gamma && CC == -1 && currIndex > 0)
            CC = lambda[currIndex-1];
        if (cdDestable >= gamma && CD == -1 && currIndex > 0)
            CD = lambda[currIndex-1];
        if (ddDestable >= gamma && DD == -1 && currIndex > 0)
            DD = lambda[currIndex-1];
	
        if (ccDestable >= CCMax) {
            CCMax = ccDestable;
            CCMaxI = lambda[currIndex];
        }
        if (cdDestable >= CDMax) {
            CDMax = cdDestable;
            CDMaxI = lambda[currIndex];
        }
        if (ddDestable >= DDMax) {
            DDMax = ddDestable;
            DDMaxI = lambda[currIndex];
        }
        if (CC != -1 && CD != -1 && DD != -1 && oneLamb != -1)
            break;
        if (currIndex == lambda.size() - 1)
            break;
        currIndex++;
    }
    // EXIT_LOOP:

    if (CC == -1)
        CC = CCMaxI;
    if (CD == -1)
        CD = CDMaxI;
    if (DD == -1)
        DD = DDMaxI;

    std::vector<double> lambda = { CC, CD, DD };
    if (verbose) Rcpp::Rcout << "Selected lambdas: { " << CC << ", " << CD << ", " << DD
			     << " }" << std::endl;
    if (oneLamb == -1)
        origLambda = allMaxI;
    else
        origLambda = oneLamb;

    MGM m(d, lambda);
    m.learnEdges(iterLimit);

    if (computeStabs) {
        if (leaveOneOut) {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads);
        } else {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads, N, b);
        }
    }
    lastLambda = lambda;

    EdgeListGraph g = m.graphFromMGM();
    g.setHyperParam("lambda", arma::vec(lambda));
    
    return g;
}


EdgeListGraph STEPS::runStepsPath() {

    // Sort in descending order
    std::sort(lambda.begin(), lambda.end(), std::greater<double>());

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
    double CC = -1;
    double CD = -1;
    double DD = -1;
    double CCMax = 0;
    double CCMaxI = -1;
    double CDMax = 0;
    double CDMaxI = -1;
    double DDMax = 0;
    double DDMaxI = -1;
    double oneLamb = -1;
    double allMax = 0;
    double allMaxI = -1;
    int p = 0;
    int q = 0;

    if (verbose) Rcpp::Rcout << "Running STEPS for " << lambda.size() << " lambdas from "
			     << lambda[lambda.size()-1] << " to " << lambda[0] << "..."
			     << std::endl;

    for (Node n : d.getVariables()) {
        if (n.isDiscrete()) {
            q++;
        } else if (n.isContinuous()) {
	    p++;
	} else {
	    throw std::runtime_error("Invalid variable type for node " + n.getName());
	}
    }

    arma::umat samps;
    
    if (leaveOneOut) {
      samps = StabilityUtils::subSampleLOO(d);
    } else {
      samps = StabilityUtils::subSampleNoReplacement(d, b, N);
    }

    // Rcpp::Rcout << "Sample indices selected\n";
    
    // BlockingQueue<int> taskQueue(samps.n_rows);
    // Rcpp::Rcout << "taskQue constructed\n";
    int numVars = d.getNumColumns();
    std::vector<MGM> mgmList;
    // Rcpp::Rcout << "mgmList constructed\n";
    arma::mat thetaMat(numVars, numVars, arma::fill::zeros);
    // Rcpp::Rcout << "thetaMat constructed\n";
    std::mutex matMutex; // For protecting thetaMat
    // Rcpp::Rcout << "matMutex constructed\n";
    // Rcpp::Rcout << "Number of threads: " << parallelism << std::endl;
    RcppThread::ThreadPool pool(parallelism);
    // Rcpp::Rcout << "pool constructed\n";
    std::vector<double> lambdaCurr;

    // auto mgmPath = [&](int sample) {
		       
    // 		       // if (RcppThread::isInterrupted()) {
    // 		       // 	   return;
    // 		       // }

    // 		       RcppThread::Rcout << "Running sample " << sample << std::endl;
		       
    // 		       DataSet dataSubSamp(d, samps.row(sample));
    // 		       MGM mgm(dataSubSamp);
    // 		       std::vector<EdgeListGraph> mgmGraphs = mgm.searchPath(lambda);
		       
    // 		       for (int j = 0; j < lambda.size(); j++) {
    // 			   arma::mat curAdj = StabilityUtils::skeletonToMatrix(mgmGraphs[j], dataSubSamp);
    // 			   std::lock_guard<std::mutex> matLock(matMutex);
    // 			   theta.slice(j) += curAdj;
    // 		       }

    // 		       return;
    // 		   };

    // RcppThread::ThreadPool pool(threads);

    // for (int i = 0; i < samps.n_rows; i++) {
    // 	pool.push(mgmPath, i);
    // }

    // pool.join();

    // theta /= samps.n_rows;

    // Rcpp::Rcout << "About to fill mgmList\n";

    for (int i = 0; i < samps.n_rows; i++) {
    	DataSet dataSubSamp(d, samps.row(i));
    	mgmList.push_back(MGM(dataSubSamp));
    }

    // Rcpp::Rcout << "mgmList Filled\n";

    // auto producer = [&]() {
    // 			for (int i = 0; i < samps.n_rows; i++) {
    // 			    taskQueue.push(i);

    // 			    if (RcppThread::isInterrupted()) {
    // 				break;
    // 			    }
    // 			}

    // 			// Poison pills
    // 			for (int i = 0; i < threads; i++) {
    // 			    taskQueue.push(-1);
    // 			}
    // 		    };

    // auto consumer = [&]() {
    // 			while (true) {
    // 			    int sample = taskQueue.pop();

    // 			    if (sample < 0) break; // Poison pill

    // 			    if (RcppThread::isInterrupted()) {
    // 				break;
    // 			    }

    // 			    mgmList[sample].setLambda(lambdaCurr);
    // 			    EdgeListGraph g = mgmList[sample].search();
    // 			    arma::mat curAdj = StabilityUtils::skeletonToMatrix(g, d);

    // 			    {
    // 				std::lock_guard<std::mutex> matLock(matMutex);
    // 				thetaMat += curAdj;
    // 			    }
    // 			}
    // 		    };

    // Rcpp::Rcout << "producer/consumer defined\n";

    // std::vector<RcppThread::Thread> threadList;

    // for (int i = 0; i < threads; i++) {
    // 	threadList.push_back(RcppThread::Thread( mgmPath, i ));
    // }

    // for (int i = 0; i < threads; i++) {
    //     threadList[i].join();
    // }

    // thetaMat /= samps.n_rows;
    

    // go until we break by having instability better than threshold
    for (currIndex = 0; currIndex < lambda.size(); currIndex++) {
	
        if (verbose) Rcpp::Rcout << "  Testing lambda = " << lambda[currIndex] << std::endl;

        lambdaCurr = { lambda[currIndex], lambda[currIndex], lambda[currIndex] };

	thetaMat.fill(0);

	pool.parallelForEach(mgmList,
			     [&] (MGM& m) {
				 m.setLambda(lambdaCurr);
				 EdgeListGraph g = m.search();
				 arma::mat curAdj = StabilityUtils::skeletonToMatrix(g, d);
				 
				 {
				     std::lock_guard<std::mutex> matLock(matMutex);
				     thetaMat += curAdj;
				 }
			     });

	pool.wait();

	// auto producer = [&]() {
	// 		    for (int i = 0; i < samps.n_rows; i++) {
	// 			taskQueue.push(i);

	// 			if (RcppThread::isInterrupted()) {
	// 			    break;
	// 			}
	// 		    }

	// 		    // Poison pills
	// 		    for (int i = 0; i < parallelism; i++) {
	// 			taskQueue.push(-1);
	// 		    }
	// 		};

	// auto consumer = [&]() {
	// 		    while (true) {
	// 			int sample = taskQueue.pop();

	// 			if (sample < 0) break; // Poison pill

	// 			if (RcppThread::isInterrupted()) {
	// 			    break;
	// 			}

	// 			mgmList[sample].setLambda(lambdaCurr);
	// 			EdgeListGraph g = mgmList[sample].search();
	// 			arma::mat curAdj = StabilityUtils::skeletonToMatrix(g, d);

	// 			{
	// 			    std::lock_guard<std::mutex> matLock(matMutex);
	// 			    thetaMat += curAdj;
	// 			}
	// 		    }
	// 		};
	
	// pool.push(producer);

	// for (int i = 0; i < parallelism; i++) {
	//     pool.push(consumer);
	// }

	// pool.join();

        arma::mat adjMat = thetaMat/samps.n_rows;

	// arma::mat adjMat = theta.slice(currIndex);
	
        // if (leaveOneOut) {
        //     adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr, threads);
        // } else {
        //     adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr, threads, N, b);
        // }

        // Rcpp::Rcout << "adjMat = " << adjMat << std::endl;

        double ccDestable = 0;
        double cdDestable = 0;
        double ddDestable = 0;
        ////////////////// TODO Decide if this is too harsh
        double numCC = 0;
        double numCD = 0;
        double numDD = 0;

        // We assume here that the subsamples have the variables in the same order
        for (int j = 0; j < d.getNumColumns(); j++) {
            for (int k = j + 1; k < d.getNumColumns(); k++) {
                Node one = d.getVariable(j);
                Node two = d.getVariable(k);

                if (one.isDiscrete() && two.isDiscrete()) {
                    numDD++;
                    ddDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                } else if (one.isDiscrete() || two.isDiscrete()) {
                    numCD++;
                    cdDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                } else {
                    numCC++;
                    ccDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                }
            }
        }

        double allDestable = ccDestable + cdDestable + ddDestable;
        allDestable = allDestable / (numCC + numCD + numDD);
	
        ccDestable = ccDestable / numCC;
        cdDestable = cdDestable / numCD;
        ddDestable = ddDestable / numDD;

	if (numCC == 0)
	    ccDestable = 0.5;
	if (numCD == 0)
	    cdDestable = 0.5;
	if (numDD == 0)
	    ddDestable = 0.5;

	if (verbose) Rcpp::Rcout << "    Instabilities for lambda = " << lambda[currIndex]
				 << ":  {" << ccDestable << ", " << cdDestable << ", "
				 << ddDestable << "}" << std::endl;

	if (allDestable >= gamma && oneLamb == -1 && currIndex > 0) {
	    oneLamb = lambda[currIndex-1];
        }
	
        if (allDestable >= allMax) {
            allMax = allDestable;
            allMaxI = lambda[currIndex];
        }

        if (ccDestable >= gamma && CC == -1 && currIndex > 0)
            CC = lambda[currIndex-1];
        if (cdDestable >= gamma && CD == -1 && currIndex > 0)
            CD = lambda[currIndex-1];
        if (ddDestable >= gamma && DD == -1 && currIndex > 0)
            DD = lambda[currIndex-1];
	
        if (ccDestable >= CCMax) {
            CCMax = ccDestable;
            CCMaxI = lambda[currIndex];
        }
        if (cdDestable >= CDMax) {
            CDMax = cdDestable;
            CDMaxI = lambda[currIndex];
        }
        if (ddDestable >= DDMax) {
            DDMax = ddDestable;
            DDMaxI = lambda[currIndex];
        }
        if (CC != -1 && CD != -1 && DD != -1 && oneLamb != -1)
            break;
        // if (currIndex == lambda.size() - 1)
        //     break;
        // currIndex++;
    }
    // EXIT_LOOP:

    if (CC == -1)
        CC = CCMaxI;
    if (CD == -1)
        CD = CDMaxI;
    if (DD == -1)
        DD = DDMaxI;

    pool.join();

    std::vector<double> lambda = { CC, CD, DD };
    if (verbose) Rcpp::Rcout << "Selected lambdas: { " << CC << ", " << CD << ", " << DD
			     << " }" << std::endl;
    if (oneLamb == -1)
        origLambda = allMaxI;
    else
        origLambda = oneLamb;

    MGM m(d, lambda);
    m.learnEdges(iterLimit);

    if (computeStabs) {
        if (leaveOneOut) {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads);
        } else {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads, N, b);
        }
    }
    lastLambda = lambda;
    
    EdgeListGraph g = m.graphFromMGM();
    g.setHyperParam("lambda", arma::vec(lambda));
    
    return g;
}
