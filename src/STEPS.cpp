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


std::vector<EdgeListGraph> STEPS::runStepsPath(arma::mat& instabs, arma::umat& samps) {

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
    double SC = -1;
    double SD = -1;
    double CCMax = 0;
    int CCMaxI = -1;
    double CDMax = 0;
    int CDMaxI = -1;
    double DDMax = 0;
    int DDMaxI = -1;
    double SCMax = 0;
    int SCMaxI = -1;
    double SDMax = 0;
    int SDMaxI = -1;
    double oneLamb = -1;
    double allMax = 0;
    int allMaxI = -1;
    int p = 0;
    int q = 0;
    int r = 0;

    if (verbose) Rcpp::Rcout << "Running STEPS for " << lambda.size() << " lambdas from "
			     << lambda[lambda.size()-1] << " to " << lambda[0] << "..."
			     << std::endl;

    for (Node n : d.getVariables()) {
        if (n.isDiscrete()) {
            q++;
        } else if (n.isContinuous()) {
	    p++;
	} else if (n.isCensored()) {
	    r++;
	} else {
	    throw std::runtime_error("Invalid variable type for node " + n.getName());
	}
    }

    // arma::umat samps;
    
    if (leaveOneOut) {
	samps = StabilityUtils::subSampleLOO(d);
    } else {
	samps = StabilityUtils::subSampleNoReplacement(d, b, N);
    }
    
    int numVars = d.getNumColumns();
    std::vector<MGM> mgmList;
    std::vector<CoxMGM> coxmgmList;
    arma::mat thetaMat(numVars, numVars, arma::fill::zeros);
    std::mutex matMutex; // For protecting thetaMat
    RcppThread::ThreadPool pool(std::max(1, std::min(parallelism, N)));
    std::vector<double> lambdaCurr;


    if (r > 0) {
	for (int i = 0; i < samps.n_rows; i++) {
	    DataSet dataSubSamp(d, samps.row(i));
	    coxmgmList.push_back(CoxMGM(dataSubSamp));
	}
    } else {
	for (int i = 0; i < samps.n_rows; i++) {
	    DataSet dataSubSamp(d, samps.row(i));
	    mgmList.push_back(MGM(dataSubSamp));
	}
    }    

    // go until we break by having instability better than threshold
    for (currIndex = 0; currIndex < lambda.size(); currIndex++) {
	
        if (verbose) Rcpp::Rcout << "  Testing lambda = " << lambda[currIndex] << std::endl;

        // lambdaCurr = { lambda[currIndex], lambda[currIndex], lambda[currIndex] };

	if (r == 0) {
	    lambdaCurr = { lambda[currIndex], lambda[currIndex], lambda[currIndex] };
	} else {
	    lambdaCurr = { lambda[currIndex], lambda[currIndex], lambda[currIndex],
			   lambda[currIndex], lambda[currIndex] };
	}

	thetaMat.fill(0);

	if (r > 0) {
	    pool.parallelForEach(coxmgmList,
				 [&] (CoxMGM& m) {
				     // RcppThread::Rcout << "MGM run started...\n";
				     // EdgeListGraph ig = m.graphFromCoxMGM();
				     // RcppThread::Rcout << "Graph starting from " << ig << "\n";
				     m.setLambda(lambdaCurr);
				     EdgeListGraph g = m.search();
				     arma::mat curAdj = StabilityUtils::skeletonToMatrix(g, d);
				 
				     {
					 std::lock_guard<std::mutex> matLock(matMutex);
					 thetaMat += curAdj;
				     }
				     // RcppThread::Rcout << "MGM run completed\n";
				 });
	} else {
	    pool.parallelForEach(mgmList,
				 [&] (MGM& m) {
				     // RcppThread::Rcout << "MGM run started...\n";
				     // EdgeListGraph ig = m.graphFromMGM();
				     // RcppThread::Rcout << "Graph starting from " << ig.getEdges().size() << " edges\n";
				     m.setLambda(lambdaCurr);
				     EdgeListGraph g = m.search();
				     arma::mat curAdj = StabilityUtils::skeletonToMatrix(g, d);
				 
				     {
					 std::lock_guard<std::mutex> matLock(matMutex);
					 thetaMat += curAdj;
				     }
				     // RcppThread::Rcout << "MGM run completed\n";
				 });
	}
	
	pool.wait();

        arma::mat adjMat = thetaMat/samps.n_rows;

        double ccDestable = 0;
        double cdDestable = 0;
        double ddDestable = 0;
	double scDestable = 0;
        double sdDestable = 0;
        ////////////////// TODO Decide if this is too harsh
        double numCC = 0;
        double numCD = 0;
        double numDD = 0;
	double numSC = 0;
        double numSD = 0;

        // We assume here that the subsamples have the variables in the same order
        for (int j = 0; j < d.getNumColumns(); j++) {
            for (int k = j + 1; k < d.getNumColumns(); k++) {
                Node one = d.getVariable(j);
                Node two = d.getVariable(k);

                if (one.isDiscrete() && two.isDiscrete()) {
                    numDD++;
                    ddDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                } else if (one.isDiscrete() || two.isDiscrete()) {
		    if (one.isCensored() || two.isCensored()) {
			numSD++;
			sdDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
		    } else {
			numCD++;
			cdDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
		    }
                } else {
		    if (one.isCensored() || two.isCensored()) {
			numSC++;
			scDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
		    } else {
			numCC++;
			ccDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
		    }
                }
            }
        }

        double allDestable = ccDestable + cdDestable + ddDestable + scDestable + sdDestable;
        allDestable = allDestable / (numCC + numCD + numDD + numSC + numSD);
	
        ccDestable = ccDestable / numCC;
        cdDestable = cdDestable / numCD;
        ddDestable = ddDestable / numDD;
	scDestable = scDestable / numSC;
	sdDestable = sdDestable / numSD;

	if (numCC == 0)
	    ccDestable = 0.5;
	if (numCD == 0)
	    cdDestable = 0.5;
	if (numDD == 0)
	    ddDestable = 0.5;
	if (numSC == 0)
	    scDestable = 0.5;	
	if (numSD == 0)
	    sdDestable = 0.5;

	instabs(currIndex, 0) = std::max(CCMax, ccDestable);
	instabs(currIndex, 1) = std::max(CDMax, cdDestable);
	instabs(currIndex, 2) = std::max(DDMax, ddDestable);
	if (r > 0) {
	    instabs(currIndex, 3) = std::max(SCMax, scDestable);
	    instabs(currIndex, 4) = std::max(SDMax, sdDestable);
	    instabs(currIndex, 5) = std::max(allMax, allDestable);
	} else {
	    instabs(currIndex, 3) = std::max(allMax, allDestable);
	}

	if (verbose) {
	    Rcpp::Rcout << "  Overall instability for lambda = " << lambda[currIndex]
			<< ":  " << allDestable << std::endl;
	    if (r > 0) {
		Rcpp::Rcout << "  Instabilities for lambda = " << lambda[currIndex]
			    << ":  {" << ccDestable << ", " << cdDestable << ", "
			    << ddDestable << ", " << scDestable << ", "
			    << sdDestable << "}" << std::endl;
	    } else {
		Rcpp::Rcout << "  Instabilities for lambda = " << lambda[currIndex]
			    << ":  {" << ccDestable << ", " << cdDestable << ", "
			    << ddDestable << "}" << std::endl;
	    }
	}

	if (allDestable >= gamma && oneLamb == -1 && currIndex > 0) {
	    oneLamb = lambda[currIndex-1];
        }
	
        if (allDestable >= allMax) {
            allMax = allDestable;
            allMaxI = currIndex;
        }

        if (ccDestable >= gamma && CC == -1 && currIndex > 0)
            CC = lambda[currIndex-1];
        if (cdDestable >= gamma && CD == -1 && currIndex > 0)
            CD = lambda[currIndex-1];
        if (ddDestable >= gamma && DD == -1 && currIndex > 0)
            DD = lambda[currIndex-1];
	if (scDestable >= gamma && SC == -1 && currIndex > 0)
            SC = lambda[currIndex-1];
        if (sdDestable >= gamma && SD == -1 && currIndex > 0)
            SD = lambda[currIndex-1];
	
        if (ccDestable >= CCMax) {
            CCMax = ccDestable;
            CCMaxI = currIndex;
        }
        if (cdDestable >= CDMax) {
            CDMax = cdDestable;
            CDMaxI = currIndex;
        }
        if (ddDestable >= DDMax) {
            DDMax = ddDestable;
            DDMaxI = currIndex;
        }
	if (scDestable >= SCMax) {
            SCMax = scDestable;
            SCMaxI = currIndex;
        }
        if (sdDestable >= SDMax) {
            SDMax = sdDestable;
            SDMaxI = currIndex;
        }

	if (r > 0) {
	    if (CC != -1 && CD != -1 && DD != -1 && SC != -1 && SD != -1 && oneLamb != -1)
	     	break;
	} else {
	    if (CC != -1 && CD != -1 && DD != -1 && oneLamb != -1)
		break;
	}
        // if (currIndex == lambda.size() - 1)
        //     break;
        // currIndex++;
    }
    // EXIT_LOOP:

    pool.join();

    if (CC == -1)
        CC = lambda[CCMaxI];
    if (CD == -1)
        CD = lambda[CDMaxI];
    if (DD == -1)
        DD = lambda[DDMaxI];
    if (SC == -1)
        SC = lambda[SCMaxI];
    if (SD == -1)
        SD = lambda[SDMaxI];
    if (oneLamb == -1)
        oneLamb = lambda[allMaxI];

    // std::vector<double> lambda = { CC, CD, DD };
    // if (verbose) Rcpp::Rcout << "Selected lambdas: { " << CC << ", " << CD << ", " << DD
    // 			     << " }" << std::endl;

    // for (arma::uword i = CCMaxI; i >= 0; i--) {
    // 	if (instabs(i,0) < gamma || i == 0) {
    // 	    CC = lambda[i];
    // 	    break;
    // 	}
    // }

    // for (arma::uword i = CDMaxI; i >= 0; i--) {
    // 	if (instabs(i,1) < gamma || i == 0) {
    // 	    CD = lambda[i];
    // 	    break;
    // 	}
    // }

    // for (arma::uword i = DDMaxI; i >= 0; i--) {
    // 	if (instabs(i,2) < gamma || i == 0) {
    // 	    DD = lambda[i];
    // 	    break;
    // 	}
    // }

    // for (arma::uword i = allMaxI; i >= 0; i--) {
    // 	if (instabs(i,3) < gamma || i == 0) {
    // 	    oneLamb = lambda[i];
    // 	    break;
    // 	}
    // }

    // if (CC == -1)
    //     CC = lambda[0];
    // if (CD == -1)
    //     CD = lambda[0];
    // if (DD == -1)
    //     DD = lambda[0];
    // if (oneLamb == -1)
    //     oneLamb = lambda[0];

    std::vector<double> lambda;
    if (r > 0) {
	lambda = { CC, CD, DD, SC, SD };
	if (verbose) {
	    Rcpp::Rcout << "Selected Lambdas: { " << CC << " " << CD << " " << DD
			<< " " << SC << " " << SD << " }" << std::endl;
	    Rcpp::Rcout << "StARS lambda: " << oneLamb << std::endl;
	}
    } else {
	lambda = { CC, CD, DD };
	if (verbose) {
	    Rcpp::Rcout << "Selected lambdas: { " << CC << ", " << CD << ", " << DD
			<< " }" << std::endl;
	    Rcpp::Rcout << "StARS lambda: " << oneLamb << std::endl;
	}
    }
    
    if (oneLamb == -1)
        origLambda = lambda[allMaxI];
    else
        origLambda = oneLamb;
    
    EdgeListGraph gSteps, gStars;
    if (d.isCensored()) {
	CoxMGM m(d, lambda);
	gSteps = m.search();
	
	lambda = { oneLamb, oneLamb, oneLamb, oneLamb, oneLamb };
	
	m = CoxMGM(d, lambda);
	gStars = m.search();
    } else {
	MGM m(d, lambda);
	gSteps = m.search();

	lambda = { oneLamb, oneLamb, oneLamb, };
	
	m = MGM(d, lambda);
	gStars = m.search();
    }

    if (computeStabs) {
        if (leaveOneOut) {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads);
        } else {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads, N, b);
        }
    }
    lastLambda = lambda;
    
    // EdgeListGraph g = m.graphFromMGM();
    // gSteps.setHyperParam("lambda", arma::vec(lambda));
    
    return { gSteps, gStars };
}
