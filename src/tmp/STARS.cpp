#include "STARS.hpp"

EdgeListGraph STARS::runStarsPar() {

    // Sort in order of sparsest to densest graphs
    if (alg == "mgm") {
	// Sort in descending order
	std::sort(params.begin(), params.end(), std::greater<double>());
    } else {
	// Sort in ascending order
	std::sort(params.begin(), params.end(), std::less<double>());
    }
    
    int currIndex = 0;
    double onePar = -1;
    double allMax = 0;
    double allMaxI = -1;
    int p = 0;
    int q = 0;
    int r = 0;

    Rcpp::Rcout << "Starting STARS" << std::endl;

    for (Variable* n : d.getVariables()) {
        if (n->isDiscrete()) {
            q++;
        } else if (n->isContinuous()) {
            p++;
        } else {
	    r++;
	}
    }

    Rcpp::Rcout << "Continuous Variables: " << p << std::endl;
    Rcpp::Rcout << "Discrete Variables: " << q << std::endl;
    Rcpp::Rcout << "Censored Variables: " << r << std::endl;

    // go until we break by having instability better than threshold
    while(true) {
        if (verbose) Rcpp::Rcout << "Testing param = " << params[currIndex] << std::endl;

	double destable = 0.0;
	
        if (leaveOneOut) {
            destable = StabilityUtils::stabilitySearchStars(d, alg, params[currIndex], initialGraph,
							    threads, adjacency);
        } else {
            destable = StabilityUtils::stabilitySearchStars(d, alg, params[currIndex], initialGraph,
							    threads, adjacency, N, b);
        }

	if (verbose) {
	    Rcpp::Rcout << "  Instability for param = " << params[currIndex]
			<< ":  " << destable << std::endl;
	}

	if (destable >= gamma && onePar == -1 && currIndex > 0) {
	    onePar = params[currIndex-1];
        }
	
        if (destable >= allMax) {
            allMax = destable;
            allMaxI = params[currIndex];
        }

	if (onePar != -1)
	    break;
		
        if (currIndex == params.size() - 1)
            break;
	
        currIndex++;
    }
    // EXIT_LOOP:

    if (onePar == -1)
	onePar = allMaxI;

    lastParam = onePar;

    
    std::vector<double> lambda;
    double alpha;
    if (alg == "mgm") {
	Rcpp::Rcout << "Lambda: { " << onePar << " }" << std::endl;
	
	if (r == 0) {
	    lambda = { onePar, onePar, onePar };
	} else {
	    lambda = { onePar, onePar, onePar, onePar, onePar };
	}
    } else {
	alpha = onePar;
    }

    EdgeListGraph g;

    if (alg == "mgm") {
	    
	MGM mgm(d, lambda);
	g = mgm.search();
	    
    } else if (alg == "pc") {
	    
	IndTestMulti itm(d, alpha);
	    
	PcStable pcs((IndependenceTest*) &itm);
	if (threads > 0) pcs.setThreads(threads);
	    
	pcs.setVerbose(false);
	pcs.setInitialGraph(initialGraph);
	g = pcs.search();
	    
    } else if (alg == "cpc") {
	    
	IndTestMulti itm(d, alpha);
	    
	CpcStable cpc((IndependenceTest*) &itm);
	if (threads > 0) cpc.setThreads(threads);
	    
	cpc.setVerbose(false);
	cpc.setInitialGraph(initialGraph);
	g = cpc.search();
	    
    } else if (alg == "pcm") {
	    
	IndTestMulti itm(d, alpha);
	    
	PcMax pcm((IndependenceTest*) &itm);
	if (threads > 0) pcm.setThreads(threads);
	    
	pcm.setVerbose(false);
	pcm.setInitialGraph(initialGraph);
	g = pcm.search();
	    
    } else if (alg == "fci") {
	    
	IndTestMulti itm(d, alpha);
	    
	Fci fci((IndependenceTest*) &itm);
	if (threads > 0) fci.setThreads(threads);
	    
	fci.setVerbose(false);
	fci.setInitialGraph(initialGraph);
	g = fci.search();
	    
    } else if (alg == "cfci") {
	    
	IndTestMulti itm(d, alpha);
	    
	Cfci cfci((IndependenceTest*) &itm);
	if (threads > 0) cfci.setThreads(threads);
	    
	cfci.setVerbose(false);
	cfci.setInitialGraph(initialGraph);
	g = cfci.search();
	    
    } else if (alg == "fcim") {
	    
	IndTestMulti itm(d, alpha);
	    
	FciMax fcim((IndependenceTest*) &itm);
	if (threads > 0) fcim.setThreads(threads);
	    
	fcim.setVerbose(false);
	fcim.setInitialGraph(initialGraph);
	g = fcim.search();
	    
    }

    // if (computeStabs) {
    //     if (leaveOneOut) {
    //         stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads);
    //     } else {
    //         stabilities = StabilityUtils::stabilitySearchPar(d, lambda, threads, N, b);
    //     }
    // }
    
    return g;
}
