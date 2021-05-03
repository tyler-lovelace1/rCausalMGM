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

    for (Variable* n : d.getVariables()) {
        if (n->isDiscrete()) {
            q++;
        } else {
            p++;
        }
    }

    // go until we break by having instability better than threshold
    while(true) {
        if (verbose) Rcpp::Rcout << "Testing lambda = " << lambda[currIndex] << std::endl;

        std::vector<double> lambdaCurr = { lambda[currIndex], lambda[currIndex], lambda[currIndex] };

        arma::mat adjMat;

        if (leaveOneOut) {
            adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr);
        } else {
            adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr, N, b);
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
                Variable* one = d.getVariable(j);
                Variable* two = d.getVariable(k);

                if (one->isDiscrete() && two->isDiscrete()) {
                    numDD++;
                    ddDestable += 2 * adjMat(j, k) * (1 - adjMat(j, k));
                } else if (one->isDiscrete() || two->isDiscrete()) {
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

	if (verbose) Rcpp::Rcout << "  Instabilities for lambda = " << lambda[currIndex]
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
    Rcpp::Rcout << "Lambdas: { " << CC << " " << CD << " " << DD << " }" << std::endl;
    if (oneLamb == -1)
        origLambda = allMaxI;
    else
        origLambda = oneLamb;

    MGM m(d, lambda);
    m.learnEdges(iterLimit);

    if (computeStabs) {
        if (leaveOneOut) {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda);
        } else {
            stabilities = StabilityUtils::stabilitySearchPar(d, lambda, N, b);
        }
    }
    lastLambda = lambda;
    return m.graphFromMGM();
}
