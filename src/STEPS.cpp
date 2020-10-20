#include "STEPS.hpp"

EdgeListGraph STEPS::runStepsPar() {
    std::sort(lambda.begin(), lambda.end());

    int currIndex = 0;
    double CC = -1;
    double CD = -1;
    double DD = -1;
    double CCMax = 100;
    double CCMaxI = -1;
    double CDMax = 100;
    double CDMaxI = -1;
    double DDMax = 100;
    double DDMaxI = -1;
    double oneLamb = -1;
    double allMax = 100;
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
        Rcpp::Rcout << "Lambda: " << lambda[currIndex] << std::endl;

        std::vector<double> lambdaCurr = { lambda[currIndex], lambda[currIndex], lambda[currIndex] };

        arma::mat adjMat;

        if (leaveOneOut) {
            adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr);
        } else {
            adjMat = StabilityUtils::stabilitySearchPar(d, lambdaCurr, N, b);
        }

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
        if (allDestable <= gamma && oneLamb == -1) {
            oneLamb = lambda[currIndex];
        }
        if (allDestable <= allMax) {
            allMax = allDestable;
            allMaxI = lambda[currIndex];
        }

        ccDestable = ccDestable / numCC;
        cdDestable = cdDestable / numCD;
        ddDestable = ddDestable / numDD;
        if (ccDestable <= gamma && CC == -1)
            CC = lambda[currIndex];
        if (cdDestable <= gamma && CD == -1)
            CD = lambda[currIndex];
        if (ddDestable <= gamma && DD == -1)
            DD = lambda[currIndex];
        if (ccDestable <= CCMax) {
            CCMax = ccDestable;
            CCMaxI = lambda[currIndex];
        }
        if (cdDestable <= CDMax) {
            CDMax = cdDestable;
            CDMaxI = lambda[currIndex];
        }
        if (ddDestable <= DDMax) {
            DDMax = ddDestable;
            DDMaxI = lambda[currIndex];
        }
        if (CC != -1 && CD != -1 && DD != -1 && oneLamb != -1)
            goto EXIT_LOOP;
        if (currIndex == lambda.size() - 1)
            goto EXIT_LOOP;
        currIndex++;
    }
    EXIT_LOOP:

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
        arma::mat stabs;
        if (leaveOneOut) {
            stabs = StabilityUtils::stabilitySearchPar(d, lambda);
        } else {
            stabs = StabilityUtils::stabilitySearchPar(d, lambda, N, b);
        }
    }
    lastLambda = lambda;
    return m.graphFromMGM();
}