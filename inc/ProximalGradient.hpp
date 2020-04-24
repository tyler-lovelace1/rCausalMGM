#ifndef PROXIMAL_GRADIENT_HPP_
#define PROXIMAL_GRADIENT_HPP_

#include <RcppArmadillo.h>
#include <ConvexProximal.hpp>
#include <chrono>
//TODO
#include "MGM.hpp" 

class ProximalGradient {

private:

    // Factors to alter Lipshitz constant estimate L, used for stepsize t = 1/L
    double beta;
    double alpha; //factor to decrease L otherwise

    bool edgeConverge; //if this is true we look to stop optimization when the edge predictions stop changing
    int noEdgeChangeTol = 3; //number of iterations in a row with no edge changes before we break

    int printIter = 100;
    double backtrackTol = 1e-10;
    

public:

    int iterComplete = 0;
    double timePerIter = 0;

    /**
     * Constructor using defaults from Becker et al 2011. beta = .5, alpha = .9
     */
    ProximalGradient() {
        beta = .5;
        alpha = .9;
        edgeConverge = false;
    }

    /**
     * Constructor, set parameters for a proximal gradient run
     *
     * @param beta (0,1) factor to increase L when Lipshitz violated, L = L_old/beta
     * @param alpha (0,1) factor to decrease L otherwise, L = L_old*alpha
     * @param edgeConverge
     */
    ProximalGradient(double beta, double alpha, bool edgeConverge);

    /**
     *  Positive edge change tolerance is the number of iterations with 0 edge changes needed to converge.
     *  Negative edge change tolerance means convergence happens when number of difference edges <= |edge change tol|.
     *  Default is 3.
     */
    void setEdgeChangeTol(int t){ noEdgeChangeTol = t; }

    arma::vec learnBackTrack(ConvexProximal *cp, arma::vec& Xin, double epsilon, int iterLimit, long time);
    arma::vec learnBackTrack(ConvexProximal *cp, arma::vec& Xin, double epsilon, int iterLimit);

};


#endif /* PROXIMAL_GRADIENT_HPP_ */
