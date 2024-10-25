#ifndef COXMGM_HPP_
#define COXMGM_HPP_

// [[Rcpp::depends(RcppThread)]]

#include "armaLapack.hpp"
#include "RcppThread.h"
#include "ConvexProximal.hpp"
#include "DataSet.hpp"
#include "CoxMGMParams.hpp"
#include "ProximalGradient.hpp"
#include "EdgeListGraph.hpp"
#include <math.h>
#include <chrono>

/**
 * Implementation of CoxMGM pseudolikelihood method for learning
 * Mixed Cox-Gaussian-Categorical Graphical Models
 * Created by Tyler Lovelace on 4/19/22.
 */
class CoxMGM : public ConvexProximal {

private:

    //Continuous Data
    arma::mat xDat;

    //Discrete Data coded as integers
    arma::mat yDat;

    //Discrete Data coded as dummy variables
    arma::mat dDat;

    //Cox Data from z = eta - grad_eta / hess_eta
    arma::mat cDat;

    //Cox Data from z = eta - grad_eta / hess_eta
    arma::mat zDat;

    //Cox Data from z = eta + grad_eta
    arma::mat wzDat;

    //Cox variables number of strata
    arma::uvec numStrata;

    //Cox variables strata index list
    std::vector<std::vector<arma::uvec>> idxList;
    
    //Cox variables failure time order
    // arma::umat orderMat;
    std::vector<std::vector<arma::uvec>> orderList;

    //Cox variables censoring indicator
    arma::umat censMat;
    std::vector<std::vector<arma::uvec>> censList;
    
    //Cox variables tied times counts
    std::vector<std::vector<arma::uvec>> HList;

    std::vector<Node> variables;
    std::vector<Node> initVariables;

    Node dummyVar;

    arma::vec lambda;

    long elapsedTime = 0;

    //Levels of Discrete variables
    std::vector<int> l;
    int lsum;
    std::vector<int> lcumsum;
    int p;
    int pDummy = 0;
    int q;
    int qDummy = 0;
    int r;
    int n;
    long timeout = -1;

    CoxMGMParams params;

    //parameter weights
    arma::vec weights;
    arma::mat weightMat;

    arma::vec nullCoxloss;
    arma::vec oldCoxloss;

    arma::mat resid;
    arma::mat catResid;

    //Cox weights
    arma::mat coxWeights;
    arma::mat coxgrad;
    arma::mat diagHess;
    arma::rowvec fitWeight;

    // bool verbose = false;

    double logsumexp(arma::vec x);

    arma::vec logsumexp(arma::mat x);

    void initParameters();  // init all parameters to zeros except for betad which is set to 1s
    void calcWeights();     // calculate parameter weights as in Lee and Hastie
    // void calcZWeights();    // calculate parameter weights for Cox features 
    void makeDummy();       // convert discrete data (in yDat) to a matrix of dummy variables (stored in dDat)
    void fixData();         // checks if yDat is zero indexed and converts to 1 index. zscores x
    arma::vec coxGradHess(arma::mat& eta, arma::mat& grad, arma::mat& diagHess);

    friend class Tests;

public:

    double timePerIter = 0;
    int iterCount = 0;

    CoxMGM() {}
    // CoxMGM(arma::mat& x, arma::mat& y, std::vector<Node>& variables, std::vector<int>& l, std::vector<double>& lambda);
    CoxMGM(DataSet& ds);
    CoxMGM(DataSet& ds, std::vector<double>& lambda);

    // CoxMGM(CoxMGM& other) = default;
    // CoxMGM& operator=(CoxMGM& other) = default;
    // CoxMGM(CoxMGM&& other) = default;
    // CoxMGM& operator=(CoxMGM&& other) = default;
    // ~CoxMGM() = default;


    CoxMGMParams getParams() {return params;}
    void setParams(CoxMGMParams& newParams) {params = newParams;}
    void setTimeout(long time) { timeout = time; }
    long getElapsedTime() { return elapsedTime; }

    std::string printParameters(arma::vec& X);

    // void setVerbose(bool v) { verbose = v; }
    void setLambda(std::vector<double> lambda) { this->lambda = arma::vec(lambda); }

    double calcLambdaMax();

    /**
     * Calculate value of g(X) and gradient of g(X) at the same time for efficiency reasons.
     *
     * @param X input Vector
     * @param Xout gradient of g(X)
     * @return value of g(X)
     */
    double smooth(arma::vec& parIn, arma::vec& gradOutVec);

    /**
     * Calculate value of smooth function g(X)
     *
     * @param X input vector
     * @return value of g(X)
     */
    double smoothValue(arma::vec& parIn);

    /**
     * Calculate value of h(X) and proxOperator of h(X) at the same time for efficiency reasons.
     *
     * @param t positive parameter for prox operator
     * @param X input vector
     * @param Xout vector solution to prox_t(X)
     * @return value of h(X)
     */
    double nonSmooth(double t, arma::vec& X, arma::vec& pX);

    /**
     * Calculate value of h(X)
     *
     * @param X input vector
     * @return value of h(X)
     */
    double nonSmoothValue(arma::vec& parIn);

    /**
     * Gradient of smooth function g(X)
     *
     * @param X input vector
     * @return vector containing gradient of g(X)
     */
    arma::vec smoothGradient(arma::vec& parIn);

    /**
     * A proximal operator is the solution to this optimization problem:
     *     prox_t(x) = argmin_z \frac{1}{2t} \|x-z\|^2_2 + h(x)
     *
     * @param t positive parameter for prox operator
     * @param X input vector
     * @return vector solution to prox_t(X)
     */
    arma::vec proximalOperator(double t, arma::vec& X);


    void iterUpdate(arma::vec& parIn);

    /**
     *  Learn CoxMGM traditional way with objective function tolerance. Recommended for inference applications that need
     *  accurate pseudolikelihood
     *
     * @param epsilon tolerance in change of objective function
     * @param iterLimit iteration limit
     */
    void learn(double epsilon, int iterLimit);

    /**
     *  Learn CoxMGM using edge convergence using default 3 iterations of no edge changes. Recommended when we only care about
     *  edge existence.
     *
     * @param iterLimit
     */
    void learnEdges(int iterlimit);

    /**
     *  Learn CoxMGM using edge convergence using edgeChangeTol (see ProximalGradient for documentation). Recommended when we only care about
     *  edge existence.
     *
     * @param iterLimit
     * @param edgeChangeTol
     */
    void learnEdges(int iterlimit, int edgeChangeTol);

    /**
     * Converts CoxMGM object to EdgeListGraph object with edges if edge parameters are non-zero. Loses all edge param information
     *
     * @return
     */
    EdgeListGraph graphFromCoxMGM();

    /**
     * Converts CoxMGM to matrix of doubles. uses 2-norm to combine c-d edge parameters into single value and f-norm for
     * d-d edge parameters.
     *
     * @return
     */
    arma::mat adjMatFromCoxMGM();

    /**
     * Simple search command for GraphSearch implementation. Uses default edge convergence, 1000 iter limit.
     *
     * @return
     */
    EdgeListGraph search();

    std::vector<EdgeListGraph> searchPath(std::vector<double> lambdas,
					  arma::vec& loglik,
					  arma::vec& nParams);

    std::vector<EdgeListGraph> searchPath(std::vector<double> lambdas);

    friend void CoxMGMTest(const Rcpp::DataFrame &df, const int maxDiscrete);

};

#endif /* COXMGM_HPP_ */
