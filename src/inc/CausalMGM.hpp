#ifndef CAUSALMGM_HPP_
#define CAUSALMGM_HPP_

#include "armaLapack.hpp"

#include "Node.hpp"
#include "ConvexProximal.hpp"
#include "DataSet.hpp"
#include "CausalMGMParams.hpp"
#include "ProximalGradient.hpp"
#include "EdgeListGraph.hpp"
#include <math.h>
#include <chrono>

/**
 * Implementation of Lee and Hastie's (2012) pseudolikelihood method
 * for learning Mixed Gaussian-Categorical Graphical Models, modified
 * to parameterize directed causal graphical models
 * Created by Tyler Lovelace on 12/16/22.
 */
class CausalMGM : public ConvexProximal {

private:

    //Continuous Data
    arma::mat xDat;

    //Discrete Data coded as integers
    arma::mat yDat;

    //Discrete Data coded as dummy variables
    arma::mat dDat;

    EdgeListGraph graph;

    std::vector<Node> variables;
    std::vector<Node> initVariables;

    std::map<Node, int> nodeIdxMap;

    Node dummyVar;

    arma::vec lambda = {0.5, 0.5, 0.5};

    long elapsedTime = 0;

    //Levels of Discrete variables
    std::vector<int> l;
    int lsum;
    std::vector<int> lcumsum;
    int p;
    int pDummy = 0;
    int q;
    int qDummy = 0;
    int n;
    long timeout = -1;

    int nParams = -1;

    CausalMGMParams params;
    CausalMGMParams active;

    //parameter weights
    arma::vec weights;

    double logsumexp(const arma::vec& x);

    void initParameters();  // init all parameters to zeros except for betad which is set to 1s
    void calcWeights();     // calculate parameter weights as in Lee and Hastie 
    void makeDummy();       // convert discrete data (in yDat) to a matrix of dummy variables (stored in dDat)
    void fixData();         // checks if yDat is zero indexed and converts to 1 index. zscores x

    void initActiveSet();   // init a binary mask for the set of active parameters

    // friend class Tests;

    

public:

    double timePerIter = 0;
    int iterCount = 0;

    CausalMGM() {}
    CausalMGM(EdgeListGraph& graph, arma::mat& x, arma::mat& y, std::vector<Node>& variables, std::vector<int>& l);
    CausalMGM(EdgeListGraph&& graph, arma::mat&& x, arma::mat&& y, std::vector<Node> variables, std::vector<int> l);
    CausalMGM(DataSet ds, EdgeListGraph& graph);
    // CausalMGM(DataSet& ds, std::vector<double>& lambda);

    // CausalMGM(const CausalMGM& other) = default;
    // CausalMGM& operator=(const CausalMGM& other) = default;
    // CausalMGM(const CausalMGM&& other) = default;
    // CausalMGM& operator=(const CausalMGM&& other) = default;
    // ~CausalMGM() = default;

    CausalMGMParams getParams() { return params; }
    void setParams(CausalMGMParams newParams) { params = newParams; }
    void setTimeout(long time) { timeout = time; }
    long getElapsedTime() { return elapsedTime; }
    arma::vec getLambda() { return lambda; }
    int getNParams() const { return nParams; }

    // void setVerbose(bool v) { verbose = v; }
    void setLambda(std::vector<double> lambda) { this->lambda = arma::vec(lambda); }

    // double calcLambdaMax();

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

    void iterUpdate(arma::vec& X) {}

    // /**
    //  *  Learn CausalMGM traditional way with objective function tolerance. Recommended for inference applications that need
    //  *  accurate pseudolikelihood
    //  *
    //  * @param epsilon tolerance in change of objective function
    //  * @param iterLimit iteration limit
    //  */
    void learn(double epsilon, int iterLimit);

    // /**
    //  *  Learn CausalMGM using edge convergence using default 3 iterations of no edge changes. Recommended when we only care about
    //  *  edge existence.
    //  *
    //  * @param iterLimit
    //  */
    // void learnEdges(int iterlimit);

    // /**
    //  *  Learn CausalMGM using edge convergence using edgeChangeTol (see ProximalGradient for documentation). Recommended when we only care about
    //  *  edge existence.
    //  *
    //  * @param iterLimit
    //  * @param edgeChangeTol
    //  */
    // void learnEdges(int iterlimit, int edgeChangeTol);

    // /**
    //  * Converts CausalMGM object to EdgeListGraph object with edges if edge parameters are non-zero. Loses all edge param information
    //  *
    //  * @return
    //  */
    // EdgeListGraph graphFromCausalMGM();

    // /**
    //  * Converts CausalMGM to matrix of doubles. uses 2-norm to combine c-d edge parameters into single value and f-norm for
    //  * d-d edge parameters.
    //  *
    //  * @return
    //  */
    // arma::mat adjMatFromCausalMGM();

    /**
     * Simple search command for GraphSearch implementation. Uses default edge convergence, 1000 iter limit.
     *
     * @return
     */
    CausalMGMParams search();

    // std::vector<EdgeListGraph> searchPath(std::vector<double> lambdas,
    // 					  arma::vec& loglik,
    // 					  arma::vec& nParams);

    // std::vector<EdgeListGraph> searchPathCV(std::vector<double> lambdas,
    // 					    int nfolds,
    // 					    arma::uvec& foldid,
    // 					    arma::mat& loglik);
    
    //  std::vector<EdgeListGraph> searchPath(std::vector<double> lambdas);
    
    friend void CausalMGMTest(const Rcpp::DataFrame &df, const int maxDiscrete);

};

#endif /* CAUSALMGM_HPP_ */
