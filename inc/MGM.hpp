#ifndef MGM_HPP_
#define MGM_HPP_

#include "ConvexProximal.hpp"
#include "DataSet.hpp"
#include "ContinuousVariable.hpp"
#include "DiscreteVariable.hpp"
#include "MGMParams.hpp"
#include <math.h>
#include <RcppArmadillo.h>

class MGM : public ConvexProximal {

private:

    //Continuous Data
    arma::mat xDat;

    //Discrete Data coded as integers
    arma::mat yDat;

    //Discrete Data coded as dummy variables
    arma::mat dDat;

    std::vector<Variable*> variables;
    std::vector<Variable*> initVariables;

    arma::vec lambda;

    long elapsedTime = 0;

    //Levels of Discrete variables
    std::vector<int> l;
    int lsum;
    std::vector<int> lcumsum;
    int p;
    int q;
    int n;
    long timeout = -1;

    //parameter weights
    arma::mat weights;

    MGMParams params;

    double logsumexp(const arma::vec& x);

    void initParameters();
    void calcWeights();
    void makeDummy();
    void fixData();

    static void runTests1();
    static void runTests2();

public:
    double timePerIter = 0;
    int iterCount = 0;

    MGM() {}
    MGM(arma::mat& x, arma::mat& y, std::vector<Variable*>& variables, std::vector<int>& l, std::vector<double>& lambda);
    MGM(DataSet& ds, std::vector<double>& lambda);
    ~MGM();

    void setParams(MGMParams& newParams) {params = newParams;}
    void setTimeout(long time) { timeout = time; }
    long getElapsedTime() { return elapsedTime; }

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

    void learn(double epsilon, int iterLimit);
    void learnEdges(int iterlimit);
    void learnEdges(int iterlimit, int edgeChangeTol);

    // TODO Graph
    // Graph graphFromMGM();

    arma::mat& adjMatFromMGM();

    // Graph search();

};

// margeSum() replaced with sum(M, dim)
// upperTri and lowerTri replaced with trimatu() / trimatl()

#endif /* MGM_HPP_ */
