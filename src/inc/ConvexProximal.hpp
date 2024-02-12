/**
 * This interface should be used for non-differentiable convex functions that are decomposable such that
 * f(x) = g(x) + h(x) where g(x) is a differentiable convex function (i.e. smooth) and h(x) is a convex but not
 * necessarily differentiable (i.e. non-smooth) and has a proximal operator prox_t(x) = argmin_z 1/(2t) norm2(x-z)^2 +
 * h(z) has a solution for any t > 0. Typically g(x) will be a likelihood, and h(x) is a penalty term (as in l_1 in the
 * lasso)
 *
 *
 * Created by ajsedgewick on 8/4/15.
 * Converted to C++ by Max Dudek on 02/2020
 */
#ifndef CONVEXPROXIMAL_HPP_
#define CONVEXPROXIMAL_HPP_

#include "armaLapack.hpp"
#include <sstream>
#include <string>

class ConvexProximal {

protected:
    bool verbose = false;

public:
    // ConvexProximal() {}
    // ConvexProximal(const ConvexProximal& other) = default;
    // ConvexProximal& operator=(const ConvexProximal& other) = default;
    // ConvexProximal(ConvexProximal&& other) = default;
    // ConvexProximal& operator=(ConvexProximal&& other) = default;
    // virtual ~ConvexProximal() = 0;
  
    /**
     * Calculate value of smooth function g(X)
     *
     * @param X input vector
     * @return value of g(X)
     */
    virtual double smoothValue(arma::vec& X) = 0;

    /**
     * Calculate value of h(X)
     *
     * @param X input vector
     * @return value of h(X)
     */
    virtual double nonSmoothValue(arma::vec& X) = 0;

    /**
     * Gradient of smooth function g(X)
     *
     * @param X input vector
     * @return vector containing gradient of g(X)
     */
    virtual arma::vec smoothGradient(arma::vec& X) = 0;

    /**
     * A proximal operator is the solution to this optimization problem:
     *     prox_t(x) = argmin_z \frac{1}{2t} \|x-z\|^2_2 + h(x)
     *
     * @param t positive parameter for prox operator
     * @param X input vector
     * @return vector solution to prox_t(X)
     */
    virtual arma::vec proximalOperator(double t, arma::vec& X) = 0;

    /**
     * Calculate value of g(X) and gradient of g(X) at the same time for efficiency reasons.
     *
     * @param X input Vector
     * @param Xout gradient of g(X)
     * @return value of g(X)
     */
    virtual double smooth(arma::vec& X, arma::vec& Xout) {
        Xout = smoothGradient(X);
        return smoothValue(X);
    }

    /**
     * Calculate value of h(X) and proxOperator of h(X) at the same time for efficiency reasons.
     *
     * @param t positive parameter for prox operator
     * @param X input vector
     * @param Xout vector solution to prox_t(X)
     * @return value of h(X)
     */
    virtual double nonSmooth(double t, arma::vec& X, arma::vec& Xout) {
        Xout = proximalOperator(t,X);
        return nonSmoothValue(X);
    }

    virtual void iterUpdate(arma::vec& X) = 0;

    virtual std::string printParameters(arma::vec& X) = 0;

    void setVerbose(bool v) { verbose = v; }
    
    bool isVerbose() { return verbose; }

};

#endif /* CONVEXPROXIMAL_HPP_ */
