#ifndef REGRESSION_HPP_
#define REGRESSION_HPP_

#include "armaLapack.hpp"
#include "Node.hpp"
#include <list>

class Regression
{
  /**
   * Sets the significance level at which coefficients are judged to be
   * significant.
   *
   * @param alpha the significance level.
   */
  virtual void setAlpha(double alpha) = 0;


  /**
   * Regresses <code>target</code> on the <code>regressors</code>, yielding
   * a regression plane.
   *
   * @param target     the target variable, being regressed.
   * @param regressors the list of variables being regressed on.
   * @return the regression plane.
   */
  virtual RegressionResult regress(const Node& target, std::list<Node>& regressors) = 0;

};

#endif /* REGRESSION_HPP_ */
