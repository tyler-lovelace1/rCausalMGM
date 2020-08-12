#ifndef INDEPENDENCETEST_HPP_
#define INDEPENDENCETEST_HPP_

#include "DataSet.hpp"
#include "Variable.hpp"
#include "armaLapack.hpp"
#include <list>
#include <string>

/**
 * Interface implemented by classes that do conditional independence testing. These classes are capable of serving as
 * conditional independence "oracles" for constraint-based searches.
 *
 * @author Don Crimbchin (djc2@andrew.cmu.edu)
 * @author Joseph Ramsey
 * @author Jack Fiore conversion to C++ 1/22
 */

class IndependenceTest
{
public:
    /**
     * @return an Independence test for a subset of the variables.
     */
    //   virtual IndependenceTest indTestSubset(std::vector<Variable*>& vars) = 0;

    /**
     * @return true if the given independence question is judged true, false if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are variables in the list returned by
     * getVariableNames().
     */
    virtual bool isIndependent(Variable* x, Variable* y, std::vector<Variable*>& z) = 0;

    /**
     * @return true if the given independence question is judged false, true if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are variables in the list returned by
     * getVariableNames().
     */
    virtual bool isDependent(Variable* x, Variable* y, std::vector<Variable*>& z) = 0;

    /**
     * @return the probability associated with the most recently executed independence test, of Double.NaN if p value is
     * not meaningful for tis test.
     */
    virtual double getPValue() = 0;

    /**
     * @return the list of variables over which this independence checker is capable of determinining independence
     * relations.
     */
    virtual std::vector<Variable*> getVariables() = 0;

    /**
     * @return the variable by the given name.
     */
    virtual Variable* getVariable(std::string name) = 0;

    /**
     * @return the list of names for the variables in getNodesInEvidence.
     */
    virtual std::vector<std::string> getVariableNames() = 0;

    /**
     * @return true if y is determined the variable in z.
     */
    virtual bool determines(std::vector<Variable*>& z, Variable* y) = 0;

    /**
     * @return the significance level of the independence test.
     * @throws UnsupportedOperationException if there is no significance level.
     */
    virtual double getAlpha() = 0;

    /**
     * Sets the significance level.
     */
    virtual void setAlpha(double alpha) = 0;

    virtual std::vector<DataSet*> getDataSets() = 0;

    virtual int getSampleSize() = 0;

    /**
     * A score that is higher with more likely models.
     */
    virtual double getScore() = 0;

};

#endif /* INDEPENDENCETEST_HPP_ */
