#ifndef SCORE_HPP_
#define SCORE_HPP_

#include "DataSet.hpp"
#include "Node.hpp"
#include "armaLapack.hpp"
#include <list>
#include <string>

/**
 * Interface implemented by classes that do model scoring. These classes are capable of serving as
 * conditional independence "oracles" for constraint-based searches.
 *
 * @author Don Crimbchin (djc2@andrew.cmu.edu)
 * @author Joseph Ramsey
 * @author Jack Fiore conversion to C++ 1/22
 */

class Score
{
public:
    /**
     * @return an Independence test for a subset of the variables.
     */
    //   virtual Score indTestSubset(std::vector<Node>& vars) = 0;

    /**
     * @return score of x | z, z = <z1,...,zn>,
     * form x | z, z = <z1,...,zn>, where x, z1,...,zn are variables in the list returned by
     * getVariableNames().
     * 
     */
    virtual double localScore(const Node& x, const std::vector<Node>& z) = 0;

    
    // /**
    //  * @return score of x | z, z = <z1,...,zn>,
    //  * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are variables in the list returned by
    //  * getVariableNames().
    //  * Optionally return the p-value into pReturn
    //  */
    // virtual arma::vec localScoreVec(const Node& x, const std::vector<Node>& y,
    // 				    const std::vector<Node>& z) = 0;

    /**
     * @return the list of variables over which this independence checker is capable of determinining independence
     * relations.
     */
    virtual std::vector<Node> getVariables() = 0;

    /**
     * @return the variable by the given name.
     */
    virtual Node getVariable(std::string name) = 0;

    /**
     * @return the list of names for the variables in getNodesInEvidence.
     */
    virtual std::vector<std::string> getVariableNames() = 0;

    /**
     * @return the significance level of the independence test.
     * @throws UnsupportedOperationException if there is no significance level.
     */
    virtual double getPenalty() = 0;

    /**
     * Sets the significance level.
     */
    virtual void setPenalty(double penalty) = 0;

    // virtual std::vector<DataSet*> getDataSets() = 0;

    // virtual DataSet getData() = 0;

    virtual int getSampleSize() = 0;

};

#endif /* SCORE_HPP_ */
