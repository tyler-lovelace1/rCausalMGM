#ifndef INDTESTMULTI_HPP_
#define INDTESTMULTI_HPP_

#include "armaLapack.hpp"
#include "IndependenceTest.hpp"
#include "LogisticRegression.hpp"
#include "LinearRegression.hpp"
#include "Variable.hpp"
#include "DataSet.hpp"

class IndTestMulti : public IndependenceTest {
private:
    int timesCalled;
    DataSet originalData;
    std::vector<Variable*> searchVariables;
    DataSet internalData;
    double alpha;
    double lastP;
    std::map<Variable*, std::vector<Variable*>> variablesPerNode;
    LogisticRegression logisticRegression;
    LinearRegression regression;
    bool verbose = false;
    bool preferLinear;

public:
    IndTestMulti() {}

    IndTestMulti(DataSet& data, double alpha);

    IndTestMulti(DataSet& data, double alpha, bool preferLinear);

    int reset();

    bool isIndependent(Variable* x, Variable* y, std::vector<Variable*>& z);

    // private synchronized List<Node> expandVariable(DataSet dataSet, Node node)
    std::vector<Variable*> expandVariable(DataSet& dataSet, Variable* var);

    bool isIndependentMultinomialLogisticRegression(Variable* x, Variable* y, std::vector<Variable*>& z);

    // This takes an inordinate amount of time. -jdramsey 20150929
    arma::uvec getNonMissingRows(Variable* x, Variable* y, std::vector<Variable*>& z);

    bool isMissing(Variable* x, int i);

    double multiLL(arma::mat& coeffs, Variable* dep, std::vector<Variable*>& indep);

    bool isIndependentRegression(Variable* x, Variable* y, std::vector<Variable*>& z);

    /**
     * @return true if the given independence question is judged false, true if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are searchVariables in the list returned by
     * getVariableNames().
     */
    bool isDependent(Variable* x, Variable* y, std::vector<Variable*>& z) { return !this->isIndependent(x, y, z); }

    /**
     * @return the probability associated with the most recently executed independence test, of Double.NaN if p value is
     * not meaningful for tis test.
     */        //STUB ????????????
    double getPValue() { return this->lastP; }

    /**
     * @return the list of searchVariables over which this independence checker is capable of determinining independence
     * relations.
     */
    // Make sure the variables from the ORIGINAL data set are returned, not the modified dataset!
    std::vector<Variable*> getVariables() {return this->searchVariables; }

    /**
     * @return the list of variable varNames.
     */
    std::vector<std::string> getVariableNames();

    Variable* getVariable(std::string name);

    /**
     * @return true if y is determined the variable in z.
     */
    bool determines(std::vector<Variable*>& z, Variable* y) { return false; /*stub ?????*/ }

    double getAlpha() { return this->alpha; /*STUB ?????*/ }

    void setAlpha(double alpha) { this->alpha = alpha; }

    DataSet getData() { return this->originalData; }

    double getScore() { return this->getPValue(); }

    bool isVerbose() { return this->verbose; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    arma::mat getSubsetData(DataSet& origData, std::vector<Variable*>& varSubset);

    int getSampleSize() { return originalData.getNumRows(); }

    std::vector<DataSet*> getDataSets() { return { &originalData }; }

    friend void indTestMultiTest(const Rcpp::DataFrame& df);

};

#endif /* INDTESTMULTI_HPP_ */
