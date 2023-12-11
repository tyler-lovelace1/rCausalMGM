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
    std::vector<Node> searchVariables;
    DataSet internalData;
    double alpha;
    double lastP;
    std::map<Node, std::vector<Node>> variablesPerNode;
    LogisticRegression logisticRegression;
    LinearRegression regression;
    bool verbose = false;
    bool preferLinear;

    // /**
    //  * Default constructor node will serve as a null node
    //  */
    // static const Node& nullNode;


public:
    IndTestMulti() {}

    IndTestMulti(DataSet& data, double alpha);

    IndTestMulti(DataSet& data, double alpha, bool preferLinear);

    IndTestMulti(const IndTestMulti& other) = default;
    IndTestMulti& operator=(const IndTestMulti& other) = default;

    IndTestMulti(IndTestMulti&& other) = default;
    IndTestMulti& operator=(IndTestMulti&& other) = default;

    ~IndTestMulti() = default;

    int reset();

    bool isIndependent(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn = NULL);

    // private synchronized List<Node> expandVariable(DataSet dataSet, const Node& node)
    std::vector<Node> expandVariable(DataSet& dataSet, const Node& var);

    bool isIndependentMultinomialLogisticRegression(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn = NULL);

    // This takes an inordinate amount of time. -jdramsey 20150929
    arma::uvec getNonMissingRows(const Node& x, const Node& y, std::vector<Node>& z);

    bool isMissing(const Node& x, int i);

    double multiLL(arma::mat& coeffs, const Node& dep, std::vector<Node>& indep);

    bool isIndependentRegression(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn = NULL);

    /**
     * @return true if the given independence question is judged false, true if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are searchVariables in the list returned by
     * getVariableNames().
     * Optionally return the p-value into pReturn
     */
    bool isDependent(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn) { return !isIndependent(x, y, z, pReturn); }

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
    std::vector<Node> getVariables() {return this->searchVariables; }

    /**
     * @return the list of variable varNames.
     */
    std::vector<std::string> getVariableNames();

    Node getVariable(std::string name);

    /**
     * @return true if y is determined the variable in z.
     */
    bool determines(std::vector<Node>& z, const Node& y) { return false; /*stub ?????*/ }

    double getAlpha() { return this->alpha; /*STUB ?????*/ }

    void setAlpha(double alpha) { this->alpha = alpha; }

    DataSet getData() { return this->originalData; }

    double getScore() { return this->getPValue(); }

    bool isVerbose() { return this->verbose; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    arma::mat getSubsetData(DataSet& origData, std::vector<Node>& varSubset);

    int getSampleSize() { return originalData.getNumRows(); }

    std::vector<DataSet*> getDataSets() { return { &originalData }; }

    friend void indTestMultiTest(const Rcpp::DataFrame& df);

};

#endif /* INDTESTMULTI_HPP_ */
