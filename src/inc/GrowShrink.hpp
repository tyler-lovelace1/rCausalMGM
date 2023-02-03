#ifndef GROWSHRINK_HPP_
#define GROWSHRINK_HPP_

// [[Rcpp::depends(BH)]]

#include <boost/functional/hash.hpp>
// #include "armaLapack.hpp"
// #include "IndependenceTest.hpp"
// #include "CoxIRLSRegression.hpp"
#include "LogisticRegression.hpp"
#include "LinearRegression.hpp"
#include "IndTestMulti.hpp"
// #include "Node.hpp"
// #include "DataSet.hpp"
#include "BlockingQueue.hpp"
#include "RcppThread.h"
#include <mutex>
#include <thread>
#include <atomic>
#include <float.h>
#include <algorithm>
#include <iterator>
#include <functional>

class GrowShrink {
private:
    DataSet originalData;
    std::vector<Node> searchVariables;
    DataSet internalData;
    std::map<Node, std::vector<Node>> variablesPerNode;
    std::map<Node, double> nullBICmap;
    IndependenceTest *test;
    // CoxIRLSRegression coxRegression;
    LogisticRegression logisticRegression;
    LinearRegression regression;
    bool verbose = false;
    double penalty = 1;
    std::map<std::pair<Node, std::vector<Node>>, double> scoreHistory;

    // Concurrency
    struct RegressionTask {
        Node x;
        Node y;
        std::list<Node> z;
	RegressionTask() : x(),
			   y(),
			   z() {}
        RegressionTask(const Node& _x, const Node& _y, std::list<Node>& _z) : x(_x),
									      y(_y),
									      z(_z) {}
        // IndependenceTask(const IndependenceTask& it) { x = it.x; y = it.y; z = it.z; }
    };

    static const int MAX_QUEUE_SIZE = 100;
    BlockingQueue<RegressionTask> taskQueue;
    
    int parallelism;
    
    std::mutex scoreMutex;
    std::mutex historyMutex;

    void producerGrow(const Node& target, std::list<Node>& active,
		      std::list<Node>& inactive);

    void producerShrink(const Node& target, std::list<Node>& active);

    void consumerGrow(std::unordered_map<Node, double>& scoreMap);

    void consumerShrink(std::unordered_map<Node, double>& scoreMap);

    double regressBIC(const Node& target, std::vector<Node>& regressors, bool history);

public:
    GrowShrink() : taskQueue(MAX_QUEUE_SIZE) {}

    GrowShrink(DataSet& data, int threads = -1);

    // GrowShrink(IndependenceTest* test, int threads = -1);

    GrowShrink(GrowShrink& other) : taskQueue(MAX_QUEUE_SIZE) {
    	originalData = DataSet(other.originalData);
    	internalData = DataSet(other.internalData);
    	searchVariables = other.searchVariables;
    	variablesPerNode = other.variablesPerNode;
	// test = other.test;
    	// coxRegression = CoxIRLSRegression(other.coxRegression);
    	logisticRegression = LogisticRegression(other.logisticRegression);
    	regression = LinearRegression(other.regression);
    	verbose = other.verbose;
    }

    GrowShrink& operator=(GrowShrink& other) {
    	taskQueue = BlockingQueue<RegressionTask>(MAX_QUEUE_SIZE);
    	originalData = DataSet(other.originalData);
    	internalData = DataSet(other.internalData);
    	searchVariables = other.searchVariables;
    	variablesPerNode = other.variablesPerNode;
	// test = other.test;
    	// coxRegression = CoxIRLSRegression(other.coxRegression);
    	logisticRegression = LogisticRegression(other.logisticRegression);
    	regression = LinearRegression(other.regression);
    	verbose = other.verbose;
    	return *this;
    }

    GrowShrink(GrowShrink&& other) : taskQueue(MAX_QUEUE_SIZE) {
    	std::swap(originalData, other.originalData);
        std::swap(internalData, other.internalData);
        std::swap(searchVariables, other.searchVariables);
        std::swap(variablesPerNode, other.variablesPerNode);
	// std::swap(test, other.test);
    	// coxRegression = std::move(other.coxRegression);
    	std::swap(logisticRegression, other.logisticRegression);
        std::swap(regression, other.regression);
        std::swap(verbose, other.verbose);
    }

    GrowShrink& operator=(GrowShrink&& other) {
    	taskQueue = BlockingQueue<RegressionTask>(MAX_QUEUE_SIZE);
        std::swap(originalData, other.originalData);
        std::swap(internalData, other.internalData);
        std::swap(searchVariables, other.searchVariables);
        std::swap(variablesPerNode, other.variablesPerNode);
	// std::swap(test, other.test);
    	// coxRegression = std::move(other.coxRegression);
    	std::swap(logisticRegression, other.logisticRegression);
        std::swap(regression, other.regression);
        std::swap(verbose, other.verbose);
    	return *this;
    }

    

    std::list<Node> search(const Node& target, std::vector<Node>& regressors, double* bicReturn = NULL);

    std::list<Node> grow(const Node& target, std::vector<Node>& regressors, double* bicReturn = NULL);

    // std::list<Node> search(const Node& target, std::list<Node>& regressors,
    // 			   std::list<Node>::iterator start, std::list<Node>::iterator stop,
    // 			   double* bicReturn = NULL);

    // std::list<Node> grow(const Node& target, std::list<Node>& regressors,
    // 			 std::list<Node>::iterator start, std::list<Node>::iterator stop,
    // 			 double* bicReturn = NULL);

    std::list<Node> shrink(const Node& target, std::list<Node>& active, double score, double* bicReturn = NULL);


    std::list<Node> searchSingle(const Node& target, std::vector<Node>& regressors, double* bicReturn = NULL);

    std::list<Node> growSingle(const Node& target, std::vector<Node>& regressors, double* bicReturn = NULL);

    std::list<Node> shrinkSingle(const Node& target, std::list<Node>& active, double score, double* bicReturn = NULL);

    // private synchronized List<Node> expandVariable(DataSet dataSet, const Node& node)
    std::vector<Node> expandVariable(DataSet& dataSet, const Node& var);

    // bool isIndependentMultinomialLogisticRegression(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn = NULL);

    // // This takes an inordinate amount of time. -jdramsey 20150929
    // arma::uvec getNonMissingRows(const Node& x, const Node& y, std::vector<Node>& z);

    // bool isMissing(const Node& x, int i);

    double multiLL(arma::mat& coeffs, const Node& dep, std::vector<Node>& indep);

    // bool isIndependentRegression(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn = NULL);

    // bool isIndependentCoxRegression(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn = NULL);

    /**
     * @return true if the given independence question is judged false, true if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are searchVariables in the list returned by
     * getVariableNames().
     * Optionally return the p-value into pReturn
     */
    // bool isDependent(const Node& x, const Node& y, std::vector<Node>& z, double* pReturn) { return !isIndependent(x, y, z, pReturn); }

    /**
     * @return the probability associated with the most recently executed independence test, of Double.NaN if p value is
     * not meaningful for tis test.
     */        //STUB ????????????
    // double getPValue() { return this->lastP; }

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
    // bool determines(std::vector<Node>& z, const Node& y) { return false; /*stub ?????*/ }

    // double getAlpha() { return this->alpha; /*STUB ?????*/ }

    // void setAlpha(double alpha) { this->alpha = alpha; }

    DataSet getData() { return this->originalData; }

    // double getScore() { return this->getPValue(); }

    bool isVerbose() { return this->verbose; }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setPenalty(double penalty) { this->penalty = penalty; }

    arma::mat getSubsetData(DataSet& origData, std::vector<Node>& varSubset);

    int getSampleSize() { return originalData.getNumRows(); }

    std::vector<DataSet*> getDataSets() { return { &originalData }; }

    void clearHistory() { scoreHistory.clear(); }

    // void scaledZ(const Node& target, std::vector<Node>& regressors);

    // friend void indTestMultiTest(const Rcpp::DataFrame& df);

    // friend void ScaledZTest(const Rcpp::DataFrame& df);
};

#endif /* GROWSHRINK_HPP_ */
