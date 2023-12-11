// [[Rcpp::depends(BH, RcppThread)]]

#include "BayesIndTestMultiCox.hpp"
#include "RcppThread.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>

#include <fstream>

// const Node& BayesIndTestMultiCox::nullconst Node& = Node();

BayesIndTestMultiCox::BayesIndTestMultiCox(DataSet& data, double alpha)
{
    this->timesCalled = 0;
    this->searchVariables = data.getVariables();
    this->originalData = DataSet(data);
    this->internalData = DataSet(data);
    this->alpha = alpha;
    this->lastP = 0;

    // Rcpp::Rcout << "dataset copied\n";

    std::vector<Node> variables = internalData.getVariables();

    for (const Node& var : variables)
    {
	// Rcpp::Rcout << var->getName() << "\n";
        std::vector<Node> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
    }

    // Rcpp::Rcout << "dataset expanded\n";

    this->coxRegression = CoxIRLSRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;
    this->preferLinear = true;

    // Rcpp::Rcout << "regressions created\n";
}

BayesIndTestMultiCox::BayesIndTestMultiCox(DataSet& data, double alpha, bool preferLinear)
{
    this->timesCalled = 0;
    this->preferLinear = preferLinear;
    this->searchVariables = data.getVariables();
    this->originalData = DataSet(data);
    this->internalData = DataSet(data);
    this->alpha = alpha;
    this->lastP = 0;

    std::vector<Node> variables = internalData.getVariables();

    for (const Node& var : variables)
    {
        //   Rcpp::Rcout << var->getName() << std::endl;
        std::vector<Node> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
    }

    this->coxRegression = CoxIRLSRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;
}

// BayesIndTestMultiCox::~BayesIndTestMultiCox() {
//     // Delete expanded variables
//     // for (Nodev : internalData.getVariables()) {
//     //     if (v->getName().find("MULTINOM.") != std::string::npos) {
//     //         delete v;
//     //     }
//     // }
// }

int BayesIndTestMultiCox::reset()
{
    this->timesCalled = 0;
    return timesCalled;
}

/**
 * @return true if the given independence question is judged true, false if not. The independence question is of the
 * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are searchVariables in the list returned by
 * getVariableNames().
 * Optionally return the p-value into pReturn
 */
bool BayesIndTestMultiCox::isIndependent(const Node& x, const Node& y, std::vector<Node> &z, double *pReturn)
{
    this->timesCalled++;

    if (x.isCensored()) {
	if (y.isCensored()) {
	    if (pReturn != NULL)
		*pReturn = alpha;
	    return true; // Assumed non-interacting
	}
	try {
	    return isIndependentCoxRegression(x, y, z);
	} catch (const std::exception &exc) {
	    if (pReturn != NULL)
		*pReturn = 0.0;
	    return false;
	}
    } else if (y.isCensored()) {
	try {
	    return isIndependentCoxRegression(y, x, z);
	} catch (const std::exception &exc) {
	    if (pReturn != NULL)
		*pReturn = 0.0;
	    return false;
	}
    } else if (x.isDiscrete()) {
	// X is discrete
	if (y.isContinuous()) {
	    // X is discrete, Y is continuous
	    try {
		if (preferLinear) {
		    return isIndependentRegression(y, x, z, pReturn);
		} else {
		    return isIndependentMultinomialLogisticRegression(x, y, z, pReturn);
		}
	    } catch (const std::exception &exc) {
		// if regression fails, try the opposite regression by reversing condition
		try {
		    if (!preferLinear) {
			return isIndependentRegression(y, x, z, pReturn);
		    } else {
			return isIndependentMultinomialLogisticRegression(x, y, z, pReturn);
		    }
		} catch (const std::exception &exc) {
		    if (pReturn != NULL)
			*pReturn = 0.0;
		    return false;
		}
		if (pReturn != NULL)
		    *pReturn = 0.0;
		return false;
	    }
	} else if (y.isDiscrete()) {
	    // X and Y are discrete
	    try {
		if (variablesPerNode.at(x).size() <= variablesPerNode.at(y).size())
		    return isIndependentMultinomialLogisticRegression(x, y, z, pReturn);
		else
		    return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
	    } catch (const std::exception &exc) {
		// if regression fails, try the opposite regression by reversing condition
		try {
		    if (variablesPerNode.at(x).size() > variablesPerNode.at(y).size())
			return isIndependentMultinomialLogisticRegression(x, y, z, pReturn);
		    else
			return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
		} catch (const std::exception &exc) {
		    if (pReturn != NULL)
			*pReturn = 0.0;
		    return false;
		}
		if (pReturn != NULL)
		    *pReturn = 0.0;
		return false;
	    }
	} else {
	    if (pReturn != NULL)
		*pReturn = 0.0;
	    return false;
	}
    } else if (y.isDiscrete()) {
	// X is continuous, Y is discrete
	try {
	    if (preferLinear) {
		return isIndependentRegression(x, y, z, pReturn);
	    }
	    else {
		return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
	    }
	} catch (const std::exception &exc) {
	    // if regression fails, try the opposite regression by reversing condition
	    try {
		if (!preferLinear) {
		    return isIndependentRegression(x, y, z, pReturn);
		}
		else {
		    return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
		}
	    } catch (const std::exception &exc) {
		if (pReturn != NULL)
		    *pReturn = 0.0;
		return false;
	    }
	    if (pReturn != NULL)
		*pReturn = 0.0;
	    return false;
	}
    } else {
	// X and Y are continuous
	try {
	    return isIndependentRegression(x, y, z, pReturn);
	}  catch (const std::exception &exc) {
	    // if regression fails, try the opposite regression
	    try {
		return isIndependentRegression(y, x, z, pReturn);
	    } catch (const std::exception &exc) {
		if (pReturn != NULL)
		    *pReturn = 0.0;
		return false;
	    }
	    if (pReturn != NULL)
		*pReturn = 0.0;
	    return false;
	}
    }
}

std::vector<Node> BayesIndTestMultiCox::expandVariable(DataSet &dataSet, const Node& var)
{
    if (var.isContinuous())
    {
        std::vector<Node> contList;
        contList.push_back(var);
        return contList;
    }

    if (var.isCensored())
    {
        std::vector<Node> censList;
        censList.push_back(var);
        return censList;
    }

    if (var.isDiscrete() && var.getNumCategories() < 3)
    {
        std::vector<Node> discList;
        discList.push_back(var);
        return discList;
    }

    if (!var.isDiscrete())
    {
        throw std::invalid_argument("*Invalid variable type*");
    }

    std::vector<std::string> varCats = var.getCategories();

    std::vector<Node> variables;
    /*********************************************************************/
    std::string temp = var.getName();
    for (auto it = varCats.begin() + 1; it != varCats.end(); it++)
    {
	const Node& newVar = Node(new DiscreteVariable(temp + "MULTINOM." + *it, 2));

        /*********************************************************************/

        variables.push_back(newVar);

        dataSet.addVariable(newVar);

        int newVarIndex = dataSet.getColumn(newVar);
        int numCases = dataSet.getNumRows();

        for (int l = 0; l < numCases; l++)
        {
            int dataCellIndex = dataSet.getInt(l, dataSet.getColumn(var));
            if (dataCellIndex == var.getIndex(*it))
            {
                dataSet.set(l, newVarIndex, 1);
            }
            else
            {
                dataSet.set(l, newVarIndex, 0);
            }
        }
    }

    return variables;
}

bool BayesIndTestMultiCox::isIndependentMultinomialLogisticRegression(const Node& x, const Node& y, std::vector<Node>& z, double *pReturn)
{
    // std::ofstream logfile;
    // logfile.open("../test_results/debug.log", std::ios_base::app);

    if (variablesPerNode.count(x) < 1)
    {
        throw std::invalid_argument("Unrecognized variable: " + x.getName());
    }

    if (variablesPerNode.count(y) < 1)
    {
        throw std::invalid_argument("Unrecognized variable: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(varZ) < 1)
        {
            throw std::invalid_argument("Unrecognized variable: " + varZ.getName());
        }
    }

    // bool debug = (x->getName() == "X5" && y->getName() == "X6");

    // std::string regressionString = "LogisticRegression::Testing " + x.getName() + " _||_ " + y.getName() + " | {";

    // if (z.size() == 0) regressionString += "}";
    // else {
    // 	for (const Node& n : z) {
    // 	    regressionString += " " + n.getName();
    // 	}
    // 	regressionString += " }";
    // }

    // RcppThread::Rcout << regressionString << std::endl;

    arma::vec pValues;

    arma::uvec rows_ = getNonMissingRows(x, y, z);
    // logisticRegression.setRows(rows_);

    double N = rows_.n_elem;

    std::vector<Node> yzList;
    std::vector<Node> zList;

    std::vector<Node> temp = variablesPerNode.at(y);
    yzList.insert(yzList.begin(), temp.begin(), temp.end());

    for (const Node& varZ : z)
    {
        temp = variablesPerNode.at(varZ);
        yzList.insert(yzList.end(), temp.begin(), temp.end());

        // temp = variablesPerNode.at(varZ);
        zList.insert(zList.end(), temp.begin(), temp.end());
    }

    LogisticRegressionResult result0, result1;
    // //double[][] coeffsDep = new double[variablesPerNode.get(x).size()][];
    arma::mat coeffsNull = arma::mat(); //zList.size()+1, variablesPerNode.at(x).size());
    arma::mat coeffsDep = arma::mat();  //yzList.size()+1, variablesPerNode.at(x).size());

    /*********************************************************************/
    for (int i = 0; i < variablesPerNode.at(x).size(); i++) {
        const Node& varX = variablesPerNode.at(x).at(i);

	try {
	    result0 = logisticRegression.regress(varX, zList, rows_);
	    result1 = logisticRegression.regress(varX, yzList, rows_);
	} catch (...) {
	    return false;
	}

        coeffsNull.insert_cols(i, result0.getCoefs());
        coeffsDep.insert_cols(i, result1.getCoefs());

    }

    double ll, ll0, bic, bic0, deltaBIC;

    ll = multiLL(coeffsDep, x, yzList);
    ll0 = multiLL(coeffsNull, x, zList);

    bic = -2 * ll + std::log(N) * coeffsDep.n_elem;
    bic0 = -2 * ll0 + std::log(N) * coeffsNull.n_elem;

    deltaBIC = bic - bic0;
    
    // RcppThread::Rcout << "ll = " << ll << std::endl;
    // RcppThread::Rcout << "bic = " << bic << std::endl;
    
    // RcppThread::Rcout << "ll0 = " << ll0 << std::endl;
    // RcppThread::Rcout << "bic0 = " << bic0 << std::endl;

    // RcppThread::Rcout << "delta bic = " << deltaBIC << std::endl;

    double p = 1.0 / (1.0 + std::exp(-deltaBIC/2));

    // RcppThread::Rcout << "pval = " << p << std::endl;
    

    if (pReturn != NULL)
        *pReturn = p;


    bool indep = R::runif(0,1) < p;


    this->lastP = p;

    return indep;
}

// This takes an inordinate amount of time. -jdramsey 20150929
arma::uvec BayesIndTestMultiCox::getNonMissingRows(const Node& x, const Node& y, std::vector<Node>& z)
{

    // std::vector<int> rows = std::vec(internalData.getNumRows());
    arma::uvec rows = arma::uvec(internalData.getNumRows());
    for (arma::uword k = 0; k < rows.size(); k++)
        rows[k] = k; // use uword or use double?

    return rows;
}

bool BayesIndTestMultiCox::isMissing(const Node& x, int i)
{
    int j = internalData.getColumn(x);

    if (x.isDiscrete())
    {
        int v = internalData.getInt(i, j);

        if (v == -99)
        {
            return true;
        }
    }

    if (x.isContinuous())
    {
        double v = internalData.getInt(i, j);

        if (std::isnan(v))
        {
            return true;
        }
    }

    return false;
}

double BayesIndTestMultiCox::multiLL(arma::mat &coeffs, const Node& dep, std::vector<Node> &indep)
{

    std::vector<Node> depList;
    depList.push_back(dep);

    int i = internalData.getColumn(dep);
    arma::mat a = internalData.getData();
    arma::vec depData = a.col(i);

    int N = depData.n_elem; // returns number of rows

    arma::mat indepData;
    if (indep.size() == 0)
        indepData = arma::mat(N, 1, arma::fill::ones); // filling it with ones
    else
    {
        indepData = getSubsetData(internalData, indep);
        indepData.insert_cols(0, arma::mat(N, 1, arma::fill::ones));
    }

    arma::mat probs = indepData * coeffs;

    probs.insert_cols(0, arma::mat(indepData.n_rows, 1, arma::fill::zeros)); // reference class

    double ll = 0;
    for (int i = 0; i < N; i++)
    {
        arma::rowvec curRow = probs.row(i); // - b;
        curRow = arma::exp(curRow);
        double sum = arma::sum(curRow);
        curRow = curRow / sum;
    	ll += std::log(curRow((int) depData.at(i)));
    }
    
    if (std::isnan(ll)) {
	ll = -std::numeric_limits<double>::infinity();
    }
    
    return ll;
}

bool BayesIndTestMultiCox::isIndependentRegression(const Node& x, const Node& y, std::vector<Node>& z, double *pReturn)
{
    if (variablesPerNode.count(x) < 1)
    {
        throw std::invalid_argument("Unrecognized node: " + x.getName());
    }

    if (variablesPerNode.count(y) < 1)
    {
        throw std::invalid_argument("Unrecognized node: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(varZ) < 1)
        {
            throw std::invalid_argument("Unrecognized node: " + varZ.getName());
        }
    }

    std::vector<Node> regressors = variablesPerNode.at(y);
    std::vector<Node> zRegressors;

    for (const Node& varZ : z) {
        std::vector<Node> temp = variablesPerNode.at(varZ);
        regressors.insert(regressors.end(), temp.begin(), temp.end());
	// if (y.isDiscrete())
	zRegressors.insert(zRegressors.end(), temp.begin(), temp.end());
    }

    arma::uvec rows = getNonMissingRows(x, y, z);

    double p, ll, ll0, bic, bic0, deltaBIC;

    RegressionResult result, result0;
    
    try {
	result = regression.regress(x, regressors, rows);
	result0 = regression.regress(x, zRegressors, rows);
    } catch (...) {
	return false;
    }

    double N = rows.n_elem;

    ll = -0.5 * N * std::log(result.getRSS() / N);
    ll0 = -0.5 * N * std::log(result0.getRSS() / N);
    
    bic = N * std::log(result.getRSS() / N) + std::log(N) * result.getCoef().n_elem;
    bic0 = N * std::log(result0.getRSS() / N) + std::log(N) * result0.getCoef().n_elem;

    deltaBIC = bic - bic0;
    
    // RcppThread::Rcout << "RSS = " << result.getRSS() << std::endl;
    // RcppThread::Rcout << "RSS0 = " << result0.getRSS() << std::endl;
    
    // RcppThread::Rcout << "ll = " << ll << std::endl;
    // RcppThread::Rcout << "bic = " << bic << std::endl;
    
    // RcppThread::Rcout << "ll0 = " << ll0 << std::endl;
    // RcppThread::Rcout << "bic0 = " << bic0 << std::endl;

    // RcppThread::Rcout << "delta bic = " << deltaBIC << std::endl;

    p = 1.0 / (1.0 + std::exp(-deltaBIC/2));

    // RcppThread::Rcout << "pval = " << p << std::endl;

    if (pReturn != NULL)
	*pReturn = p;
    
    bool indep = R::runif(0,1) < p;

    this->lastP = p;

    return indep;
}

bool BayesIndTestMultiCox::isIndependentCoxRegression(const Node& x, const Node& y, std::vector<Node> &z, double *pReturn)
{
    if (variablesPerNode.count(internalData.getVariable(x.getName())) < 1)
    {
        throw std::invalid_argument("Unrecognized node: " + x.getName());
    }

    if (variablesPerNode.count(internalData.getVariable(y.getName())) < 1)
    {
        throw std::invalid_argument("Unrecognized node: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)
        {
            throw std::invalid_argument("Unrecognized node: " + varZ.getName());
        }
    }
    
    std::vector<Node> zList;
    std::vector<Node> regressors;
    std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(y.getName()));
    regressors.insert(regressors.end(), temp.begin(), temp.end());

    for (const Node& varZ : z) {
        temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));
	if (varZ.isCensored()) {
	    // regressors.insert(regressors.end(), temp.begin()+1, temp.end());
	    // if (y.isDiscrete()) {
	    // 	zList.insert(zList.end(), temp.begin()+1, temp.end());
	    // }
	    throw std::runtime_error("Related censored variable: " + varZ.getName());
	} else {
	    regressors.insert(regressors.end(), temp.begin(), temp.end());
	    if (y.isDiscrete()) {
		zList.insert(zList.end(), temp.begin(), temp.end());
	    }
	}
    }

    arma::uvec rows = getNonMissingRows(x, y, z);

    double p, ll, ll0, bic, bic0, deltaBIC;

    CoxRegressionResult result, result0;

    try {
	result = coxRegression.regress(x, regressors, rows);
	result0 = coxRegression.regress(x, zList, rows);
    } catch (...) {
	return false;
    }

    double N = x.getNEvents();

    ll = result.getLoglikelihood();
    ll0 = result0.getLoglikelihood();

    bic = -2 * ll + std::log(N) * result.getCoef().n_elem;
    bic0 = -2 * ll0 + std::log(N) * result0.getCoef().n_elem;

    deltaBIC = bic - bic0;
    
    // RcppThread::Rcout << "ll = " << ll << std::endl;
    // RcppThread::Rcout << "bic = " << bic << std::endl;
    
    // RcppThread::Rcout << "ll0 = " << ll0 << std::endl;
    // RcppThread::Rcout << "bic0 = " << bic0 << std::endl;

    // RcppThread::Rcout << "delta bic = " << deltaBIC << std::endl;

    p = 1.0 / (1.0 + std::exp(-deltaBIC/2));

    // RcppThread::Rcout << "pval = " << p << std::endl;

    if (pReturn != NULL)
	*pReturn = p;
    
    bool indep = R::runif(0,1) < p;


    this->lastP = p;

    return indep;
}

/**
 * @return the list of variable varNames.
 */
std::vector<std::string> BayesIndTestMultiCox::getVariableNames()
{
    std::vector<Node> variables = getVariables();
    std::vector<std::string> variableNames;

    for (const Node& var1 : variables)
    {
        variableNames.push_back(var1.getName());
    }
    return variableNames;
}

Node BayesIndTestMultiCox::getVariable(std::string name)
{
    std::vector<Node> variables = getVariables();
    for (int i = 0; i < variables.size(); i++)
    {
        const Node& var = variables.at(i);
        if (var.getName() == name)
        {
            return var;
        }
    }
    Node emptyVar;
    return emptyVar;
}

arma::mat BayesIndTestMultiCox::getSubsetData(DataSet &origData, std::vector<Node> &varSubset)
{
    arma::mat origMat = origData.getData();
    arma::uvec colIndices(varSubset.size());
    arma::uvec rowIndices(origMat.n_rows);
    // for (const Node& var : varSubset){
    for (int i = 0; i < varSubset.size(); i++)
    {
        const Node& var = varSubset.at(i);
        arma::uword j = origData.getColumn(var);
        colIndices[i] = j;
    }
    for (int i = 0; i < origMat.n_rows; i++)
    {
        rowIndices[i] = i;
    }
    return origMat.submat(rowIndices, colIndices);
}

void BayesIndTestMultiCox::resetWZ(Node target, std::vector<Node>& neighbors) {
    if (variablesPerNode.count(internalData.getVariable(target.getName())) < 1)
    {
        throw std::invalid_argument("Unrecognized node: " + target.getName());
    }

    for (const Node& var : neighbors)
    {
        if (variablesPerNode.count(internalData.getVariable(var.getName())) < 1)
        {
            throw std::invalid_argument("Unrecognized node: " + var.getName());
        }
    }

    std::vector<Node> temp;
    std::vector<Node> regressors;
    std::vector<Node> emptySet = {};

    for (const Node& var : neighbors) {
        temp = variablesPerNode.at(internalData.getVariable(var.getName()));
	regressors.insert(regressors.end(), temp.begin(), temp.end());
    }

    // Rcpp::Rcout << "Reset WZ for " + target.getName() << std::endl;

    CoxRegressionResult result;

    try {
	result = coxRegression.regress(target, regressors);
    } catch (...) {
	result = coxRegression.regress(target, emptySet);
    }
    // Rcpp::Rcout << result << std::endl;

    arma::vec WZ(result.getResid());

    std::vector<std::string> _neighbors;

    for (const Node& var : neighbors) {
	_neighbors.push_back(var.getName());
    }

    target.setNeighbors(_neighbors);
    target.setWZ(WZ);

    // std::map<std::pair<Node,Node>, arma::vec> WZmap;

    // for (int i = 0; i < neighbors.size(); i++) {
    // 	regressors.clear();
    // 	for (const Node& var : neighbors) {
    // 	    if (var == neighbors[i]) continue;
    // 	    temp = variablesPerNode.at(internalData.getVariable(var.getName()));
    // 	    regressors.insert(regressors.end(), temp.begin(), temp.end());
    // 	}

    // 	try {
    // 	    result = coxRegression.regress(target, regressors);
    // 	} catch (...) {
    // 	    result = coxRegression.regress(target, emptySet);
    // 	}

    // 	WZmap.insert(std::pair<std::pair<Node,Node>, arma::vec>(std::minmax(target, neighbors[i]), result.getResid()));
    // }

    // Rcpp::Rcout << "WZmap complete for Node " + target.getName() + "\n";

    if (internalData.updateNode(target)) {
	this->coxRegression = CoxIRLSRegression(internalData);
	this->logisticRegression = LogisticRegression(internalData);
	// this->logisticRegression.setWZmap(WZmap);
	this->regression = LinearRegression(internalData);
	// this->regression.setWZmap(WZmap);
    }
}

// [[Rcpp::export]]
void BayesIndTestMultiCoxTest(const std::string& x,
			      const std::string& y,
			      const Rcpp::StringVector& z,
			      const Rcpp::DataFrame& df) {
  Rcpp::Rcout << "*******Start******* \n";
  DataSet data(df, 5);
  Rcpp::Rcout << "Dataset Constructed \n";
  BayesIndTestMultiCox itm(data, 0.05);
  Rcpp::Rcout << "Independence Test Constructed \n";
  
  // std::vector<std::string> varNames = itm.getVariableNames(); // unnecessary but this would be how you would normally retrieve the data
  Node _x = itm.getVariable(x);
  Node _y = itm.getVariable(y);
  Rcpp::Rcout << "FIRST X: " << _x.getName() << " Continuous: " << _x.isContinuous() << std::endl;
  Rcpp::Rcout << "FIRST Y: " << _y.getName() << " Continuous: " << _y.isContinuous() << std::endl;
  double pval;
  std::vector<Node> _z;
  if (z.size() > 0) {
      for (int i = 0; i < z.size(); i++) {
	  std::string zVarName = Rcpp::as<std::string>(z(i));
	  _z.push_back(itm.getVariable(zVarName));
      }
  }
  Rcpp::Rcout << "Reject Null: " << itm.isIndependent(_x,_y,_z,&pval) << std::endl;
  Rcpp::Rcout << "p value: " << pval << std::endl;
  // z.push_back(itm.getVariable("Y5"));
  // Rcpp::Rcout << "SECOND CASE: " << itm.isIndependent(x,y,z,&pval) << std::endl;
  // Rcpp::Rcout << "Second P: " << pval << std::endl;
  // z.push_back(itm.getVariable("X2"));
  // z.push_back(itm.getVariable("Y3"));
  // Rcpp::Rcout << "THIRD CASE: " << itm.isIndependent(x,y,z,&pval) << std::endl;
  // Rcpp::Rcout << "Third P: " << pval << std::endl;
  // x = itm.getVariable("X2");
  // Rcpp::Rcout << "X IS DISCRETE CASE: " << itm.isIndependent(x,y,z) << std::endl;
  // Rcpp::Rcout << "DISCRETE CASE P: " << itm.getPValue() << std::endl;

  // try the case where there are values for x and y but z is empty, or when z has 1 or more variables
  //when x and y are mixed between    at least one where x is discrete
}
