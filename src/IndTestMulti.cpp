// [[Rcpp::depends(BH, RcppThread)]]

#include "IndTestMulti.hpp"
#include "RcppThread.h"
#include <boost/math/distributions/chi_squared.hpp>

#include <fstream>

// const Node& IndTestMulti::nullconst Node& = Node();

IndTestMulti::IndTestMulti(DataSet& data, double alpha)
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
    
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;
    this->preferLinear = false;

    // Rcpp::Rcout << "regressions created\n";
}

IndTestMulti::IndTestMulti(DataSet& data, double alpha, bool preferLinear)
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

    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;
}

// IndTestMulti::~IndTestMulti() {
//     // Delete expanded variables
//     // for (Nodev : internalData.getVariables()) {
//     //     if (v->getName().find("MULTINOM.") != std::string::npos) {
//     //         delete v;
//     //     }
//     // }
// }

int IndTestMulti::reset()
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
bool IndTestMulti::isIndependent(const Node& x, const Node& y, std::vector<Node> &z, double *pReturn)
{
    this->timesCalled++;

    // Rcpp::Rcout << "X: " << x->getName() << std::endl;
    // Rcpp::Rcout << "Y: " << y->getName() << std::endl;

    // Rcpp::Rcout << "Z: ";
    // for (const Node& zVar : z)
    //     Rcpp::Rcout << zVar->getName() << " ";
    // Rcpp::Rcout << std::endl;
    // bool debug = (x->getName() == "X5" && y->getName() == "X6");

    if (x.isDiscrete())
    {
        // if (debug) Rcpp::Rcout << "Path 1" << std::endl;
	try {
	    return isIndependentMultinomialLogisticRegression(x, y, z, pReturn);
	}  catch (...) {
	    if (y.isDiscrete()) {
		try {
		    // Rcpp::Rcout << "Trying reverse MultinomialLogisticRegression test" << std::endl;
		    return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
		} catch (...) {
		    return false;
		}
	    }
	    return false;
	}
    }
    else if (y.isDiscrete())
    {
	try {
	    if (preferLinear) {
		// if (debug) Rcpp::Rcout << "Path 2" << std::endl;
		return isIndependentRegression(x, y, z, pReturn);
	    }
	    else {
		// if (debug) Rcpp::Rcout << "Path 3" << std::endl;
		return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
	    }
	} catch (...) {
	    return false;
	}
    }
    else
    {
        // if (debug) Rcpp::Rcout << "Path 4" << std::endl;
	try {
	    return isIndependentRegression(x, y, z, pReturn);
	}  catch (...) {
	    try {
		// Rcpp::Rcout << "Trying reverse LinearRegression test" << std::endl;
		return isIndependentRegression(y, x, z, pReturn);
	    } catch (...) {
		return false;
	    }
	    return false;
	}
    }
}

std::vector<Node> IndTestMulti::expandVariable(DataSet &dataSet, const Node& var)
{
    if (var.isContinuous())
    {
        std::vector<Node> contList;
        contList.push_back(var);
        return contList;
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

bool IndTestMulti::isIndependentMultinomialLogisticRegression(const Node& x, const Node& y, std::vector<Node>& z, double *pReturn)
{
    // std::ofstream logfile;
    // logfile.open("../test_results/debug.log", std::ios_base::app);

    if (variablesPerNode.count(x) < 1)
    {
        throw std::invalid_argument("Unrecogized variable: " + x.getName());
    }

    if (variablesPerNode.count(y) < 1)
    {
        throw std::invalid_argument("Unrecogized variable: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(varZ) < 1)
        {
            throw std::invalid_argument("Unrecogized variable: " + varZ.getName());
        }
    }

    // bool debug = (x->getName() == "X5" && y->getName() == "X6");

    arma::vec pValues;

    arma::uvec rows_ = getNonMissingRows(x, y, z);
    // logisticRegression.setRows(rows_);

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

    // //double[][] coeffsDep = new double[variablesPerNode.get(x).size()][];
    arma::mat coeffsNull = arma::mat(); //zList.size()+1, variablesPerNode.at(x).size());
    arma::mat coeffsDep = arma::mat();  //yzList.size()+1, variablesPerNode.at(x).size());

    /*********************************************************************/
    for (int i = 0; i < variablesPerNode.at(x).size(); i++)
    {
        const Node& varX = variablesPerNode.at(x).at(i);
	// LogisticRegressionResult result0;
	// LogisticRegressionResult result1;
	// try {
	LogisticRegressionResult result0 = logisticRegression.regress(varX, zList, rows_);
	LogisticRegressionResult result1 = logisticRegression.regress(varX, yzList, rows_);
	// } catch (...) {
	//     return false;
	// }

        coeffsNull.insert_cols(i, result0.getCoefs());
        coeffsDep.insert_cols(i, result1.getCoefs());

        // delete result0;
        // delete result1;
    }
    /*********************************************************************/

    // Rcpp::Rcout << "coeffsDep:" << std::endl << "Intercept\t";
    // for (const Node& v1 : yzList)
    //     for (const Node& v2 : variablesPerNode.at(v1))
    //         Rcpp::Rcout << v2->getName() << "\t";
    // Rcpp::Rcout << std::endl;
    // for (int i = 0; i < coeffsDep.n_cols; i++) {
    //     for (int j = 0; j < coeffsDep.n_rows; j++) {
    // 	    Rcpp::Rcout << coeffsDep(j,i) << "\t";
    // 	}
    // 	Rcpp::Rcout << std::endl;
    // }
    // Rcpp::Rcout << std::endl;

    // Rcpp::Rcout << "coeffsNull:" << std::endl << "Intercept\t";
    // for (const Node& v1 : zList)
    //     for (const Node& v2 : variablesPerNode.at(v1))
    //         Rcpp::Rcout << v2->getName() << "\t";
    // Rcpp::Rcout << std::endl;
    // for (int i = 0; i < coeffsNull.n_cols; i++) {
    //     for (int j = 0; j < coeffsNull.n_rows; j++) {
    // 	    Rcpp::Rcout << coeffsNull(j,i) << "\t";
    // 	}
    // 	Rcpp::Rcout << std::endl;
    // }
    // Rcpp::Rcout << std::endl;

    double ll = multiLL(coeffsDep, x, yzList);
    double ll0 = multiLL(coeffsNull, x, zList);
    double chisq = std::max(2*(ll - ll0), 1e-15);

    // if ((std::isinf(ll) && std::isinf(ll0)) || (ll0 > ll))
    // {
    //     chisq = 1e-10;
    // }
    // else if (std::isinf(ll0))
    // {
    //     chisq = 1e20;
    // }
    // else
    // {
    //     chisq = 2 * (ll - ll0); // Need to make multiLL
    // }

    int df = variablesPerNode.at(y).size() * variablesPerNode.at(x).size();
    boost::math::chi_squared dist(df);

    // if (std::isnan(chisq)) {
    // 	logfile << "IND TEST" << std::endl;
    // 	// // logfile << "dist.df = " << dist.degrees_of_freedom() << std::endl;
    // 	logfile << "chisq = " << chisq << std::endl << std::endl;
    // }

    double p = 1.0 - cdf(dist, chisq);

    if (pReturn != NULL)
        *pReturn = p;

    // if(debug) {
    //     Rcpp::Rcout << "DUDEK_DEBUG_MultiLogistic" << std::endl;
    //     Rcpp::Rcout << "chisq value = " << chisq << std::endl;
    //     Rcpp::Rcout << "dist.df = " << dist.degrees_of_freedom() << std::endl;
    // }

    // // double p = 1.0;
    //
    // // Choose the minimum of the p-values
    // // This is only one method that can be used, this requires every coefficient to be significant
    // //for (double val : pValues) {
    // //    if (val < p) p = val;
    // //}

    bool indep = p > alpha;

    // TODO - write this to a log file
    // if (indep)
    //     logfile << x->getName() << " _||_ " << y->getName() << " | { ";
    // else
    //     logfile << x->getName() << " _|/_ " << y->getName() << " | { ";
    // for (const Node& zVar : z)
    //     logfile << zVar->getName() << " ";
    // logfile << "}";

    // logfile << " with p = " << p << std::endl;

    // logfile.close();

    this->lastP = p;

    // //t.println(x + " is independent of " + y + " given " + z + ": " + indep);
    return indep;
}

// This takes an inordinate amount of time. -jdramsey 20150929
arma::uvec IndTestMulti::getNonMissingRows(const Node& x, const Node& y, std::vector<Node>& z)
{

    // std::vector<int> rows = std::vec(internalData.getNumRows());
    arma::uvec rows = arma::uvec(internalData.getNumRows());
    for (arma::uword k = 0; k < rows.size(); k++)
        rows[k] = k; // use uword or use double?

    return rows;
}

bool IndTestMulti::isMissing(const Node& x, int i)
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

double IndTestMulti::multiLL(arma::mat &coeffs, const Node& dep, std::vector<Node> &indep)
{

    // if (dep == nullNode)
    //     throw std::invalid_argument("must have a dependent node to regress on!");
    // std::ofstream logfile;
    // logfile.open("../test_results/debug.log", std::ios_base::app);

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

    // double p;
    // for (arma::uword i = 0; i < probs.n_cols; i++) {
    // 	for (arma::uword j = 0; j < probs.n_rows; j++) {
    // 	    p = probs(i,j);
    // 	    probs(i,j) = std::abs(p) < 200 ? p : ((p > 0) - (p < 0)) * 200;
    // 	}
    // }

    probs.insert_cols(0, arma::mat(indepData.n_rows, 1, arma::fill::ones));
    // probs = arma::exp(probs);

    // Rcpp::Rcout << "probs = " << probs << std::endl;
    double ll = 0;
    for (int i = 0; i < N; i++)
    {
        double b = probs.row(i).max();
        arma::rowvec curRow = probs.row(i) - b;
        //b = curRow.max();
        // double p;
        // for (arma::uword j = 0; j < curRow.n_elem; j++) {
        //     // p = curRow(j);
        //     curRow(j) = curRow(j) > -200 ? curRow(j) : -200;
        // }
        curRow = arma::exp(curRow);
        double sum = arma::sum(curRow);
        // Rcpp::Rcout << "curRow " << i << " = " << curRow << std::endl;
        // Rcpp::Rcout << "sum(curRow) " << i << " = " << sum << std::endl;
        curRow = curRow / sum;
        // Rcpp::Rcout << "curRow " << i << " = " << curRow << std::endl;
        // Rcpp::Rcout << "loglikelihood " << i << " = " << std::log(curRow.at((int)depData.at(i))) << std::endl;
        ll += std::log(curRow.at((int)depData.at(i)));
    }
    // Rcpp::Rcout << "multiLL loglikelihood = " << ll << std::endl;
    // logfile << "ll = " << ll << std::endl;
    // if (std::isinf(ll)) {
    // 	ll = std::numeric_limits<double>::lowest();
    // }
    if (std::isnan(ll))
    {
        ll = -std::numeric_limits<double>::infinity();
        // logfile << "Full separation found: ll = " << ll << std::endl;
        // for (auto it = indep.begin(); it != indep.end(); it++) {
        //     logfile << (*it)->getName() << "\t";
        // }
        // logfile << std::endl;
    }
    //   logfile << "Error in multiLL: ll = " << ll << std::endl;
    //   for (auto it = indep.begin(); it != indep.end(); it++) {
    // 	logfile << (*it)->getName() << "\t";
    //   }
    //   logfile << std::endl;
    //   logfile << "probs\n";
    //   logfile << probs << std::endl;
    //   // ll = std::numeric_limits<double>::min();
    // }
    // logfile.close();
    return ll;
}

bool IndTestMulti::isIndependentRegression(const Node& x, const Node& y, std::vector<Node>& z, double *pReturn)
{
    if (variablesPerNode.count(x) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + x.getName());
    }

    if (variablesPerNode.count(y) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(varZ) < 1)
        {
            throw std::invalid_argument("Unrecogized node: " + varZ.getName());
        }
    }

    // std::ofstream logfile;
    // logfile.open("itm_debug.log", std::ios_base::app);

    std::vector<Node> regressors = variablesPerNode.at(y);

    // logfile << regressors[0]->getName() << std::endl;

    for (const Node& varZ : z)
    {
	// logfile << varZ->getName() << std::endl;
        std::vector<Node> temp = variablesPerNode.at(varZ);
        regressors.insert(regressors.end(), temp.begin(), temp.end());
    }

    // logfile.close();

    arma::uvec rows = getNonMissingRows(x, y, z);
    // LinearRegression regression = this->regression;
    // regression.setRows(rows); // regression was never declared

    double p;

    // try
    // {
    RegressionResult result;

    try {
	result = regression.regress(x, regressors, rows);
    } catch (...) {
	return false;
    }

    // logfile.open("debug.log", std::ios_base::app);
    // logfile << result << std::endl;
    // logfile.close();
    
    p = result.getP().at(1); // double check on .at(1)
        // delete result;
    // }
    // catch (std::exception e)
    // {
    // 	// Rcpp::Rcout << e.what() << std::endl;
    //     return false;
    // }

    if (pReturn != NULL)
	*pReturn = p;
    
    bool indep = p > alpha;

    this->lastP = p;
    // TODO - write this to a log file
    // if (indep)
    //     Rcpp::Rcout << x->getName() << " _||_ " << y->getName() << " | { ";
    // else
    //     Rcpp::Rcout << x->getName() << " _|/_ " << y->getName() << " | { ";
    // for (const Node& zVar : z)
    //     Rcpp::Rcout << zVar->getName() << " ";
    // Rcpp::Rcout << "}";

    // Rcpp::Rcout << " with p = " << p << std::endl;

    // if (indep)
    //     logfile << x->getName() << " _||_ " << y->getName() << " | { ";
    // else
    //     logfile << x->getName() << " _|/_ " << y->getName() << " | { ";
    // for (const Node& zVar : z)
    //     logfile << zVar->getName() << " ";
    // logfile << "}";

    // logfile << " with p = " << p << std::endl;

    // logfile.close();

    return indep;
}

/**
 * @return the list of variable varNames.
 */
std::vector<std::string> IndTestMulti::getVariableNames()
{
    std::vector<Node> variables = getVariables();
    std::vector<std::string> variableNames;

    for (const Node& var1 : variables)
    {
        variableNames.push_back(var1.getName());
    }
    return variableNames;
}

Node IndTestMulti::getVariable(std::string name)
{
    for (int i = 0; i < getVariables().size(); i++)
    {
        const Node& var = getVariables().at(i);
        if (var.getName() == name)
        {
            return var;
        }
    }
    Node emptyVar;
    return emptyVar;
}

arma::mat IndTestMulti::getSubsetData(DataSet &origData, std::vector<Node> &varSubset)
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

// void indTestMultiTest(const Rcpp::DataFrame& df) {
//   Rcpp::Rcout << "*******Start******* \n";
//   DataSet data(df, 5);
//   Rcpp::Rcout << "Dataset Constructed \n";
//   IndTestMulti itm(data, 0.05);
//   Rcpp::Rcout << "Independence Test Constructed \n";
//   // std::vector<std::string> varNames = itm.getVariableNames(); // unnecessary but this would be how you would normally retrieve the data
//   const Node& x = itm.getVariable("X1");
//   const Node& y = itm.getVariable("X6");
//   Rcpp::Rcout << "FIRST X: " << x->getName() << " Continuous: " << x->isContinuous() << std::endl;
//   Rcpp::Rcout << "FIRST Y: " << y->getName() << " Continuous: " << y->isContinuous() << std::endl;
//   std::vector<Node> z;
//   Rcpp::Rcout << "FIRST CASE: " << itm.isIndependent(x,y,z) << std::endl;
//   Rcpp::Rcout << "First P: " << itm.getPValue() << std::endl;
//   z.push_back(itm.getVariable("X9"));
//   Rcpp::Rcout << "SECOND CASE: " << itm.isIndependent(x,y,z) << std::endl;
//   Rcpp::Rcout << "Second P: " << itm.getPValue() << std::endl;
//   z.push_back(itm.getVariable("X10"));
//   z.push_back(itm.getVariable("X8"));
//   Rcpp::Rcout << "THIRD CASE: " << itm.isIndependent(x,y,z) << std::endl;
//   Rcpp::Rcout << "Third P: " << itm.getPValue() << std::endl;
//   x = itm.getVariable("X2");
//   Rcpp::Rcout << "X IS DISCRETE CASE: " << itm.isIndependent(x,y,z) << std::endl;
//   Rcpp::Rcout << "DISCRETE CASE P: " << itm.getPValue() << std::endl;

//   // try the case where there are values for x and y but z is empty, or when z has 1 or more variables
//   //when x and y are mixed between    at least one where x is discrete
// }
