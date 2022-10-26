// [[Rcpp::depends(BH,RcppThread)]]

#include "IndTestMultiScaledZ.hpp"
#include "RcppThread.h"
#include <boost/math/distributions/chi_squared.hpp>

#include <fstream>

IndTestMultiScaledZ::IndTestMultiScaledZ(DataSet& data, double alpha)
{
    this->timesCalled = 0;
    this->searchVariables = data.getVariables();
    this->originalData = DataSet(data);
    this->internalData = DataSet(data);
    // DataSet originData(data);
    // this->originalData = originData;
    // DataSet internData(data);
    // this->internalData = internData;
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
    this->preferLinear = false;

    // Rcpp::Rcout << "regressions created\n";
}

IndTestMultiScaledZ::IndTestMultiScaledZ(DataSet& data, double alpha, bool preferLinear)
{
    this->timesCalled = 0;
    this->preferLinear = preferLinear;
    this->searchVariables = data.getVariables();
    DataSet originalData(data);
    this->originalData = originalData;
    DataSet internalData(data);
    this->alpha = alpha;
    this->lastP = 0;

    std::vector<Node> variables = internalData.getVariables();

    for (const Node& var : variables)
    {
        //   Rcpp::Rcout << var->getName() << std::endl;
        std::vector<Node> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
    }

    this->internalData = internalData;
    this->coxRegression = CoxIRLSRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;
}

IndTestMultiScaledZ::~IndTestMultiScaledZ() {
    // Delete expanded variables
    // for (Nodev : internalData.getVariables()) {
    //     if (v->getName().find("MULTINOM.") != std::string::npos) {
    //         delete v;
    //     }
    // }
}

int IndTestMultiScaledZ::reset()
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
bool IndTestMultiScaledZ::isIndependent(const Node& x, const Node& y, std::vector<Node> &z, double *pReturn)
{
    this->timesCalled++;

    if (x.isCensored()) {
	if (y.isCensored()) {
	    return true; // Assumed non-interacting
	}
	try {
	    return isIndependentCoxRegression(x, y, z);
	} catch (...) {
	    return false;
	}
    } else if (y.isCensored()) {
	try {
	    return isIndependentCoxRegression(y, x, z);
	} catch (...) {
	    return false;
	}
    } else if (x.isDiscrete()) {
	try {
	    return isIndependentMultinomialLogisticRegression(x, y, z, pReturn);
	}  catch (...) {
	    if (y.isDiscrete()) {
		try {
		    return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
		} catch (...) {
		    return false;
		}
	    }
	    return false;
	}
    }
    else if (y.isDiscrete()) {
	try {
	    if (preferLinear) {
		return isIndependentRegression(x, y, z, pReturn);
	    }
	    else {
		return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
	    }
	} catch (...) {
	    return false;
	}
    }
    else {
	try {
	    return isIndependentRegression(x, y, z, pReturn);
	}  catch (...) {
	    try {
		return isIndependentRegression(y, x, z, pReturn);
	    } catch (...) {
		return false;
	    }
	    return false;
	}
    }
}

std::vector<Node> IndTestMultiScaledZ::expandVariable(DataSet &dataSet, const Node& var)
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
	
	std::string temp = var.getName();
	Node newVar = Node(new ContinuousVariable(temp + ".Z"));
	censList.push_back(newVar);
	dataSet.addVariable(newVar);

	int newVarIndex = dataSet.getColumn(newVar);
        int numCases = dataSet.getNumRows();

	CoxIRLSRegression coxIRLS(dataSet);

	std::vector<Node> regressors = {};

	CoxRegressionResult result = coxIRLS.regress(var, regressors);

	arma::vec scaledZ = result.getResid();

	// RcppThread::Rcout << "Scaled Z " << var.getName() << ": " << scaledZ.t() << std::endl;

        for (int i = 0; i < numCases; i++) {
            dataSet.set(i, newVarIndex, scaledZ[i]);
        }
	
        return censList;
    }

    if (var.isDiscrete() && var.getNumCategories() < 3)
    {
        std::vector<Node> discList;
        discList.push_back(var);
        return discList;
    }

    if (!(var.isDiscrete()))
    {
        throw std::invalid_argument("*Invalid variable type*");
    }

    std::vector<std::string> varCats = var.getCategories();

    std::vector<Node> variables;
    /*********************************************************************/
    std::string temp = var.getName();
    for (auto it = varCats.begin() + 1; it != varCats.end(); it++)
    {
	Node newVar = Node(new DiscreteVariable(temp + "MULTINOM." + *it, 2));

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

bool IndTestMultiScaledZ::isIndependentMultinomialLogisticRegression(const Node& x, const Node& y, std::vector<Node> &z, double *pReturn)
{
    // std::ofstream logfile;
    // logfile.open("../test_results/debug.log", std::ios_base::app);

    if (variablesPerNode.count(internalData.getVariable(x.getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized variable: " + x.getName());
    }

    if (variablesPerNode.count(internalData.getVariable(y.getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized variable: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)
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

    std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(y.getName()));
    yzList.insert(yzList.begin(), temp.begin(), temp.end());

    for (const Node& varZ : z) {
        temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));

	if (varZ.isCensored()) {
	    yzList.insert(yzList.end(), temp.begin()+1, temp.end());

	    // temp = variablesPerNode.at(varZ);
	    zList.insert(zList.end(), temp.begin()+1, temp.end());	    
	} else {
	    yzList.insert(yzList.end(), temp.begin(), temp.end());

	    // temp = variablesPerNode.at(varZ);
	    zList.insert(zList.end(), temp.begin(), temp.end());
	}
    }

    // //double[][] coeffsDep = new double[variablesPerNode.get(x).size()][];
    arma::mat coeffsNull = arma::mat(); //zList.size()+1, variablesPerNode.at(x).size());
    arma::mat coeffsDep = arma::mat();  //yzList.size()+1, variablesPerNode.at(x).size());

    /*********************************************************************/
    for (int i = 0; i < variablesPerNode.at(internalData.getVariable(x.getName())).size(); i++)
    {
        const Node& varX = variablesPerNode.at(internalData.getVariable(x.getName())).at(i);
	// LogisticRegressionResult result0;
	// LogisticRegressionResult result1;
	// try {
	LogisticRegressionResult result0 = logisticRegression.regress(varX,
								      zList, rows_);
	LogisticRegressionResult result1 = logisticRegression.regress(varX,
								      yzList, rows_);
	// } catch (std::runtime_error& e) {
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

    int df = variablesPerNode.at(internalData.getVariable(y.getName())).size()
	* variablesPerNode.at(internalData.getVariable(x.getName())).size();
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
arma::uvec IndTestMultiScaledZ::getNonMissingRows(const Node& x, const Node& y, std::vector<Node> &z)
{

    // std::vector<int> rows = std::vec(internalData.getNumRows());
    arma::uvec rows = arma::uvec(internalData.getNumRows());
    for (arma::uword k = 0; k < rows.size(); k++)
        rows[k] = k; // use uword or use double?

    // std::vector<Node> yzList;
    // std::vector<Node> zList;

    // std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(y->getName()));
    // yzList.insert(yzList.begin(), temp.begin(), temp.end());

    // for (NodevarZ : z) {
    //     temp = variablesPerNode.at(internalData.getVariable(varZ->getName()));
    //     yzList.insert(yzList.end(), temp.begin(), temp.end());
    // }

    // for (Nodevar : yzList)
    // 	if (var->isCensored())
    // 	    rows = arma::intersect(rows, ((CensoredNode)var)->getCC());

    return rows;
}

bool IndTestMultiScaledZ::isMissing(const Node& x, int i)
{
    int j = internalData.getColumn(x);

    if (x.isDiscrete())
    {
        int v = (int)internalData.getInt(i, j);

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

double IndTestMultiScaledZ::multiLL(arma::mat &coeffs, const Node& dep, std::vector<Node> &indep)
{

    if (dep.getName() == "??")
        throw std::invalid_argument("must have a dependent node to regress on!");
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
    else {
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

bool IndTestMultiScaledZ::isIndependentRegression(const Node& x, const Node& y, std::vector<Node> &z, double *pReturn)
{
    if (variablesPerNode.count(internalData.getVariable(x.getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + x.getName());
    }

    if (variablesPerNode.count(internalData.getVariable(y.getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)
        {
            throw std::invalid_argument("Unrecogized node: " + varZ.getName());
        }
    }

    // std::ofstream logfile;
    // logfile.open("itm_debug.log", std::ios_base::app);

    std::vector<Node> regressors;
    regressors.push_back(internalData.getVariable(y.getName()));

    // logfile << regressors[0]->getName() << std::endl;

    for (const Node& varZ : z)
    {
	// logfile << varZ->getName() << std::endl;
        std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));

	if (varZ.isCensored()) {
	    regressors.insert(regressors.end(), temp.begin()+1, temp.end());	    
	} else {
	    regressors.insert(regressors.end(), temp.begin(), temp.end());
	}
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


bool IndTestMultiScaledZ::isIndependentCoxRegression(const Node& x, const Node& y, std::vector<Node> &z, double *pReturn)
{
    if (variablesPerNode.count(internalData.getVariable(x.getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + x.getName());
    }

    if (variablesPerNode.count(internalData.getVariable(y.getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + y.getName());
    }

    for (const Node& varZ : z)
    {
        if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)
        {
            throw std::invalid_argument("Unrecogized node: " + varZ.getName());
        }
    }

    // std::ofstream logfile;
    // logfile.open("../itmsz_debug.log", std::ios_base::app);
    
    std::vector<Node> zList;
    std::vector<Node> regressors;
    std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(y.getName()));
    regressors.insert(regressors.end(), temp.begin(), temp.end());

    for (const Node& varZ : z) {
        temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));
	if (varZ.isCensored()) {
	    regressors.insert(regressors.end(), temp.begin()+1, temp.end());
	    if (y.isDiscrete()) {
		zList.insert(zList.end(), temp.begin()+1, temp.end());
	    }	    
	} else {
	    regressors.insert(regressors.end(), temp.begin(), temp.end());
	    if (y.isDiscrete()) {
		zList.insert(zList.end(), temp.begin(), temp.end());
	    }
	}
    }

    arma::uvec rows = getNonMissingRows(x, y, z);

    // logfile << "Testing:\n" << x->getName() << "  _||_  " << y->getName() << "  |  {";

    // for (NodevarZ : z) {
    // 	logfile << varZ->getName() << ", ";
    // }
    
    // logfile << "}" << std::endl;

    double p;

    CoxRegressionResult result, result0;

    try {
	result = coxRegression.regress(x, regressors, rows);
	if (y.isDiscrete()) {
	    result0 = coxRegression.regress(x, zList, rows);
	}
    } catch (...) {
	return false;
    }

    // logfile << "Cox IRLS Result:\n" << result << std::endl;

    // if (y->isDiscrete()) {
    // 	logfile << "Cox IRLS Null Result:\n" << result0 << std::endl;
    // }

    if (y.isDiscrete()) {
	double ll = result.getLoglikelihood();
	double ll0 = result0.getLoglikelihood();
	double chisq = std::max(2*(ll - ll0), 1e-15);

	// logfile << result << std::endl;
	// logfile << result0 << std::endl;

	// logfile << "2 ( " << ll << " - " << ll0 << " ) = " << chisq << std::endl;

	int df = variablesPerNode.at(internalData.getVariable(y.getName())).size();
	boost::math::chi_squared dist(df);

        p = 1.0 - cdf(dist, chisq);

	// logfile << "p = " << p << std::endl;
    } else {
      p = result.getP().at(0); 
    }
    
    if (pReturn != NULL)
	*pReturn = p;
    
    bool indep = p > alpha;

    this->lastP = p;

    // logfile.close();

    return indep;
}


/**
 * @return the list of variable varNames.
 */
std::vector<std::string> IndTestMultiScaledZ::getVariableNames()
{
    std::vector<Node> variables = getVariables();
    std::vector<std::string> variableNames;

    for (const Node& var1 : variables)
    {
        variableNames.push_back(var1.getName());
    }
    return variableNames;
}

Node IndTestMultiScaledZ::getVariable(std::string name)
{
    for (int i = 0; i < getVariables().size(); i++)
    {
        Node var = getVariables().at(i);
        if (var.getName() == name)
        {
            return var;
        }
    }
    Node emptyVar;
    return emptyVar;
}

arma::mat IndTestMultiScaledZ::getSubsetData(DataSet &origData, std::vector<Node>& varSubset) {
    arma::mat origMat = origData.getData();
    arma::uvec colIndices(varSubset.size());
    arma::uvec rowIndices(origMat.n_rows);
    
    // for (const Node& var : varSubset){
    for (int i = 0; i < varSubset.size(); i++) {
        Node var = varSubset.at(i);
        arma::uword j = origData.getColumn(var);
        colIndices[i] = j;
    }
    for (int i = 0; i < origMat.n_rows; i++) {
        rowIndices[i] = i;
    }
    
    return origMat.submat(rowIndices, colIndices);
}

void IndTestMultiScaledZ::scaledZ(const Node& target,
				  std::vector<Node>& regressors) {

    if (variablesPerNode.count(internalData.getVariable(target.getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + target.getName());
    }

    for (const Node& var : regressors)
    {
        if (variablesPerNode.count(internalData.getVariable(var.getName())) < 1)
        {
            throw std::invalid_argument("Unrecogized node: " + var.getName());
        }
    }

    // std::ofstream logfile;
    // logfile.open("../itmcc_debug.log", std::ios_base::app);

    std::vector<Node> temp;
    std::vector<Node> _regressors;
    
    for (const Node& var : regressors) {
        temp = variablesPerNode.at(internalData.getVariable(var.getName()));
	if (var.isCensored()) {
	    _regressors.insert(_regressors.end(), temp.begin()+1, temp.end());	    
	} else {
	    _regressors.insert(_regressors.end(), temp.begin(), temp.end());
	}
    }

    // RcppThread::Rcout << "Computing Efron Scaled Z for " << target->getName()
    // 		      << "..." << std::endl;
    
    std::vector<Node> tempZ(_regressors.begin()+1, _regressors.end());
    arma::uvec rows = getNonMissingRows(target, _regressors.at(0), tempZ);

    // for (Nodevar : _regressors) {
    // 	RcppThread::Rcout << var->getName() << ", ";
    // }
    // RcppThread::Rcout << std::endl;
    
    // CoxRegressionResult result = coxRegression.regress(target, _regressors, rows);

    CoxIRLSRegression coxIRLS(internalData);

    _regressors.clear();

    CoxRegressionResult result = coxIRLS.regress(target, _regressors, rows);

    // Rcpp::Rcout << result << std::endl;

    // arma::vec beta = result.getCoef();

    // int n = rows.n_elem;
    
    // arma::uvec regressors_ = arma::uvec(_regressors.size());

    // for (int i = 0; i < regressors.size(); i++) {
    //     regressors_[i] = internalData.getColumn(_regressors[i]);
    // }

    // arma::mat X = internalData.getData().submat(rows, regressors_);
    // X.insert_cols(0, arma::vec(n, arma::fill::ones));

    arma::vec scaledZ = result.getResid();

    // Rcpp::Rcout << scaledZ.t() << std::endl;

    // arma::vec beta = result.getCoef();

    int n = rows.size();
    
    // arma::uvec regressors_ = arma::uvec(_regressors.size());

    // for (int i = 0; i < regressors.size(); i++) {
    //     regressors_[i] = internalData.getColumn(_regressors[i]);
    // }

    // arma::mat X = internalData.getData().submat(rows, regressors_);

    // arma::uvec order = target->getOrder();
    // arma::uvec H = target->getH();
    // arma::uvec censor = target->getCensor();
    // double sub;

    // arma::vec eta = X * beta;

    // eta -= arma::mean(eta);

    // RcppThread::Rcout << "Cox regression fit, beta = " << beta.t() << std::endl;

    // arma::vec scaledZ(eta);
    // arma::vec grad(arma::conv_to<arma::vec>::from(censor));
  
    // arma::vec theta = arma::exp(eta);
    // arma::vec grad_theta_weight(n, arma::fill::zeros);
    // double rs_sum = arma::accu(theta);
    // double HsumTheta, d, phi, grad_theta_weight_sum = 0;
    // int rs_size = n;
    // int m, l;

    std::string tempName = target.getName() + ".Z";
    int scaledZIdx = internalData.getColumn(internalData.getVariable(tempName));

    // int i = 0;
    // for (int j = 0; j < H.n_elem; j++) {
    //     sub = 0;
    // 	m = 0;
    // 	HsumTheta = 0.0;
    // 	for (int k = 0; k < H[j]; k++) {
    // 	    if (censor[order[i+k]]) {
    // 		m++;
    // 		HsumTheta += theta[order[i+k]];
    // 	    }
    // 	    sub += theta[order[i+k]];
    // 	}
	    
    // 	// Rcpp::Rcout << "rs_sum = " << rs_sum << std::endl
    // 	// 	    << "HsumTheta = " << HsumTheta << std::endl
    // 	// 	    << "H_j = " << H[j] << std::endl
    // 	// 	    << "m = " << m << std::endl
    // 	// 	    << "i = " << i << std::endl
    // 	// 	    << "first ii = " << i + H[j] << std::endl;

    // 	if (m > 0) { 

    // 	    // for (int k = 0; k < H[j]; k++) {
    // 	    //     for (int l = 0; l < m; l++) {
    // 	    // 	d = l / ((double) m);
    // 	    // 	phi = rs_sum - d * HsumTheta;
    // 	    // 	grad[order[i+k]] -= (1-d) * theta[order[i+k]] / phi;
    // 	    //     }
    // 	    // }

    // 	    if (sub - rs_sum > 1e-5) {
    // 		if ((H[j] + i) != n) {
    // 		    rs_sum = arma::accu(theta(order.subvec(i, n-1)));
    // 		    if (sub - rs_sum > 1e-5) {
    // 			throw std::runtime_error("Error in Cox regression loss, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
    // 		    }
    // 		} else {
    // 		    sub = rs_sum;
    // 		}
    // 	    }

    // 	    for (int l = 0; l < m; l++) {
    // 		d = l / ((double) m);
    // 		phi = rs_sum - d * HsumTheta;
    // 		grad_theta_weight_sum += (1-d) / phi;
    // 	    }
    // 	}

    // 	for (int k = 0; k < H[j]; k++) {
    // 	    grad_theta_weight[order[i+k]] = grad_theta_weight_sum;
    // 	}

    // 	i += H[j];
    // 	rs_sum -= sub;

    // 	// if (i + H[j] < n) {
    // 	//     for (int ii = i + H[j]; ii < n; ii++) {
    // 	// 	for (int l = 0; l < m; l++) {
    // 	// 	    d = l / ((double) m);
    // 	// 	    phi = rs_sum - d * HsumTheta;
    // 	// 	    grad[order[ii]] -= (1-d) * theta[order[ii]] / phi;
    // 	// 	}
    // 	//     }
    // 	// }
    // }

    // // Rcpp::Rcout << "grad_theta_weight: " << grad_theta_weight.t() << std::endl;

    // grad -= theta % grad_theta_weight;

    // scaledZ = eta + grad;
	
    for (int i = 0; i < n; i++) {
	internalData.set(i, scaledZIdx, scaledZ[i]);
    }
    
    // RcppThread::Rcout << "Eta:  " << eta.t() << std::endl;
    // RcppThread::Rcout << "Gradient:  " << grad.t() << std::endl;
    // RcppThread::Rcout << "ScaledZ:  " << scaledZ.t() << std::endl;

    this->coxRegression = CoxIRLSRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);

}

// [[Rcpp::export]]
void ScaledZTest(const Rcpp::DataFrame& df) {
    Rcpp::Rcout << "*******Start******* \n";
    DataSet data(df, 5);
    Rcpp::Rcout << "Dataset Constructed \n";
    IndTestMultiScaledZ itm(data, 0.05);
    Rcpp::Rcout << "Independence Test Constructed \n";

    Node target = itm.getVariable("Survival3");
    std::vector<Node> regressors;
//    regressors.push_back(itm.getVariable("X9"));
    regressors.push_back(itm.getVariable("Y22"));

    if (itm.variablesPerNode.count(itm.internalData.getVariable(target.getName())) < 1) {
	throw std::invalid_argument("Unrecogized node: " + target.getName());
    }

    for (const Node& var : regressors)    {
        if (itm.variablesPerNode.count(itm.internalData.getVariable(var.getName())) < 1)        {
            throw std::invalid_argument("Unrecogized node: " + var.getName());
        }
    }

    // std::ofstream logfile;
    // logfile.open("../itmcc_debug.log", std::ios_base::app);

    std::vector<Node> temp;
    std::vector<Node> _regressors;
    
    for (const Node& var : regressors) {
        temp = itm.variablesPerNode.at(itm.internalData.getVariable(var.getName()));
	if (var.isCensored()) {
	    _regressors.insert(_regressors.end(), temp.begin()+1, temp.end());	    
	} else {
	    _regressors.insert(_regressors.end(), temp.begin(), temp.end());
	}
    }

    RcppThread::Rcout << "Computing Effron Scaled Z..." << std::endl;

    std::vector<Node> tempZ(_regressors.begin()+1, _regressors.end());
    arma::uvec rows = itm.getNonMissingRows(target, _regressors.at(0), tempZ);

    for (const Node& var : _regressors) {
	RcppThread::Rcout << var.getName() << ", ";
    }
    RcppThread::Rcout << std::endl;

    // CoxRegression cox(itm.internalData);

    CoxRegressionResult result; // = cox.regress(target, _regressors);

    // Rcpp::Rcout << result << std::endl;

    result = itm.coxRegression.regress(target, _regressors);

    Rcpp::Rcout << result << std::endl;

    // int n = rows.size();
    
    // arma::uvec regressors_ = arma::uvec(_regressors.size());

    // for (int i = 0; i < regressors.size(); i++) {
    //     regressors_[i] = itm.internalData.getColumn(_regressors[i]);
    // }

    // arma::mat X = itm.internalData.getData().submat(rows, regressors_);

    // arma::mat xSub = itm.internalData.getData().cols(regressors_);
    
    // arma::mat x = arma::mat(xSub.n_rows, xSub.n_cols + 1);
    
    // for (arma::uword i = 0; i < x.n_rows; i++) {
    // 	for (arma::uword j = 0; j < x.n_cols; j++) {
    // 	    if (j == 0)	{
    // 		x(i, j) = 1;
    // 	    }
    // 	    else {
    // 		x(i, j) = xSub(i, j - 1);
    // 	    }
    // 	}
    // }

    // arma::uvec order = target->getOrder();
    // arma::uvec H = target->getH();
    // arma::uvec censor = target->getCensor();
    // double sub;

    // arma::vec beta(_regressors.size()+1, arma::fill::zeros);
    // arma::vec betaOld(_regressors.size()+1, arma::fill::ones);
    // arma::vec betaGrad(_regressors.size()+1, arma::fill::ones);
    
    // int iter = 0;

    // while (arma::norm(beta.subvec(1, _regressors.size()) - betaOld.subvec(1, _regressors.size())) / arma::norm(beta.subvec(1, _regressors.size())) > 1e-8) {
    // 	iter++;
    // 	if (iter > 100) break;

    // 	Rcpp::Rcout << "Iter: " << iter << std::endl;

    // 	betaOld = beta;

    // 	arma::vec eta = x.cols(1, _regressors.size()) * beta.subvec(1, _regressors.size());

    // 	eta -= arma::mean(eta);

    // 	RcppThread::Rcout << "Cox regression fit, beta = " << beta.t() << std::endl;

    // 	arma::vec scaledZ(eta);
    // 	arma::vec Z(eta);
    // 	arma::vec grad(arma::conv_to<arma::vec>::from(censor));
    // 	arma::vec diag_hessian(n, arma::fill::zeros);
  
    // 	arma::vec theta = arma::exp(eta);
    // 	arma::vec theta_weight(n, arma::fill::zeros);
    // 	arma::vec theta_weight2(n, arma::fill::zeros);
    // 	double rs_sum = arma::accu(theta);
    // 	double HsumTheta, d, phi, theta_weight_sum = 0, theta_weight2_sum = 0;
    // 	int rs_size = n;
    // 	int m, l;

    // 	std::string tempName = target->getName() + ".Z";
    // 	int scaledZIdx = itm.internalData.getColumn(itm.internalData.getVariable(tempName));

    // 	int i = 0;
    // 	for (int j = 0; j < H.n_elem; j++) {
    // 	    sub = 0;
    // 	    m = 0;
    // 	    HsumTheta = 0.0;
    // 	    for (int k = 0; k < H[j]; k++) {
    // 		if (censor[order[i+k]]) {
    // 		    m++;
    // 		    HsumTheta += theta[order[i+k]];
    // 		}
    // 		sub += theta[order[i+k]];
    // 	    }
	    
    // 	    // Rcpp::Rcout << "rs_sum = " << rs_sum << std::endl
    // 	    // 	    << "HsumTheta = " << HsumTheta << std::endl
    // 	    // 	    << "H_j = " << H[j] << std::endl
    // 	    // 	    << "m = " << m << std::endl
    // 	    // 	    << "i = " << i << std::endl
    // 	    // 	    << "first ii = " << i + H[j] << std::endl;

    // 	    if (m > 0) { 

    // 		// for (int k = 0; k < H[j]; k++) {
    // 		//     for (int l = 0; l < m; l++) {
    // 		// 	d = l / ((double) m);
    // 		// 	phi = rs_sum - d * HsumTheta;
    // 		// 	grad[order[i+k]] -= (1-d) * theta[order[i+k]] / phi;
    // 		//     }
    // 		// }

    // 		if (sub - rs_sum > 1e-5) {
    // 		    if ((H[j] + i) != n) {
    // 			rs_sum = arma::accu(theta(order.subvec(i, n-1)));
    // 			if (sub - rs_sum > 1e-5) {
    // 			    throw std::runtime_error("Error in Cox regression loss, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
    // 			}
    // 		    } else {
    // 			sub = rs_sum;
    // 		    }
    // 		}

    // 		for (int l = 0; l < m; l++) {
    // 		    d = l / ((double) m);
    // 		    phi = rs_sum - d * HsumTheta;
    // 		    theta_weight_sum += (1-d) / phi;
    // 		    theta_weight2_sum += std::pow((1-d) / phi, 2);
    // 		}
    // 	    }

    // 	    for (int k = 0; k < H[j]; k++) {
    // 		theta_weight[order[i+k]] = theta_weight_sum;
    // 		theta_weight2[order[i+k]] = theta_weight2_sum;
    // 	    }

    // 	    i += H[j];
    // 	    rs_sum -= sub;

    // 	    // if (i + H[j] < n) {
    // 	    //     for (int ii = i + H[j]; ii < n; ii++) {
    // 	    // 	for (int l = 0; l < m; l++) {
    // 	    // 	    d = l / ((double) m);
    // 	    // 	    phi = rs_sum - d * HsumTheta;
    // 	    // 	    grad[order[ii]] -= (1-d) * theta[order[ii]] / phi;
    // 	    // 	}
    // 	    //     }
    // 	    // }
    // 	}

    // 	// Rcpp::Rcout << "grad_theta_weight: " << grad_theta_weight.t() << std::endl;

    // 	grad -= theta % theta_weight;

    // 	diag_hessian = arma::square(theta) % theta_weight2 - theta % theta_weight;

    // 	scaledZ = eta + grad;

    // 	Z = eta - grad / diag_hessian;

    // 	for (int i = 0; i < n; i++) {
    // 	    if (diag_hessian[i]==0) Z[i] = 0;
    // 	}

    // 	// if (iter==1) {
    // 	//   double nullloss = arma::as_scalar(Z.t() * arma::diagmat(-diag_hessian) * Z) / 2;
    // 	//   Rcpp::Rcout << "Null Cox loss ~= " << nullloss << std::endl;
    // 	// }
	
    // 	// arma::vec eta = X * beta;

    // 	// RcppThread::Rcout << "Cox regression fit, beta = " << beta.t() << std::endl;

    // 	// arma::vec scaledZ(eta);
    // 	// arma::vec grad(arma::conv_to<arma::vec>::from(censor));
  
    // 	// arma::vec theta = arma::exp(eta);
    // 	// double rs_sum = 0; // arma::accu(theta);
    // 	// double HsumTheta, d, phi;
    // 	// // int rs_size = n;
    // 	// int m, l;

    // 	// std::string tempName = target->getName() + ".Z";
    // 	// int scaledZIdx = itm.internalData.getColumn(itm.internalData.getVariable(tempName));

    // 	// int i = n;
    // 	// for (int j = H.n_elem-1; j >= 0; j--) {
    // 	//     i -= H[j];
    // 	//     add = 0;
    // 	//     m = 0;
    // 	//     HsumTheta = 0.0;
    // 	//     for (int k = 0; k < H[j]; k++) {
    // 	// 	if (censor[order[i+k]]) {
    // 	// 	    m++;
    // 	// 	    HsumTheta += theta[order[i+k]];
    // 	// 	}
    // 	// 	add += theta[order[i+k]];
    // 	//     }
	    
    // 	//     rs_sum += add;

    // 	//     Rcpp::Rcout << "rs_sum = " << rs_sum << std::endl
    // 	// 		<< "HsumTheta = " << HsumTheta << std::endl
    // 	// 		<< "H_j = " << H[j] << std::endl
    // 	// 		<< "m = " << m << std::endl
    // 	// 		<< "i = " << i << std::endl
    // 	// 		<< "first ii = " << i + H[j] << std::endl;

    // 	//     if (m == 0) continue;

    // 	//     for (int k = 0; k < H[j]; k++) {
    // 	// 	for (int l = 0; l < m; l++) {
    // 	// 	    d = l / ((double) m);
    // 	// 	    phi = rs_sum - d * HsumTheta;
    // 	// 	    grad[order[i+k]] -= (1-d) * theta[order[i+k]] / phi;
    // 	// 	}
    // 	//     }

    // 	//     if (i + H[j] < n) {
    // 	// 	for (int ii = i + H[j]; ii < n; ii++) {
    // 	// 	    for (int l = 0; l < m; l++) {
    // 	// 		d = l / ((double) m);
    // 	// 		phi = rs_sum - d * HsumTheta;
    // 	// 		grad[order[ii]] -= (1-d) * theta[order[ii]] / phi;
    // 	// 	    }
    // 	// 	}
    // 	//     }
    // 	// }

    // 	// scaledZ = eta + grad;
	
    // 	for (int i = 0; i < n; i++) {
    // 	    itm.internalData.set(order[i], scaledZIdx, scaledZ[order[i]]);
    // 	}

    // 	// RcppThread::Rcout << "Eta:  " << eta.t() << std::endl;
    // 	// RcppThread::Rcout << "Gradient:  " << grad.t() << std::endl;
    // 	RcppThread::Rcout << "Weight:  " << -diag_hessian.t() << std::endl;
    // 	RcppThread::Rcout << "ScaledZ:  " << scaledZ.t() << std::endl;

    // 	RcppThread::Rcout << "Z:  " << Z.t() << std::endl;

    // 	Rcpp::Rcout << "E[grad] = " << arma::mean(grad) << std::endl;

    // 	Rcpp::Rcout << "E[sqrt(-diag_hesian)] = " << arma::mean(arma::sqrt(-1/(2*diag_hessian))) << std::endl;
	
    // 	Rcpp::Rcout << "E[-diag_hessian] = " << arma::mean(-1/(2*diag_hessian)) << std::endl;
	
    // 	Rcpp::Rcout << "Linear regression target: " << itm.internalData.getVariable(scaledZIdx)->getName() << std::endl;

    // 	itm.regression = LinearRegression(itm.internalData);
	
    // 	RegressionResult result = itm.regression.regress(itm.internalData.getVariable(scaledZIdx),
    // 						     _regressors, rows);

    // 	Rcpp::Rcout << result << std::endl;

    // 	// arma::mat xSub = itm.internalData.getData().cols(regressors_);

    // 	// arma::mat x = arma::mat(xSub.n_rows, xSub.n_cols + 1);
    
    // 	// for (arma::uword i = 0; i < x.n_rows; i++) {
    // 	//     for (arma::uword j = 0; j < x.n_cols; j++) {
    // 	// 	if (j == 0)	{
    // 	// 	    x(i, j) = 1;
    // 	// 	}
    // 	// 	else {
    // 	// 	    x(i, j) = xSub(i, j - 1);
    // 	// 	}
    // 	//     }
    // 	// }

    // 	arma::vec b = arma::solve(arma::diagmat(arma::sqrt(-diag_hessian)) * x, arma::sqrt(-diag_hessian) % Z);

    // 	// arma::vec b = arma::inv_sympd(arma::trans(x) * arma::diagmat(-diag_hessian) * x) * arma::trans(x) * arma::diagmat(-diag_hessian) * Z;

	
    // 	arma::mat xTwxInv = arma::inv_sympd(x.t() * arma::diagmat(-diag_hessian) * x);

    // 	// arma::mat M = n*arma::cov(arma::trans(x));
	
    // 	// arma::mat Mb = xTwxInv * arma::trans(x) * arma::diagmat(-diag_hessian) * M * arma::diagmat(-diag_hessian) * x * xTwxInv;

    // 	arma::vec newEta = x * b;
    // 	arma::vec resid = Z - newEta;
    // 	scaledZ = -diag_hessian % resid + newEta;

    // 	double S = arma::as_scalar(resid.t() * arma::diagmat(-diag_hessian) * resid);
	
    // 	arma::vec se = arma::sqrt(S / ((double) n - 3.0) * arma::diagvec(xTwxInv));
	
    // 	// beta = result.getCoef().subvec(1,_regressors.size());

    // 	betaGrad = x.t() * arma::diagmat(-diag_hessian) * (eta - Z) / n;

    // 	Rcpp::Rcout << "Beta: " << b.t() << std::endl;

    // 	Rcpp::Rcout << "Beta v2: " << (beta - 0.5 * betaGrad).t() << std::endl;

    // 	Rcpp::Rcout << "t-statistic: " << b.t() / se.t() << std::endl;	

    // 	Rcpp::Rcout << "BetaGrad: " << betaGrad.t() << std::endl;

    // 	// arma::vec newEta = x * b;
    // 	// arma::vec resid = Z - newEta;
    // 	// scaledZ = -diag_hessian % resid + newEta;

    // 	double loss = arma::as_scalar(resid.t() * arma::diagmat(-diag_hessian) * resid) / 2;

    // 	double nullloss = arma::as_scalar(Z.t() * arma::diagmat(-diag_hessian) * Z) / 2;
	
    // 	Rcpp::Rcout << "Null Cox loss ~= " << nullloss << std::endl;
	
    // 	Rcpp::Rcout << "Cox loss ~= " << loss << std::endl;

    // 	Rcpp::Rcout << "ResidScaledZ:  " << scaledZ.t() << std::endl;

    // 	beta = b; // eta - 0.5 * betaGrad; // .subvec(1,_regressors.size());
    // }
}

// void indTestMultiTest(const Rcpp::DataFrame& df) {
//   Rcpp::Rcout << "*******Start******* \n";
//   DataSet data(df, 5);
//   Rcpp::Rcout << "Dataset Constructed \n";
//   IndTestMultiScaledZ itm(data, 0.05);
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
