// [[Rcpp::depends(BH)]]

#include "IndTestMultiCC.hpp"
#include <boost/math/distributions/chi_squared.hpp>

#include <fstream>

IndTestMultiCC::IndTestMultiCC(DataSet& data, double alpha)
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

    std::vector<Variable*> variables = internalData.getVariables();

    for (Variable* var : variables)
    {
	// Rcpp::Rcout << var->getName() << "\n";
        std::vector<Variable*> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Variable*, std::vector<Variable*>>(var, vars));
    }

    // Rcpp::Rcout << "dataset expanded\n";
    
    this->coxRegression = CoxRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;
    this->preferLinear = false;

    // Rcpp::Rcout << "regressions created\n";
}

IndTestMultiCC::IndTestMultiCC(DataSet& data, double alpha, bool preferLinear)
{
    this->timesCalled = 0;
    this->preferLinear = preferLinear;
    this->searchVariables = data.getVariables();
    DataSet originalData(data);
    this->originalData = originalData;
    DataSet internalData(data);
    this->alpha = alpha;
    this->lastP = 0;

    std::vector<Variable *> variables = internalData.getVariables();

    for (Variable *var : variables)
    {
        //   Rcpp::Rcout << var->getName() << std::endl;
        std::vector<Variable *> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Variable *, std::vector<Variable *>>(var, vars));
    }

    this->internalData = internalData;
    this->coxRegression = CoxRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;
}

IndTestMultiCC::~IndTestMultiCC() {
    // Delete expanded variables
    // for (Variable *v : internalData.getVariables()) {
    //     if (v->getName().find("MULTINOM.") != std::string::npos) {
    //         delete v;
    //     }
    // }
}

int IndTestMultiCC::reset()
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
bool IndTestMultiCC::isIndependent(Variable *x, Variable *y, std::vector<Variable *> &z, double *pReturn)
{
    this->timesCalled++;

    if (x->isCensored()) {
	if (y->isCensored()) {
	    return true; // Assumed non-interacting at first
	}
	try {
	    return isIndependentCoxRegression(x, y, z);
	} catch (std::runtime_error& e) {
	    return false;
	}
    } else if (y->isCensored()) {
	try {
	    return isIndependentCoxRegression(y, x, z);
	} catch (std::runtime_error& e) {
	    return false;
	}
    } else if (x->isDiscrete()) {
	try {
	    return isIndependentMultinomialLogisticRegression(x, y, z, pReturn);
	}  catch (std::runtime_error& e) {
	    if (y->isDiscrete()) {
		try {
		    return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
		} catch (std::runtime_error& e) {
		    return false;
		}
	    }
	    return false;
	}
    }
    else if (y->isDiscrete()) {
	try {
	    if (preferLinear) {
		return isIndependentRegression(x, y, z, pReturn);
	    }
	    else {
		return isIndependentMultinomialLogisticRegression(y, x, z, pReturn);
	    }
	} catch (std::runtime_error& e) {
	    return false;
	}
    }
    else {
	try {
	    return isIndependentRegression(x, y, z, pReturn);
	}  catch (std::runtime_error& e) {
	    try {
		return isIndependentRegression(y, x, z, pReturn);
	    } catch (std::runtime_error& e) {
		return false;
	    }
	    return false;
	}
    }
}

std::vector<Variable *> IndTestMultiCC::expandVariable(DataSet &dataSet, Variable *var)
{
    if (var->isContinuous())
    {
        std::vector<Variable*> contList;
        contList.push_back(var);
        return contList;
    }

    if (var->isCensored())
    {
        std::vector<Variable*> censList;
        censList.push_back(var);
        return censList;
    }

    if (var->isDiscrete() && ((DiscreteVariable *)var)->getNumCategories() < 3)
    {
        std::vector<Variable*> discList;
        discList.push_back(var);
        return discList;
    }

    if (!(var->isDiscrete()))
    {
        throw std::invalid_argument("*Invalid variable type*");
    }

    std::vector<std::string> varCats = ((DiscreteVariable *)var)->getCategories();

    std::vector<Variable*> variables;
    /*********************************************************************/
    std::string temp = var->getName();
    for (auto it = varCats.begin() + 1; it != varCats.end(); it++)
    {
        DiscreteVariable* newVar = new DiscreteVariable(temp + "MULTINOM." + *it, 2);

        /*********************************************************************/

        variables.push_back(newVar);

        dataSet.addVariable(newVar);

        int newVarIndex = dataSet.getColumn(newVar);
        int numCases = dataSet.getNumRows();

        for (int l = 0; l < numCases; l++)
        {
            int dataCellIndex = dataSet.getInt(l, dataSet.getColumn(var));
            if (dataCellIndex == ((DiscreteVariable *)var)->getIndex(*it))
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

bool IndTestMultiCC::isIndependentMultinomialLogisticRegression(Variable *x, Variable *y, std::vector<Variable *> &z, double *pReturn)
{
    // std::ofstream logfile;
    // logfile.open("../test_results/debug.log", std::ios_base::app);

    if (variablesPerNode.count(internalData.getVariable(x->getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized variable: " + x->getName());
    }

    if (variablesPerNode.count(internalData.getVariable(y->getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized variable: " + y->getName());
    }

    for (Variable *varZ : z)
    {
        if (variablesPerNode.count(internalData.getVariable(varZ->getName())) < 1)
        {
            throw std::invalid_argument("Unrecogized variable: " + varZ->getName());
        }
    }

    // bool debug = (x->getName() == "X5" && y->getName() == "X6");

    arma::vec pValues;

    arma::uvec rows_ = getNonMissingRows(x, y, z);
    // logisticRegression.setRows(rows_);

    std::vector<Variable*> yzList;
    std::vector<Variable*> zList;

    std::vector<Variable*> temp = variablesPerNode.at(internalData.getVariable(y->getName()));
    yzList.insert(yzList.begin(), temp.begin(), temp.end());

    for (Variable *varZ : z) {
        temp = variablesPerNode.at(internalData.getVariable(varZ->getName()));
        yzList.insert(yzList.end(), temp.begin(), temp.end());

        // temp = variablesPerNode.at(varZ);
        zList.insert(zList.end(), temp.begin(), temp.end());
    }

    // //double[][] coeffsDep = new double[variablesPerNode.get(x).size()][];
    arma::mat coeffsNull = arma::mat(); //zList.size()+1, variablesPerNode.at(x).size());
    arma::mat coeffsDep = arma::mat();  //yzList.size()+1, variablesPerNode.at(x).size());

    /*********************************************************************/
    for (int i = 0; i < variablesPerNode.at(internalData.getVariable(x->getName())).size(); i++)
    {
        Variable *varX = variablesPerNode.at(internalData.getVariable(x->getName())).at(i);
	// LogisticRegressionResult result0;
	// LogisticRegressionResult result1;
	// try {
	LogisticRegressionResult result0 = logisticRegression.regress((DiscreteVariable *)varX,
								      zList, rows_);
	LogisticRegressionResult result1 = logisticRegression.regress((DiscreteVariable *)varX,
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
    // for (Variable* v1 : yzList)
    //     for (Variable* v2 : variablesPerNode.at(v1))
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
    // for (Variable* v1 : zList)
    //     for (Variable* v2 : variablesPerNode.at(v1))
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

    int df = variablesPerNode.at(internalData.getVariable(y->getName())).size()
	* variablesPerNode.at(internalData.getVariable(x->getName())).size();
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
    // for (Variable* zVar : z)
    //     logfile << zVar->getName() << " ";
    // logfile << "}";

    // logfile << " with p = " << p << std::endl;

    // logfile.close();

    this->lastP = p;

    // //t.println(x + " is independent of " + y + " given " + z + ": " + indep);
    return indep;
}

// This takes an inordinate amount of time. -jdramsey 20150929
arma::uvec IndTestMultiCC::getNonMissingRows(Variable *x, Variable *y, std::vector<Variable *> &z)
{

    // std::vector<int> rows = std::vec(internalData.getNumRows());
    arma::uvec rows = arma::uvec(internalData.getNumRows());
    for (arma::uword k = 0; k < rows.size(); k++)
        rows[k] = k; // use uword or use double?

    std::vector<Variable*> yzList;
    std::vector<Variable*> zList;

    std::vector<Variable*> temp = variablesPerNode.at(internalData.getVariable(y->getName()));
    yzList.insert(yzList.begin(), temp.begin(), temp.end());

    for (Variable *varZ : z) {
        temp = variablesPerNode.at(internalData.getVariable(varZ->getName()));
        yzList.insert(yzList.end(), temp.begin(), temp.end());
    }

    for (Variable *var : yzList)
	if (var->isCensored())
	    rows = arma::intersect(rows, ((CensoredVariable*)var)->getCC());

    return rows;
}

bool IndTestMultiCC::isMissing(Variable *x, int i)
{
    int j = internalData.getColumn(x);

    if (x->isDiscrete())
    {
        int v = (int)internalData.getInt(i, j);

        if (v == -99)
        {
            return true;
        }
    }

    if (x->isContinuous())
    {
        double v = internalData.getInt(i, j);

        if (std::isnan(v))
        {
            return true;
        }
    }

    return false;
}

double IndTestMultiCC::multiLL(arma::mat &coeffs, Variable *dep, std::vector<Variable *> &indep)
{

    if (dep->getName() == "??")
        throw std::invalid_argument("must have a dependent node to regress on!");
    // std::ofstream logfile;
    // logfile.open("../test_results/debug.log", std::ios_base::app);

    std::vector<Variable *> depList;
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

bool IndTestMultiCC::isIndependentRegression(Variable *x, Variable *y, std::vector<Variable *> &z, double *pReturn)
{
    if (variablesPerNode.count(internalData.getVariable(x->getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + x->getName());
    }

    if (variablesPerNode.count(internalData.getVariable(y->getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + y->getName());
    }

    for (Variable *varZ : z)
    {
        if (variablesPerNode.count(internalData.getVariable(varZ->getName())) < 1)
        {
            throw std::invalid_argument("Unrecogized node: " + varZ->getName());
        }
    }

    // std::ofstream logfile;
    // logfile.open("itm_debug.log", std::ios_base::app);

    std::vector<Variable*> regressors;
    regressors.push_back(internalData.getVariable(y->getName()));

    // logfile << regressors[0]->getName() << std::endl;

    for (Variable *varZ : z)
    {
	// logfile << varZ->getName() << std::endl;
        std::vector<Variable*> temp = variablesPerNode.at(internalData.getVariable(varZ->getName()));
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
    } catch (std::runtime_error& e) {
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
    // for (Variable* zVar : z)
    //     Rcpp::Rcout << zVar->getName() << " ";
    // Rcpp::Rcout << "}";

    // Rcpp::Rcout << " with p = " << p << std::endl;

    // if (indep)
    //     logfile << x->getName() << " _||_ " << y->getName() << " | { ";
    // else
    //     logfile << x->getName() << " _|/_ " << y->getName() << " | { ";
    // for (Variable* zVar : z)
    //     logfile << zVar->getName() << " ";
    // logfile << "}";

    // logfile << " with p = " << p << std::endl;

    // logfile.close();

    return indep;
}


bool IndTestMultiCC::isIndependentCoxRegression(Variable *x, Variable *y, std::vector<Variable *> &z, double *pReturn)
{
    if (variablesPerNode.count(internalData.getVariable(x->getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + x->getName());
    }

    if (variablesPerNode.count(internalData.getVariable(y->getName())) < 1)
    {
        throw std::invalid_argument("Unrecogized node: " + y->getName());
    }

    for (Variable *varZ : z)
    {
        if (variablesPerNode.count(internalData.getVariable(varZ->getName())) < 1)
        {
            throw std::invalid_argument("Unrecogized node: " + varZ->getName());
        }
    }

    // std::ofstream logfile;
    // logfile.open("../itmcc_debug.log", std::ios_base::app);
    
    std::vector<Variable*> zList;
    std::vector<Variable*> regressors;
    std::vector<Variable*> temp = variablesPerNode.at(internalData.getVariable(y->getName()));
    regressors.insert(regressors.end(), temp.begin(), temp.end());

    for (Variable *varZ : z)
    {
        temp = variablesPerNode.at(internalData.getVariable(varZ->getName()));
        regressors.insert(regressors.end(), temp.begin(), temp.end());
	if (y->isDiscrete()) {
	    zList.insert(zList.end(), temp.begin(), temp.end());
	}
    }

    arma::uvec rows = getNonMissingRows(x, y, z);

    double p;

    CoxRegressionResult result, result0;

    try {
	result = coxRegression.regress(((CensoredVariable*)x), regressors, rows);
	if (y->isDiscrete()) {
	    result0 = coxRegression.regress(((CensoredVariable*)x), zList, rows);
	}
    } catch (std::runtime_error& e) {
	return false;
    }

    if (y->isDiscrete()) {
	double ll = result.getLoglikelihood();
	double ll0 = result0.getLoglikelihood();
	double chisq = std::max(2*(ll - ll0), 1e-15);

	// logfile << result << std::endl;
	// logfile << result0 << std::endl;

	// logfile << "2 ( " << ll << " - " << ll0 << " ) = " << chisq << std::endl;

	int df = variablesPerNode.at(internalData.getVariable(y->getName())).size();
	boost::math::chi_squared dist(df);

        p = 1.0 - cdf(dist, chisq);

	// logfile << "p = " << p << std::endl;
    }
    p = result.getP().at(0); 

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
std::vector<std::string> IndTestMultiCC::getVariableNames()
{
    std::vector<Variable *> variables = getVariables();
    std::vector<std::string> variableNames;

    for (Variable *var1 : variables)
    {
        variableNames.push_back(var1->getName());
    }
    return variableNames;
}

Variable *IndTestMultiCC::getVariable(std::string name)
{
    for (int i = 0; i < getVariables().size(); i++)
    {
        Variable *var = getVariables().at(i);
        if (var->getName() == name)
        {
            return var;
        }
    }
    Variable *emptyVar;
    return emptyVar;
}

arma::mat IndTestMultiCC::getSubsetData(DataSet &origData, std::vector<Variable *>& varSubset) {
    arma::mat origMat = origData.getData();
    arma::uvec colIndices(varSubset.size());
    arma::uvec rowIndices(origMat.n_rows);
    // for (Variable* var : varSubset){
    for (int i = 0; i < varSubset.size(); i++) {
        Variable *var = varSubset.at(i);
        arma::uword j = origData.getColumn(var);
        colIndices[i] = j;
    }
    for (int i = 0; i < origMat.n_rows; i++) {
        rowIndices[i] = i;
    }
    
    return origMat.submat(rowIndices, colIndices);
}

// void indTestMultiTest(const Rcpp::DataFrame& df) {
//   Rcpp::Rcout << "*******Start******* \n";
//   DataSet data(df, 5);
//   Rcpp::Rcout << "Dataset Constructed \n";
//   IndTestMultiCC itm(data, 0.05);
//   Rcpp::Rcout << "Independence Test Constructed \n";
//   // std::vector<std::string> varNames = itm.getVariableNames(); // unnecessary but this would be how you would normally retrieve the data
//   Variable* x = itm.getVariable("X1");
//   Variable* y = itm.getVariable("X6");
//   Rcpp::Rcout << "FIRST X: " << x->getName() << " Continuous: " << x->isContinuous() << std::endl;
//   Rcpp::Rcout << "FIRST Y: " << y->getName() << " Continuous: " << y->isContinuous() << std::endl;
//   std::vector<Variable*> z;
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
