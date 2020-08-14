#include "IndTestMulti.hpp"
#include <boost/math/distributions/chi_squared.hpp>


IndTestMulti::IndTestMulti(DataSet& data, double alpha){
  Rcpp::Rcout << "Step 0 \n";
  this->timesCalled = 0;
  this->searchVariables = data.getVariables();
  this->originalData = DataSet(data);
  this->internalData = DataSet(data);
  this->alpha = alpha;
  this->lastP = 0;
  Rcpp::Rcout << "Step 1 \n";

  std::vector<Variable*> variables = internalData.getVariables();

  for (Variable* var : variables) {
      std::vector<Variable*> vars = expandVariable(internalData, var); // See expandVariable function below
      variablesPerNode.insert(std::pair<Variable*, std::vector<Variable*>> (var, vars));
  }
  Rcpp::Rcout << "Step 3 \n";
  Rcpp::Rcout << "DUDEK TEST 0" << std::endl;

  this->logisticRegression = LogisticRegression(internalData);
  Rcpp::Rcout << "DUDEK TEST 1" << std::endl;
  this->regression = LinearRegression(internalData);
  this->verbose = false;
  this->preferLinear = false;
  Rcpp::Rcout << "Step 4 \n";
}


IndTestMulti::IndTestMulti(DataSet& data, double alpha, bool preferLinear) {
  Rcpp::Rcout << "Step 0 \n";

  this->timesCalled = 0;
  this->preferLinear = preferLinear;
  this->searchVariables = data.getVariables();
  DataSet originalData(data);
  this->originalData = originalData;
  DataSet internalData(data);
  this->alpha = alpha;
  this->lastP = 0;
  Rcpp::Rcout << "Step 1 \n";


  std::vector<Variable*> variables = internalData.getVariables();

  for (Variable* var : variables) {
      Rcpp::Rcout << var->getName() << std::endl;
      std::vector<Variable*> vars = expandVariable(internalData, var); // See expandVariable function below
      variablesPerNode.insert(std::pair<Variable*, std::vector<Variable*>> (var, vars));
  }
  Rcpp::Rcout << "Step 3 \n";

  this->internalData = internalData;
  LogisticRegression logReg (internalData);
  this->logisticRegression = logReg;
  LinearRegression linReg (internalData);
  this->regression = linReg;
  this->verbose = false;
  Rcpp::Rcout << "Step 4 \n";

}

int IndTestMulti::reset() {
  this->timesCalled = 0;
  return timesCalled;
}

/**
 * @return true if the given independence question is judged true, false if not. The independence question is of the
 * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are searchVariables in the list returned by
 * getVariableNames().
 */
bool IndTestMulti::isIndependent(Variable* x, Variable* y, std::vector<Variable*>& z) {
        this->timesCalled++;

        if (x->isDiscrete()) {
            return isIndependentMultinomialLogisticRegression(x, y, z);
        } else if (y->isDiscrete()) {
            if(preferLinear)
            {
                return isIndependentRegression(x,y,z);
            }
            else {
                return isIndependentMultinomialLogisticRegression(y, x, z);
            }
        } else {
            return isIndependentRegression(x, y, z);
        }
}

std::vector<Variable*> IndTestMulti::expandVariable(DataSet& dataSet, Variable* var) {
    if (var->isContinuous()) {
      std::vector<Variable*> contList;
      contList.push_back(var);
      return contList;
    }

    if (var->isDiscrete() && ((DiscreteVariable*) var)->getNumCategories() < 3) {
      std::vector<Variable*> discList;
      discList.push_back(var);
      return discList;
    }

    if (!(var->isDiscrete())) {
      throw std::invalid_argument("*Invalid variable type*");
    }

    std::vector<std::string> varCats = ((DiscreteVariable*) var)->getCategories();


    std::vector<Variable*> variables;
    /*********************************************************************/
    std::string temp = var->getName();
    for (auto it = varCats.begin()+1; it != varCats.end(); it++) {
        DiscreteVariable* newVar = new DiscreteVariable(temp + "MULTINOM." + *it, 2);

    /*********************************************************************/

        variables.push_back(newVar);

        dataSet.addVariable(newVar);

        int newVarIndex = dataSet.getColumn(newVar);
        int numCases = dataSet.getNumRows();

        for (int l = 0; l < numCases; l++) {
            int dataCellIndex = dataSet.getInt(l, dataSet.getColumn(var));
            if (dataCellIndex == ((DiscreteVariable*) var)->getIndex(*it)) {
	        dataSet.set(l, newVarIndex, 1);
	    }
            else {
		dataSet.set(l, newVarIndex, 0);
	    }
        }
    }

    return variables;
}

bool IndTestMulti::isIndependentMultinomialLogisticRegression(Variable* x, Variable* y, std::vector<Variable*>& z) {
    if (variablesPerNode.count(x) < 1) {
      throw std::invalid_argument("Unrecogized variable: " + x->getName());
    }

    if (variablesPerNode.count(y) < 1) {
      throw std::invalid_argument("Unrecogized variable: " + y->getName());
    }

    for (Variable* varZ : z) {
        if (variablesPerNode.count(varZ) < 1) {
          throw std::invalid_argument("Unrecogized variable: " + varZ->getName());
        }
    }

    arma::vec pValues;

    arma::uvec rows_ = getNonMissingRows(x, y, z);
    logisticRegression.setRows(rows_);

    std::vector<Variable*> yzList;
    std::vector<Variable*> zList;

    std::vector<Variable*> temp = variablesPerNode.at(y);
    yzList.insert(yzList.begin(), temp.begin(), temp.end());

    for (Variable* varZ : z) {
        temp = variablesPerNode.at(varZ);
        yzList.insert(yzList.end(), temp.begin(), temp.end());

        temp = variablesPerNode.at(varZ);
        zList.insert(zList.end(), temp.begin(), temp.end());
    }

    // //double[][] coeffsDep = new double[variablesPerNode.get(x).size()][];
    arma::mat coeffsNull = arma::mat(zList.size()+1, variablesPerNode.at(x).size());
    arma::mat coeffsDep = arma::mat(yzList.size()+1, variablesPerNode.at(x).size());

    /*********************************************************************/
    for (int i = 0; i < variablesPerNode.at(x).size(); i++) {
        Variable* varX = variablesPerNode.at(x).at(i);
        LogisticRegressionResult* result0 = logisticRegression.regress((DiscreteVariable*)varX, zList);
        LogisticRegressionResult* result1 = logisticRegression.regress((DiscreteVariable*)varX, yzList);

        coeffsNull.insert_cols(i, result0->getCoefs());
        coeffsDep.insert_cols(i, result1->getCoefs());
    }
    /*********************************************************************/

    double chisq = 2*(multiLL(coeffsDep, x, yzList) - multiLL(coeffsNull, x, zList)); // Need to make multiLL
    int df = variablesPerNode.at(y).size()*variablesPerNode.at(x).size();
    boost::math::chi_squared dist(df);
    Rcpp::Rcout << "CDF CALL 1" << std::endl;
    Rcpp::Rcout << "dist.df = " << dist.degrees_of_freedom() << std::endl;
    Rcpp::Rcout << "chisq = " << chisq << std::endl;
    double p = 1.0 - cdf(dist, chisq);

    // // double p = 1.0;
    //
    // // Choose the minimum of the p-values
    // // This is only one method that can be used, this requires every coefficient to be significant
    // //for (double val : pValues) {
    // //    if (val < p) p = val;
    // //}

    bool indep = p > alpha;

    this->lastP = p;

    // //t.println(x + " is independent of " + y + " given " + z + ": " + indep);
    return indep;
}

    // This takes an inordinate amount of time. -jdramsey 20150929
arma::uvec IndTestMulti::getNonMissingRows(Variable* x, Variable* y, std::vector<Variable*>& z) {


    // std::vector<int> rows = std::vec(internalData.getNumRows());
    arma::uvec rows = arma::uvec(internalData.getNumRows());
    for (arma::uword k = 0; k < rows.size(); k++) rows[k] = k; // use uword or use double?

    return rows;
}

bool IndTestMulti::isMissing(Variable* x, int i) {
    int j = internalData.getColumn(x);

    if (x->isDiscrete()) {
        int v = (int) internalData.getInt(i, j);

        if (v == -99) {
            return true;
        }
    }

    if (x->isContinuous()) {
        double v = internalData.getInt(i, j);

        if (std::isnan(v)) {
            return true;
        }
    }

    return false;
}


double IndTestMulti::multiLL(arma::mat& coeffs, Variable* dep, std::vector<Variable*>& indep){

    if(dep->getName() == "??") throw std::invalid_argument("must have a dependent node to regress on!");
    std::vector<Variable*> depList;
    depList.push_back(dep);

    int i = internalData.getColumn(dep);
    arma::mat a = internalData.getData();
    arma::vec depData = a.col(i);

    int N = depData.n_elem; // returns number of rows

    arma::mat indepData;
    if(indep.size()==0)
        indepData = arma::mat(N,1,arma::fill::ones); // filling it with ones
    else {
        indepData = getSubsetData(internalData, indep);
        indepData.insert_cols(0, arma::mat(N,1,arma::fill::ones));
    }

    arma::mat probs = indepData * coeffs;

    probs.insert_cols(0,arma::mat(indepData.n_rows, 1,arma::fill::ones));
    probs = arma::exp(probs);
    double ll = 0;
    for(int i = 0; i < N; i++){
        arma::vec curRow = probs.row(i);
        curRow = curRow/(arma::sum(curRow));
        ll += std::log(curRow.at((int)depData.at(i)));
    }
    return ll;
}

bool IndTestMulti::isIndependentRegression(Variable* x, Variable* y, std::vector<Variable*>& z) {    
    if (variablesPerNode.count(x) < 1) {
        throw std::invalid_argument("Unrecogized node: " + x->getName());
    }

    if (variablesPerNode.count(y) < 1) {
        throw std::invalid_argument("Unrecogized node: " + y->getName());
    }

    for (Variable* varZ : z) {
        if (variablesPerNode.count(varZ) < 1) {
            throw std::invalid_argument("Unrecogized node: " + varZ->getName());
        }
    }

    std::vector<Variable*> regressors;
    regressors.push_back(y);

    for (Variable* varZ : z) {
        std::vector<Variable*> temp = variablesPerNode.at(varZ);
        regressors.insert(regressors.end(), temp.begin(), temp.end());
    }

    arma::uvec rows = getNonMissingRows(x, y, z);
    LinearRegression regression = this->regression;
    regression.setRows(rows); // regression was never declared

    RegressionResult* result;

    try {
        result = regression.regress(x, regressors);
    } catch (std::exception e) {
        return false;
    }

    double p = result->getP().at(1); // double check on .at(1)
    this->lastP = p;

    bool indep = p > alpha;

    return indep;
}



/**
 * @return the list of variable varNames.
 */
std::vector<std::string> IndTestMulti::getVariableNames() {
    std::vector<Variable*> variables = getVariables();
    std::vector<std::string> variableNames;

    for (Variable* var1 : variables) {
        variableNames.push_back(var1->getName());
    }
    return variableNames;
}

Variable* IndTestMulti::getVariable(std::string name) {
    for (int i = 0; i < getVariables().size(); i++) {
        Variable* var = getVariables().at(i);
        if (var->getName() == name) {
            return var;
        }
    }
    Variable* emptyVar;
    return emptyVar;
}

arma::mat IndTestMulti::getSubsetData(DataSet& origData, std::vector<Variable*>& varSubset) {
  arma::mat origMat = origData.getData();
  arma::uvec colIndices (varSubset.size());
  arma::uvec rowIndices (origMat.n_rows);
  // for (Variable* var : varSubset){
  for (int i = 0; i < varSubset.size(); i++){
    Variable* var = varSubset.at(i);
    arma::uword j = origData.getColumn(var);
    colIndices[i] = j;
  }
  for (int i = 0; i < origMat.n_rows; i++){
    rowIndices[i] = i;
  }
  return origMat.submat(rowIndices, colIndices);
}

// [[Rcpp::export]]
void indTestMultiTest(const Rcpp::DataFrame& df) {
  Rcpp::Rcout << "*******Start******* \n";
  DataSet data(df, 5);
  Rcpp::Rcout << "Dataset Constructed \n";
  IndTestMulti itm(data, 0.05);
  Rcpp::Rcout << "Independence Test Constructed \n";
  // std::vector<std::string> varNames = itm.getVariableNames(); // unnecessary but this would be how you would normally retrieve the data
  Variable* x = itm.getVariable("X1");
  Variable* y = itm.getVariable("X6");
  Rcpp::Rcout << "FIRST X: " << x->getName() << " Continuous: " << x->isContinuous() << std::endl;
  Rcpp::Rcout << "FIRST Y: " << y->getName() << " Continuous: " << y->isContinuous() << std::endl;
  std::vector<Variable*> z;
  Rcpp::Rcout << "FIRST CASE: " << itm.isIndependent(x,y,z) << std::endl;
  Rcpp::Rcout << "First P: " << itm.getPValue() << std::endl;
  z.push_back(itm.getVariable("X9"));
  Rcpp::Rcout << "SECOND CASE: " << itm.isIndependent(x,y,z) << std::endl;
  Rcpp::Rcout << "Second P: " << itm.getPValue() << std::endl;
  z.push_back(itm.getVariable("X10"));
  z.push_back(itm.getVariable("X8"));
  Rcpp::Rcout << "THIRD CASE: " << itm.isIndependent(x,y,z) << std::endl;
  Rcpp::Rcout << "Third P: " << itm.getPValue() << std::endl;
  x = itm.getVariable("X2");
  Rcpp::Rcout << "X IS DISCRETE CASE: " << itm.isIndependent(x,y,z) << std::endl;
  Rcpp::Rcout << "DISCRETE CASE P: " << itm.getPValue() << std::endl;

  // try the case where there are values for x and y but z is empty, or when z has 1 or more variables
  //when x and y are mixed between    at least one where x is discrete
}
