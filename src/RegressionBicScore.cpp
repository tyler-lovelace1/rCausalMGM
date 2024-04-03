// [[Rcpp::depends(RcppThread)]]

#include "RegressionBicScore.hpp"
#include "RcppThread.h"

RegressionBicScore::RegressionBicScore(DataSet& data, double penalty) {
    // if (data.isCensored()) {
    // 	throw std::invalid_argument("The Degenerate Gaussian Score is not supported for censored variables.");
    // }
    
    this->searchVariables = data.getVariables();
    this->originalData = DataSet(data);
    this->internalData = DataSet(data);
    this->penalty = penalty;
    this->N = originalData.getNumRows();
    this->logN = std::log(this->N);

    std::vector<Node> variables = internalData.getVariables();

    for (Node var : variables) {
        std::vector<Node> vars = expandVariable(internalData, var);
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
        if (var.isCensored()) {
	    double nEvents = var.getNEvents();
	    this->logNevents[var] = std::log(nEvents);
	}
    }

    this->coxRegression = CoxRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);

    this->verbose = false;

    std::vector<Node> emptySet = {};

    // Rcpp::Rcout << "Null Likelihoods:\n";

    for (Node n : variables) {
	if (n.isContinuous()) {
	    nullLL[n] = logLikLinearRegression(n, emptySet);
	} else if (n.isDiscrete()) {
	    nullLL[n] = logLikMultinomialLogisticRegression(n, emptySet);
	} else if (n.isCensored()) {
	    nullLL[n] = logLikCoxRegression(n, emptySet);
	    resetWZ(n, emptySet);
	} else {
	    throw std::invalid_argument(n.getName() + " is an unrecognized variable type");
	}
	// Rcpp::Rcout << "  " << n << " : " << nullLL.at(n) << "\n";
    }
    
}

std::vector<Node> RegressionBicScore::expandVariable(DataSet& dataSet, const Node& var) {
    if (var.isContinuous())
    {
        std::vector<Node> contList;
        contList.push_back(var);
        return contList;
    }

    if (var.isCensored()) {
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

void RegressionBicScore::resetWZ(Node target, std::vector<Node>& neighbors) {
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

    CoxRegressionResult result;

    try {
	result = coxRegression.regress(target, regressors);
    } catch (...) {
	result = coxRegression.regress(target, emptySet);
    }

    arma::vec WZ(result.getResid());

    std::vector<std::string> _neighbors;

    for (const Node& var : neighbors) {
	_neighbors.push_back(var.getName());
    }

    target.setNeighbors(_neighbors);
    target.setWZ(WZ);

    if (internalData.updateNode(target)) {
	this->coxRegression = CoxRegression(internalData);
	this->logisticRegression = LogisticRegression(internalData);
	this->regression = LinearRegression(internalData);
    }
}

arma::mat RegressionBicScore::getSubsetData(DataSet& origData, std::vector<Node>& varSubset) {
    arma::mat origMat = origData.getData();
    arma::uvec colIndices(varSubset.size());
    arma::uvec rowIndices(origMat.n_rows);
    // for (const Node& var : varSubset){
    for (int i = 0; i < varSubset.size(); i++) {
        const Node& var = varSubset.at(i);
        arma::uword j = origData.getColumn(var);
        colIndices[i] = j;
    }
    for (int i = 0; i < origMat.n_rows; i++) {
        rowIndices[i] = i;
    }
    return origMat.submat(rowIndices, colIndices);
}

double RegressionBicScore::multiLL(arma::mat &coeffs, const Node& dep, std::vector<Node>& indep) {

    if (dep.getName() == "??")
        throw std::invalid_argument("must have a dependent node to regress on!");
    
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

    probs.insert_cols(0, arma::mat(indepData.n_rows, 1, arma::fill::zeros)); // reference class

    double ll = 0;
    for (int i = 0; i < N; i++) {
        arma::rowvec curRow = probs.row(i);
	curRow -= curRow.max();
        double logsumexp = std::log(arma::sum(arma::exp(curRow)));
	curRow -= logsumexp;
	ll += curRow((int) depData.at(i));	
    }
    
    // probs = arma::exp(probs);

    // probs.each_row( [](arma::rowvec& r) { r /= arma::sum(r); } );

    // double ll = 0;
    // for (int i = 0; i < N; i++) {
    //     arma::rowvec curRow = probs.row(i); // - b;
    //     ll += std::log(curRow.at((int)depData.at(i)));
    // }

    if (std::isnan(ll)) {
        ll = -std::numeric_limits<double>::infinity();
    }
    
    return ll;
}

double RegressionBicScore::logLikLinearRegression(const Node& x, std::vector<Node>& regressors) {
    double ll = 1e20;
    try {
	RegressionResult result;
	result = regression.regress(x, regressors);
	ll = -N / 2.0 * std::log(result.getRSS() / N);
    } catch (...) {
	ll = 1e20;
    }
    return ll;
}

double RegressionBicScore::logLikMultinomialLogisticRegression(const Node& x, std::vector<Node>& regressors) {
    double ll = 1e20;
    try {
	int numCats = variablesPerNode.at(x).size();

	LogisticRegressionResult result;
	
	if (numCats == 1) {
	    result = logisticRegression.regress(x, regressors);
	    ll = -result.getLogLikelihood() / 2.0;
	} else { 
	    arma::mat coeffs = arma::mat();

	    /*********************************************************************/
	    for (int i = 0; i < numCats; i++) {
		const Node& varX = variablesPerNode.at(x).at(i);
		result = logisticRegression.regress(varX, regressors);
		coeffs.insert_cols(i, result.getCoefs());
	    }
	    
	    ll = multiLL(coeffs, x, regressors);

	    // if (numCats == 1) {
	    // 	Rcpp::Rcout << "multiLL = " << ll << std::endl;
	    // 	Rcpp::Rcout << "logisticRegression ll = " << result.getLogLikelihood() << std::endl;
	    // }
	
	}
	
	} catch (...) {
	    ll = 1e20;
	}
	return ll;
    }

double RegressionBicScore::logLikCoxRegression(const Node& x, std::vector<Node>& regressors) {
    double ll = 1e20;
    try {
	CoxRegressionResult result;
	result = coxRegression.regress(x, regressors);
	ll = result.getLoglikelihood();
    } catch (...) {
	ll = 1e20;
    }
    return ll;
}


double RegressionBicScore::localScore(const Node& x, const std::vector<Node>& z) {

    if (variablesPerNode.count(x) < 1) {
        throw std::invalid_argument("Unrecogized node: " + x.getName());
    }

    for (const Node& varZ : z) {
        if (variablesPerNode.count(varZ) < 1) {
            throw std::invalid_argument("Unrecogized node: " + varZ.getName());
        }
    }
    
    // arma::uvec xIdx(getIndices(variablesPerNode.at(x)));

    std::vector<Node> regressors;

    for (const Node& varZ : z) {
        std::vector<Node> temp(variablesPerNode.at(varZ));
        regressors.insert(regressors.end(), temp.begin(), temp.end());
    }

    double score = 1e20;
    int K;

    if (x.isContinuous()) {
	K = regressors.size() + 1;
	score = -2 * (logLikLinearRegression(x, regressors) - nullLL.at(x)) + penalty * logN * K;
    } else if (x.isDiscrete()) {
	K = (regressors.size() + 1) * variablesPerNode.at(x).size();
	score = -2 * (logLikMultinomialLogisticRegression(x, regressors) - nullLL.at(x))
	    + penalty * logN * K;
    } else if (x.isCensored()) {
	K = regressors.size();
	score = - 2 * (logLikCoxRegression(x, regressors) - nullLL.at(x))
	    + penalty * logNevents.at(x) * K;
    } else {
	throw std::invalid_argument(x.getName() + " is an unrecognized variable type");
    }

    return score;

    // arma::uvec zIdx(getIndices(regressors));

    // Rcpp::Rcout << "xIdx:\n" << xIdx.t();
    
    // Rcpp::Rcout << "zIdx:\n" << zIdx.t();
    
    // arma::mat sigmae = covMat(xIdx,xIdx) - covMat(xIdx,zIdx) * arma::solve(covMat(zIdx,zIdx), covMat(zIdx,xIdx));

    // Rcpp::Rcout << "sigma epsilon:\n" << sigmae << std::endl;

    // Rcpp::Rcout << "sigma x:\n" << covMat(xIdx,xIdx) << std::endl;

    // Rcpp::Rcout << "inv(sigmax)*sigmae:\n" << arma::solve(covMat(xIdx,xIdx), sigmae) << std::endl;

    // double scaledLogDet;

    // if (xIdx.n_elem > 1) {
    // 	scaledLogDet = std::real(arma::log_det(arma::solve(covMat(xIdx,xIdx), sigmae)));
    // } else {
    // 	scaledLogDet = std::log(sigmae(0)) - std::log(covMat(xIdx(0),xIdx(0)));
    // }

    // Rcpp::Rcout << "scaledLogDet = " << scaledLogDet << std::endl;

    // Rcpp::Rcout << "-2 * ll = " << N * scaledLogDet << std::endl;

    // Rcpp::Rcout << "prior term = " << penalty * logN * xIdx.n_elem * zIdx.n_elem << std::endl;

    // return N * ll + penalty * logN * nParams;
}


/**
 * @return the list of variable varNames.
 */
std::vector<std::string> RegressionBicScore::getVariableNames() {
    std::vector<Node> variables = getVariables();
    std::vector<std::string> variableNames;

    for (const Node& var1 : variables)
    {
        variableNames.push_back(var1.getName());
    }
    return variableNames;
}

Node RegressionBicScore::getVariable(std::string name) {
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

// no export // [[Rcpp::export]]
double RegrBicScoreTest(const Rcpp::DataFrame& df,
		   std::string targetName,
		   std::vector<std::string>& regressorNames) {
    DataSet data = DataSet(df);
    data.dropMissing();
    
    Rcpp::Rcout << "-----START----- \n";
    Node target = data.getVariable(targetName);
    std::vector<Node> inputRegressors;
    for (std::string varName : regressorNames) {
	inputRegressors.push_back(data.getVariable(varName));
    }

    RegressionBicScore regrBicScore(data, 1.0);
    double score = regrBicScore.localScore(target, inputRegressors);
    Rcpp::Rcout << "Regression BIC Score = " << score << std::endl;
    Rcpp::Rcout << "-----END----- \n";

    return score;
}


// no export // [[Rcpp::export]]
void RegrBicScoreTest2(const Rcpp::DataFrame& df) {
    DataSet data = DataSet(df);
    data.dropMissing();

    auto startTime = std::chrono::high_resolution_clock::now();

    long elapsedTime;
    
    RegressionBicScore regrBicScore(data, 1.0);
    double score;

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    Rcpp::Rcout << "Initialize elapsed time: " << elapsedTime / 1000.0 << " s\n";
    
    Rcpp::Rcout << "-----START----- \n";
    int m = data.getNumColumns();

    for (int i = 0; i < 10000; i++) {
	int targetIdx = std::floor(m * R::runif(0,1));
	Node target = data.getVariable(targetIdx);
	std::set<Node> regressors;
	for (int ii = 0; ii < 5; ii++) {
	    int regrIdx = std::floor(m * R::runif(0,1));
	    if (regrIdx != targetIdx) {
		regressors.insert(data.getVariable(regrIdx));
	    }
	}
	std::vector<Node> z(regressors.begin(), regressors.end());

	// Rcpp::Rcout << target << " | ";
	// for (const Node& zVar : z) {
	//     Rcpp::Rcout << zVar << " ";
	// }
	
	score = regrBicScore.localScore(target, z);
	
	// Rcpp::Rcout << ": " << score << "\n";
    }

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    Rcpp::Rcout << "Total elapsed time: " << elapsedTime / 1000.0 << " s\n";
    Rcpp::Rcout << "Average time per score: " << elapsedTime / 10000.0 << " ms\n";
    
    Rcpp::Rcout << "-----END----- \n";
}
