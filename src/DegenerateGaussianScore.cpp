// [[Rcpp::depends(RcppThread)]]

#include "DegenerateGaussianScore.hpp"
#include "RcppThread.h"

DegenerateGaussianScore::DegenerateGaussianScore(DataSet data, double penalty) {
    if (data.isCensored()) {
	throw std::invalid_argument("The Degenerate Gaussian Score is not supported for censored variables.");
    }
    
    this->searchVariables = data.getVariables();
    // this->originalData = DataSet(data);
    // DataSet internalData(data);
    this->penalty = penalty;
    this->N = data.getNumRows();
    this->logN = std::log(this->N);

    std::vector<Node> variables = data.getVariables();

    for (const Node& var : variables) {
        std::vector<Node> vars = expandVariable(data, var);
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
	if (vars.size()>1) {
	    data.removeVariable(var);
	}
    }

    variables = data.getVariables();
    
    for (int i = 0; i < variables.size(); i++) {
	indexMap.insert(std::pair<Node, int>(variables[i], i));
    }

    this->verbose = false;

    this->covMat = arma::cov(data.getData());

    // Rcpp::Rcout << "covMat:\n" << covMat << std::endl;
}

std::vector<Node> DegenerateGaussianScore::expandVariable(DataSet& dataSet, const Node& var) {
    if (var.isContinuous()) {
	std::vector<Node> contList;
	contList.push_back(var);
	return contList;
    }

    if (var.isDiscrete() && var.getNumCategories() < 3)	{
	std::vector<Node> discList;
	discList.push_back(var);
	return discList;
    }

    if (!var.isDiscrete()) {
	throw std::invalid_argument("*Invalid variable type*");
    }

    std::vector<std::string> varCats = var.getCategories();

    std::vector<Node> variables;
    /*********************************************************************/
    std::string temp = var.getName();
    for (auto it = varCats.begin() + 1; it != varCats.end(); it++) {
	const Node& newVar = Node(new DiscreteVariable(temp + "MULTINOM." + *it, 2));

	/*********************************************************************/

	variables.push_back(newVar);

	dataSet.addVariable(newVar);

	int newVarIndex = dataSet.getColumn(newVar);
	int numCases = dataSet.getNumRows();

	for (int l = 0; l < numCases; l++) {
	    int dataCellIndex = dataSet.getInt(l, dataSet.getColumn(var));
	    if (dataCellIndex == var.getIndex(*it)) {
		dataSet.set(l, newVarIndex, 1);
	    }
	    else {
		dataSet.set(l, newVarIndex, 0);
	    }
	}
    }

    return variables;
}


arma::uvec DegenerateGaussianScore::getIndices(const std::vector<Node>& nodes) {
    arma::uvec idx(nodes.size(), arma::fill::zeros);

    for (arma::uword i = 0; i < nodes.size(); i++) {
	idx(i) = (uint) indexMap[nodes.at(i)];
    }

    return idx;
}


double DegenerateGaussianScore::localScore(const Node& x, const std::vector<Node>& z) {

    if (variablesPerNode.count(x) < 1) {
        throw std::invalid_argument("Unrecogized node : " + x.getName());
    }

    for (const Node& varZ : z) {
        if (variablesPerNode.count(varZ) < 1) {
            throw std::invalid_argument("Unrecogized node : " + varZ.getName());
        }
    }

    if (z.size() == 0)
	return 0;
    
    arma::uvec xIdx(getIndices(variablesPerNode.at(x)));

    std::vector<Node> regressors;

    for (const Node& varZ : z) {
        std::vector<Node> temp(variablesPerNode.at(varZ));
        regressors.insert(regressors.end(), temp.begin(), temp.end());
    }

    arma::uvec zIdx(getIndices(regressors));

    // Rcpp::Rcout << "xIdx:\n" << xIdx.t();
    
    // Rcpp::Rcout << "zIdx:\n" << zIdx.t();
    
    arma::mat sigmae = covMat(xIdx,xIdx) - covMat(xIdx,zIdx) * arma::solve(covMat(zIdx,zIdx), covMat(zIdx,xIdx));

    // Rcpp::Rcout << "sigma epsilon:\n" << sigmae << std::endl;

    // Rcpp::Rcout << "sigma x:\n" << covMat(xIdx,xIdx) << std::endl;

    // Rcpp::Rcout << "inv(sigmax)*sigmae:\n" << arma::solve(covMat(xIdx,xIdx), sigmae) << std::endl;

    double scaledLogDet;

    if (xIdx.n_elem > 1) {
	scaledLogDet = std::real(arma::log_det(arma::solve(covMat(xIdx,xIdx), sigmae)));
    } else {
	scaledLogDet = std::log(sigmae(0)) - std::log(covMat(xIdx(0),xIdx(0)));
    }

    // Rcpp::Rcout << "scaledLogDet = " << scaledLogDet << std::endl;

    // Rcpp::Rcout << "-2 * ll = " << N * scaledLogDet << std::endl;

    // Rcpp::Rcout << "prior term = " << penalty * logN * xIdx.n_elem * zIdx.n_elem << std::endl;

    return N * scaledLogDet + penalty * logN * xIdx.n_elem * zIdx.n_elem;
}


/**
 * @return the list of variable varNames.
 */
std::vector<std::string> DegenerateGaussianScore::getVariableNames()
{
    std::vector<Node> variables = getVariables();
    std::vector<std::string> variableNames;

    for (const Node& var1 : variables)
	{
	    variableNames.push_back(var1.getName());
	}
    return variableNames;
}

Node DegenerateGaussianScore::getVariable(std::string name)
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

// no export // [[Rcpp::export]]
double DGScoreTest(const Rcpp::DataFrame& df,
		   std::string targetName,
		   std::vector<std::string>& regressorNames,
		   double penalty=1.0) {
    DataSet data = DataSet(df);
    data.dropMissing();
    
    Rcpp::Rcout << "-----START----- \n";
    Node target = data.getVariable(targetName);
    std::vector<Node> inputRegressors;
    for (std::string varName : regressorNames) {
	inputRegressors.push_back(data.getVariable(varName));
    }

    DegenerateGaussianScore dgScore(data, penalty);
    double score = dgScore.localScore(target, inputRegressors);
    Rcpp::Rcout << "DG Score = " << score << std::endl;
    Rcpp::Rcout << "-----END----- \n";

    return score;
}


// no export // [[Rcpp::export]]
void DGScoreTest2(const Rcpp::DataFrame& df) {
    DataSet data = DataSet(df);
    data.dropMissing();

    auto startTime = std::chrono::high_resolution_clock::now();

    long elapsedTime;
    
    DegenerateGaussianScore dgScore(data, 1.0);
    double score;

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    Rcpp::Rcout << "Covariance Matrix elapsed time: " << elapsedTime / 1000.0 << " s\n";
    
    Rcpp::Rcout << "-----START----- \n";
    int m = data.getNumColumns();

    for (int i = 0; i < 1000000; i++) {
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
	
	score = dgScore.localScore(target, z);
	
	// Rcpp::Rcout << ": " << score << "\n";
    }

    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-startTime).count();

    Rcpp::Rcout << "Total elapsed time: " << elapsedTime / 1000.0 << " s\n";
    Rcpp::Rcout << "Average time per score: " << elapsedTime / 1000000.0 << " ms\n";
    
    Rcpp::Rcout << "-----END----- \n";
}
