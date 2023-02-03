// [[Rcpp::depends(BH,RcppThread)]]

#include "GrowShrink.hpp"

GrowShrink::GrowShrink(DataSet& data, int threads) : taskQueue(MAX_QUEUE_SIZE) {
    
    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }
  
    this->searchVariables = data.getVariables();
    this->originalData = DataSet(data);
    this->internalData = DataSet(data);

    std::vector<Node> variables = internalData.getVariables();

    for (const Node& var : variables)
    {
	// Rcpp::Rcout << var->getName() << "\n";
        std::vector<Node> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
    }

    // Rcpp::Rcout << "dataset expanded\n";
    
    // this->coxRegression = CoxIRLSRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;

    std::vector<Node> emptySet = {};

    for (const Node& n : variables) {
	nullBICmap[n] = regressBIC(n, emptySet, true);
	// Rcpp::Rcout << n.getName() << " Null BIC = " << nullBICmap[n] << std::endl;
    }
    
    // Rcpp::Rcout << "regressions created\n";
}


// GrowShrink(IndependenceTest* test, int threads = -1) : taskQueue(MAX_QUEUE_SIZE) {
    
//     if (threads > 0) parallelism = threads;
//     else {
//         parallelism = std::thread::hardware_concurrency();
//         if (parallelism == 0) {
//             parallelism = 4;
//             Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
//         }
//     }

//     this->test = test;
//     this->searchVariables = test->getVariables();
//     this->originalData = test->getData();
//     this->internalData = test->getData();

//     std::vector<Node> variables = internalData.getVariables();

//     for (const Node& var : variables)
//     {
// 	// Rcpp::Rcout << var->getName() << "\n";
//         std::vector<Node> vars = expandVariable(internalData, var); // See expandVariable function below
//         variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
//     }

//     // Rcpp::Rcout << "dataset expanded\n";
    
//     // this->coxRegression = CoxIRLSRegression(internalData);
//     this->logisticRegression = LogisticRegression(internalData);
//     this->regression = LinearRegression(internalData);
//     this->verbose = false;

//     // Rcpp::Rcout << "regressions created\n";
// }


std::vector<Node> GrowShrink::expandVariable(DataSet &dataSet, const Node& var) {
    if (var.isContinuous()) {
        std::vector<Node> contList;
        contList.push_back(var);
        return contList;
    }

    // if (var.isCensored()) {
    //     std::vector<Node> censList;
    //     censList.push_back(var);
	
    // 	std::string temp = var.getName();
    // 	Node newVar = Node(new ContinuousVariable(temp + ".Z"));
    // 	censList.push_back(newVar);
    // 	dataSet.addVariable(newVar);

    // 	int newVarIndex = dataSet.getColumn(newVar);
    //     int numCases = dataSet.getNumRows();

    // 	CoxIRLSRegression coxIRLS(dataSet);

    // 	std::vector<Node> regressors = {};

    // 	CoxRegressionResult result = coxIRLS.regress(var, regressors);

    // 	arma::vec scaledZ = result.getResid();

    // 	// RcppThread::Rcout << "Scaled Z " << var.getName() << ": " << scaledZ.t() << std::endl;

    //     for (int i = 0; i < numCases; i++) {
    //         dataSet.set(i, newVarIndex, scaledZ[i]);
    //     }
	
    //     return censList;
    // }

    if (var.isDiscrete() && var.getNumCategories() < 3) {
        std::vector<Node> discList;
        discList.push_back(var);
        return discList;
    }

    if (!(var.isDiscrete())) {
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

double GrowShrink::multiLL(arma::mat &coeffs, const Node& dep, std::vector<Node> &indep) {

    if (dep.getName() == "??")
        throw std::invalid_argument("must have a dependent node to regress on!");
    
    std::vector<Node> depList;
    depList.push_back(dep);

    int i = internalData.getColumn(dep);
    arma::mat a = internalData.getData();
    arma::vec depData = a.col(i);

    // RcppThread::Rcout << depData.t() << std::endl;
    
    int N = depData.n_elem; // returns number of rows

    arma::mat indepData;
    if (indep.size() == 0)
        indepData = arma::mat(N, 1, arma::fill::ones); // filling it with ones
    else {
        indepData = getSubsetData(internalData, indep);
        indepData.insert_cols(0, arma::mat(N, 1, arma::fill::ones));
    }

    arma::mat probs = indepData * coeffs;

    // RcppThread::Rcout << "Coefficients:\n" << coeffs << std::endl;

    probs.insert_cols(0, arma::mat(indepData.n_rows, 1, arma::fill::zeros)); // reference class

    probs = arma::exp(probs);

    probs.each_row( [](arma::rowvec& r) { r /= arma::sum(r); } );

    // RcppThread::Rcout << "Probabilities\n" << probs << std::endl;

    double ll = 0;
    for (int i = 0; i < N; i++) {
        // double b = probs.row(i).max();
        arma::rowvec curRow = probs.row(i); // - b;
        // curRow = arma::exp(curRow);
        // double sum = arma::sum(curRow);
        // curRow = curRow / sum;
        ll += std::log(curRow.at((int)depData.at(i)));
    }
    if (std::isnan(ll)) {
        ll = -std::numeric_limits<double>::infinity();
    }

    // RcppThread::Rcout << "logLikelihood\n" << ll << std::endl;
    
    return ll;
}


void GrowShrink::producerGrow(const Node& target, std::list<Node>& active,
			      std::list<Node>& inactive) {

    for (const Node& y : inactive) {
	// std::vector<Node> regressors(active.begin(), active.end());

	// regressors.push_back(y);
	
	// std::sort(regressors.begin(), regressors.end());

	// std::pair<Node, std::vector<Node>> key(target, regressors);

	// if (regressors.size() < 4 && scoreHistory.count(key)) {
	//     {
	// 	std::lock_guard<std::mutex> scoreLock(scoreMutex);
	// 	scoreMap[y] = scoreHistory[key];
	//     }
	//     continue;
	// }
	
	taskQueue.push(RegressionTask(target, y, active));
	if (RcppThread::isInterrupted()) {
	    break;
	}
    }

    RegressionTask poisonPill;

    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }
}

void GrowShrink::producerShrink(const Node& target, std::list<Node>& active) {
    
    for (const Node& y : active) {

	// std::vector<Node> regressors(active.begin(), active.end());

	// regressors.erase(std::find(regressors.begin(), regressors.end(), y));
	
	// std::sort(regressors.begin(), regressors.end());

	// std::pair<Node, std::vector<Node>> key(target, regressors);

	// if (regressors.size() < 4 && scoreHistory.count(key)) {
	//     {
	// 	std::lock_guard<std::mutex> scoreLock(scoreMutex);
	// 	scoreMap[y] = scoreHistory[key];
	//     }
	//     continue;
	// }
	
	taskQueue.push(RegressionTask(target, y, active));
	if (RcppThread::isInterrupted()) {
	    break;
	}
    }

    RegressionTask poisonPill;

    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }
}


double GrowShrink::regressBIC(const Node& target, std::vector<Node>& regressors, bool history) {
    
    double n = (double) internalData.getNumRows();
    double score = 1e20;

    std::sort(regressors.begin(), regressors.end());

    std::pair<Node, std::vector<Node>> key(target, regressors);

    // RcppThread::Rcout << "Regressing " << target.getName() << " on ";
    // for (const Node& n : regressors) RcppThread::Rcout << n.getName() << " ";
    // RcppThread::Rcout << std::endl;

    if (history && scoreHistory.count(key)) {
    	// RcppThread::Rcout << "Regression in history\n";
    	return scoreHistory[key];
    }

    if (target.isContinuous()) {
	try {
	    RegressionResult result;
	    result = regression.regress(target, regressors);
	    score = n * std::log(result.getRSS() / n) + penalty * std::log(n) * (regressors.size()+1);	    
	} catch (...) {
	    score = 1e20;
	}

    } else if (target.isDiscrete()) {

	try {
	    int numCats = variablesPerNode.at(internalData.getVariable(target.getName())).size();
	    double ll;

	    if (numCats == 1) {
		LogisticRegressionResult result = logisticRegression.regress(target,
									     regressors);
		ll = -result.getLogLikelihood();
    
		
	    } else { 
		arma::mat coeffs = arma::mat();

		/*********************************************************************/
		for (int i = 0; i < numCats; i++) {
		    const Node& varX = variablesPerNode.at(internalData.getVariable(target.getName())).at(i);
		    LogisticRegressionResult result = logisticRegression.regress(varX,
										 regressors);
		    coeffs.insert_cols(i, result.getCoefs());
		}

		ll = multiLL(coeffs, target, regressors);
		
	    }

	    // Check if model is saturated
	    // if (abs(ll) < 1e-5) {
	    // 	RcppThread::Rcout << "Error: saturated multinomial logistic regression model\n";
	    // 	throw std::runtime_error("saturated multinomial logistic regression model\n");
	    // }

	    score = -2 * ll + penalty * std::log(n) * numCats * (regressors.size() + 1);
		
	} catch (...) {
	    score = 1e20;
	}

    // } else if (target.isCensored()) {
    // 	try {
    // 	    CoxRegressionResult result;
    // 	    result = coxRegression.regress(target, regressors);
    // 	    score = -2 * result.getLoglikelihood() + penalty * std::log(n) * regressors.size();	    
    // 	} catch (...) {
    // 	    score = 1e20;
    // 	}

    } else {
	throw std::invalid_argument("Unrecognized variable type");
    }

    if (history) {
    	std::lock_guard<std::mutex> historyLock(historyMutex);
    	scoreHistory[key] = score;
    }
    
    return score;
}


void GrowShrink::consumerGrow(std::unordered_map<Node, double>& scoreMap) {
    while(true) {
        RegressionTask task = taskQueue.pop();

        //If poison, return
        if (task.x.isNull() && task.y.isNull()) return;

	if (variablesPerNode.count(internalData.getVariable(task.x.getName())) < 1) {
	    throw std::invalid_argument("Unrecogized variable: " + task.x.getName());
	}

	if (variablesPerNode.count(internalData.getVariable(task.y.getName())) < 1) {
	    throw std::invalid_argument("Unrecogized variable: " + task.y.getName());
	}

	for (const Node& varZ : task.z) {
	    if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)  {
		throw std::invalid_argument("Unrecogized variable: " + varZ.getName());
	    }
	}

	std::vector<Node> regressors(variablesPerNode.at(internalData.getVariable(task.y.getName())));

	for (const Node& varZ : task.z) {
	    std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));
	    // if (varZ.isCensored()) {
	    // 	if (!task.x.isCensored()) {
	    // 	    regressors.insert(regressors.end(), temp.begin()+1, temp.end());
	    // 	}	    
	    // } else {
	    // 	regressors.insert(regressors.end(), temp.begin(), temp.end());
	    // }

	    regressors.insert(regressors.end(), temp.begin(), temp.end());
	}

	double score = regressBIC(task.x, regressors, task.z.size() < 4);

	// RcppThread::Rcout << "    Consuming " << task.y.getName()
	// 		  << ", BIC = " << score << std::endl;
	
	{
	    std::lock_guard<std::mutex> scoreLock(scoreMutex);
	    scoreMap[task.y] = score;
	}
    }
}

void GrowShrink::consumerShrink(std::unordered_map<Node, double>& scoreMap) {
    while(true) {
        RegressionTask task = taskQueue.pop();

        //If poison, return
        if (task.x.isNull() && task.y.isNull()) return;

	if (variablesPerNode.count(internalData.getVariable(task.x.getName())) < 1) {
	    throw std::invalid_argument("Unrecogized variable: " + task.x.getName());
	}

	if (variablesPerNode.count(internalData.getVariable(task.y.getName())) < 1) {
	    throw std::invalid_argument("Unrecogized variable: " + task.y.getName());
	}

	for (const Node& varZ : task.z) {
	    if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)  {
		throw std::invalid_argument("Unrecogized variable: " + varZ.getName());
	    }
	}

	std::vector<Node> zList;
	
	for (const Node& varZ : task.z) {
	    std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));
	    // if (varZ.isCensored()) {
	    // 	if (!task.x.isCensored()) {
	    // 	    zList.insert(zList.end(), temp.begin()+1, temp.end());
	    // 	}	    
	    // } else {
	    // 	zList.insert(zList.end(), temp.begin(), temp.end());
	    // }
	    zList.insert(zList.end(), temp.begin(), temp.end());
	}

	std::vector<Node> yList = variablesPerNode.at(internalData.getVariable(task.y.getName()));
	// std::vector<Node> regressors;

	// RcppThread::Rcout << "    yList : { ";
	// for (const Node& y : yList) {
	//     RcppThread::Rcout << y.getName() << " ";
	// }
	// RcppThread::Rcout << "}\n";

	// RcppThread::Rcout << "    zList : { ";
	// for (const Node& z : zList) {
	//     RcppThread::Rcout << z.getName() << " ";
	// }
	// RcppThread::Rcout << "}\n";
	
	// std::set_difference(zList.begin(), zList.end(),
	// 		    yList.begin(), yList.end(),
	// 		    std::inserter(regressors, regressors.end()));

	for (const Node& y : yList) {
	    auto it = std::remove(zList.begin(), zList.end(), y);
	    zList.erase(it, zList.end());
	}

	std::vector<Node> regressors(zList);

	// RcppThread::Rcout << "    " << task.y.getName() << " : { ";
	// for (const Node& n : regressors) {
	//     RcppThread::Rcout << n.getName() << " ";
	// }
	// RcppThread::Rcout << "}\n";

	double score = regressBIC(task.x, regressors, task.z.size() <= 4);

	// RcppThread::Rcout << "    Consuming " << task.y.getName()
	// 		  << ", BIC = " << score << std::endl;

	{
	    std::lock_guard<std::mutex> scoreLock(scoreMutex);
	    scoreMap[task.y] = score;
	}
    }
}


std::list<Node> GrowShrink::search(const Node& target, std::vector<Node>& regressors,
				   double* bicReturn) {

    if (parallelism == 1) {
	return searchSingle(target, regressors, bicReturn);
    }

    double score;
    if (verbose) RcppThread::Rcout << "Searching for Markov Boundary of " << target.getName() << "...\n";
    std::list<Node> active = grow(target, regressors, &score);
    active = shrink(target, active, score, &score);
    if (verbose) RcppThread::Rcout << "Finished. \n";

    if (scoreHistory.size() > 15000) scoreHistory.clear();

    // Rcpp::Rcout << "Score History size: " << scoreHistory.size() << std::endl;
    
    // RcppThread::Rcout << "  BIC = " << score << " : { ";
    // for (const Node& n : active) {
    // 	RcppThread::Rcout << n.getName() << " ";
    // }
    // RcppThread::Rcout << "}\n";

    // score -= nullBICmap[target];

    if (bicReturn != NULL)
	*bicReturn = score;
    
    return active;
}

std::list<Node> GrowShrink::grow(const Node& target, std::vector<Node>& regressors,
				 double* bicReturn) {
    std::vector<Node> emptySet = {};
    
    double oldScore = 1e20;
    // double curScore = regressBIC(target, emptySet, true);

    double curScore = nullBICmap[target];

    std::list<Node> active;
    std::list<Node> inactive(regressors.begin(), regressors.end());

    std::unordered_map<Node, double> scoreMap;

    for (Node n : inactive) {
	scoreMap[n] = 1e20;
    }

    // RcppThread::ThreadPool pool(parallelism);

    if (verbose) RcppThread::Rcout << "  Growing...\n";
    if (verbose) RcppThread::Rcout << "    BIC = " << curScore << " : { }\n";

    while (curScore < oldScore) {
	oldScore = curScore;

	std::vector<RcppThread::Thread> threads;
	
	threads.push_back(RcppThread::Thread( [&]() { producerGrow(target, active, inactive); } ));

	for (int i = 0; i < parallelism; i++) {
	    threads.push_back(RcppThread::Thread( [&]() { consumerGrow(scoreMap); } ));
	}

	// pool.wait();

	for (int i = 0; i < threads.size(); i++) {
	    if (threads[i].joinable()) {
		threads[i].join();
	    } else {
		Rcpp::Rcout << "#### THREAD " << i << " NOT JOINABLE ####\n";
	    }
	}

	Node bestNode;

	for (const Node& n : inactive) {
	    if (scoreMap[n] < curScore) {
		curScore = scoreMap[n];
		bestNode = n;
	    }
	}

	if (!bestNode.isNull()) {
	    auto it = std::remove(inactive.begin(), inactive.end(), bestNode);
	    inactive.erase(it, inactive.end());
	    active.push_back(bestNode);
	}

	if (verbose) {
	    RcppThread::Rcout << "    BIC = " << curScore << " : { ";
	    for (const Node& n : active) {
		RcppThread::Rcout << n.getName() << " ";
	    }
	    RcppThread::Rcout << "}\n";
	}
	// scoreMap.clear();

	// for (Node n : inactive) {
	//     scoreMap[n] = 1e20;
	// }
    }

    // pool.join();
    
    if (verbose) RcppThread::Rcout << "  Finished\n";

    if (bicReturn != NULL)
	*bicReturn = curScore;
    
    return active;
}

std::list<Node> GrowShrink::shrink(const Node& target, std::list<Node>& active,
				   double score, double* bicReturn) {
    double oldScore = 1e20;
    double curScore = score;

    std::unordered_map<Node, double> scoreMap;

    for (Node n : active) {
	scoreMap[n] = 1e20;
    }

    // RcppThread::ThreadPool pool(parallelism);

    if (verbose) RcppThread::Rcout << "  Shrinking...\n";
    
    while (curScore < oldScore) {
	oldScore = curScore;

	// auto producer = std::bind(&GrowShrink::producerShrink, this);
	// auto consumer = std::bind(&GrowShrink::consumerShrink, this);

	std::vector<RcppThread::Thread> threads;
	
	threads.push_back(RcppThread::Thread( [&]() { producerShrink(target, active); } ));

	for (int i = 0; i < parallelism; i++) {
	    threads.push_back(RcppThread::Thread( [&]() { consumerShrink(scoreMap); } ));
	}

	// pool.wait();

	for (int i = 0; i < threads.size(); i++) {
	    if (threads[i].joinable()) {
		threads[i].join();
	    } else {
		Rcpp::Rcout << "#### THREAD " << i << " NOT JOINABLE ####\n";
	    }
	}

	Node bestNode;

	for (const Node& n : active) {
	    if (scoreMap[n] < curScore) {
		curScore = scoreMap[n];
		bestNode = n;
	    }
	}

	if (!bestNode.isNull()) {
	    auto it = std::remove(active.begin(), active.end(), bestNode);
	    active.erase(it, active.end());
	}

	if (verbose) {
	    RcppThread::Rcout << "    BIC = " << curScore << " : { ";
	    for (const Node& n : active) {
		RcppThread::Rcout << n.getName() << " ";
	    }
	    RcppThread::Rcout << "}\n";
	}
	
    }

    // pool.join();
    
    if (verbose) RcppThread::Rcout << "  Finished\n";

    if (bicReturn != NULL)
	*bicReturn = curScore;
    
    return active;
}


std::list<Node> GrowShrink::searchSingle(const Node& target, std::vector<Node>& regressors,
					 double* bicReturn) {

    double score;
    if (verbose) RcppThread::Rcout << "Searching for Markov Boundary of " << target.getName() << "...\n";
    std::list<Node> active = growSingle(target, regressors, &score);
    active = shrinkSingle(target, active, score, &score);
    if (verbose) RcppThread::Rcout << "Finished. \n";

    if (scoreHistory.size() > 15000) scoreHistory.clear();

    if (bicReturn != NULL)
	*bicReturn = score;
    
    return active;
}


std::list<Node> GrowShrink::growSingle(const Node& target, std::vector<Node>& regressors,
				       double* bicReturn) {
    std::vector<Node> emptySet = {};
    
    double oldScore = 1e20;
    // double curScore = regressBIC(target, emptySet, true);

    double curScore = nullBICmap[target];

    std::list<Node> active;
    std::list<Node> inactive(regressors.begin(), regressors.end());

    std::unordered_map<Node, double> scoreMap;

    for (Node n : inactive) {
	scoreMap[n] = 1e20;
    }

    // RcppThread::ThreadPool pool(parallelism);

    if (verbose) RcppThread::Rcout << "  Growing...\n";
    if (verbose) RcppThread::Rcout << "    BIC = " << curScore << " : { }\n";

    while (curScore < oldScore) {
	oldScore = curScore;

	for (const Node& n : inactive) {

	    if (variablesPerNode.count(internalData.getVariable(target.getName())) < 1) {
		throw std::invalid_argument("Unrecogized variable: " + target.getName());
	    }

	    if (variablesPerNode.count(internalData.getVariable(n.getName())) < 1) {
		throw std::invalid_argument("Unrecogized variable: " + n.getName());
	    }

	    for (const Node& varZ : active) {
		if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)  {
		    throw std::invalid_argument("Unrecogized variable: " + varZ.getName());
		}
	    }

	    std::vector<Node> regressors(variablesPerNode.at(internalData.getVariable(n.getName())));

	    for (const Node& varZ : active) {
		std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));
		// if (varZ.isCensored()) {
		// 	if (!target.isCensored()) {
		// 	    regressors.insert(regressors.end(), temp.begin()+1, temp.end());
		// 	}	    
		// } else {
		// 	regressors.insert(regressors.end(), temp.begin(), temp.end());
		// }

		regressors.insert(regressors.end(), temp.begin(), temp.end());
	    }

	    double score = regressBIC(target, regressors, active.size() < 4);

	    // RcppThread::Rcout << "    Consuming " << n.getName()
	    // 		  << ", BIC = " << score << std::endl;
	
	
	    scoreMap[n] = score;
	
	}
	
	Node bestNode;

	for (const Node& n : inactive) {
	    if (scoreMap[n] < curScore) {
		curScore = scoreMap[n];
		bestNode = n;
	    }
	}

	if (!bestNode.isNull()) {
	    auto it = std::remove(inactive.begin(), inactive.end(), bestNode);
	    inactive.erase(it, inactive.end());
	    active.push_back(bestNode);
	}

	if (verbose) {
	    RcppThread::Rcout << "    BIC = " << curScore << " : { ";
	    for (const Node& n : active) {
		RcppThread::Rcout << n.getName() << " ";
	    }
	    RcppThread::Rcout << "}\n";
	}
	// scoreMap.clear();

	// for (Node n : inactive) {
	//     scoreMap[n] = 1e20;
	// }
    }

    // pool.join();
    
    if (verbose) RcppThread::Rcout << "  Finished\n";

    if (bicReturn != NULL)
	*bicReturn = curScore;
    
    return active;
}

std::list<Node> GrowShrink::shrinkSingle(const Node& target, std::list<Node>& active,
					 double score, double* bicReturn) {
    double oldScore = 1e20;
    double curScore = score;

    std::unordered_map<Node, double> scoreMap;

    for (Node n : active) {
	scoreMap[n] = 1e20;
    }

    // RcppThread::ThreadPool pool(parallelism);

    if (verbose) RcppThread::Rcout << "  Shrinking...\n";
    
    while (curScore < oldScore) {
	oldScore = curScore;

	// auto producer = std::bind(&GrowShrink::producerShrink, this);
	// auto consumer = std::bind(&GrowShrink::consumerShrink, this);

        for (const Node& n : active) {

	    if (variablesPerNode.count(internalData.getVariable(target.getName())) < 1) {
		throw std::invalid_argument("Unrecogized variable: " + target.getName());
	    }

	    if (variablesPerNode.count(internalData.getVariable(n.getName())) < 1) {
		throw std::invalid_argument("Unrecogized variable: " + n.getName());
	    }

	    for (const Node& varZ : active) {
		if (variablesPerNode.count(internalData.getVariable(varZ.getName())) < 1)  {
		    throw std::invalid_argument("Unrecogized variable: " + varZ.getName());
		}
	    }

	    std::vector<Node> zList;
	
	    for (const Node& varZ : active) {
		std::vector<Node> temp = variablesPerNode.at(internalData.getVariable(varZ.getName()));
		// if (varZ.isCensored()) {
		// 	if (!target.isCensored()) {
		// 	    zList.insert(zList.end(), temp.begin()+1, temp.end());
		// 	}	    
		// } else {
		// 	zList.insert(zList.end(), temp.begin(), temp.end());
		// }
		zList.insert(zList.end(), temp.begin(), temp.end());
	    }

	    std::vector<Node> yList = variablesPerNode.at(internalData.getVariable(n.getName()));
	    // std::vector<Node> regressors;

	    // RcppThread::Rcout << "    yList : { ";
	    // for (const Node& y : yList) {
	    //     RcppThread::Rcout << y.getName() << " ";
	    // }
	    // RcppThread::Rcout << "}\n";

	    // RcppThread::Rcout << "    zList : { ";
	    // for (const Node& z : zList) {
	    //     RcppThread::Rcout << z.getName() << " ";
	    // }
	    // RcppThread::Rcout << "}\n";
	
	    // std::set_difference(zList.begin(), zList.end(),
	    // 		    yList.begin(), yList.end(),
	    // 		    std::inserter(regressors, regressors.end()));

	    for (const Node& y : yList) {
		auto it = std::remove(zList.begin(), zList.end(), y);
		zList.erase(it, zList.end());
	    }

	    std::vector<Node> regressors(zList);

	    // RcppThread::Rcout << "    " << n.getName() << " : { ";
	    // for (const Node& n : regressors) {
	    //     RcppThread::Rcout << n.getName() << " ";
	    // }
	    // RcppThread::Rcout << "}\n";

	    double score = regressBIC(target, regressors, active.size() <= 4);

	    // RcppThread::Rcout << "    Consuming " << n.getName()
	    // 		  << ", BIC = " << score << std::endl;

	
	    scoreMap[n] = score;
	
	}

	Node bestNode;

	for (const Node& n : active) {
	    if (scoreMap[n] < curScore) {
		curScore = scoreMap[n];
		bestNode = n;
	    }
	}

	if (!bestNode.isNull()) {
	    auto it = std::remove(active.begin(), active.end(), bestNode);
	    active.erase(it, active.end());
	}

	if (verbose) {
	    RcppThread::Rcout << "    BIC = " << curScore << " : { ";
	    for (const Node& n : active) {
		RcppThread::Rcout << n.getName() << " ";
	    }
	    RcppThread::Rcout << "}\n";
	}
	
    }

    // pool.join();
    
    if (verbose) RcppThread::Rcout << "  Finished\n";

    if (bicReturn != NULL)
	*bicReturn = curScore;
    
    return active;
}



arma::mat GrowShrink::getSubsetData(DataSet &origData, std::vector<Node>& varSubset) {
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

/**
 * @return the list of variable varNames.
 */
std::vector<std::string> GrowShrink::getVariableNames()
{
    std::vector<Node> variables = getVariables();
    std::vector<std::string> variableNames;

    for (const Node& var1 : variables)
    {
        variableNames.push_back(var1.getName());
    }
    return variableNames;
}

Node GrowShrink::getVariable(std::string name)
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
