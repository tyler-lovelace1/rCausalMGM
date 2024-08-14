// [[Rcpp::depends(RcppThread)]]

#include "SearchCV.hpp"

SearchCV::SearchCV(DataSet& data, std::string alg, uint nfolds, int threads) {
    this->scoreNodes = data.getVariables();
    this->originalData = DataSet(data);
    this->internalData = DataSet(data);
    this->nfolds = nfolds;
    this->alg = alg;

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            RcppThread::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }
    
    std::vector<Node> variables = internalData.getVariables();

    for (const Node& var : variables) {
        std::vector<Node> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
    }

    int N = data.getNumRows();

    foldid = arma::linspace<arma::uvec>(1, N, N) - 1;

    foldid.transform([nfolds](arma::uword val) { return val % nfolds + 1; });

    foldid = arma::shuffle(foldid);

    this->coxRegression = CoxRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;

    for (Node n : scoreNodes) {
	if (n.isCensored()) {
	    std::vector<Node> emptySet;
	    CoxRegressionResult result = coxRegression.regress(n, emptySet);
	    arma::vec WZ(result.getResid());
	    n.setWZ(WZ);

	    if (internalData.updateNode(n)) {
		this->coxRegression = CoxRegression(internalData);
		this->logisticRegression = LogisticRegression(internalData);
		this->regression = LinearRegression(internalData);
	    }
	}
    }
}

SearchCV::SearchCV(DataSet& data, std::string alg, const arma::uvec& foldid, int threads) {

    if (!checkFoldID(foldid)) {
	throw std::invalid_argument("Invalid input for foldid. Values must be in the range 1:nfolds, with folds of approximately equal size.");
    }

    this->alg = alg;

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            RcppThread::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }
    
    this->scoreNodes = data.getVariables();
    this->originalData = DataSet(data);
    this->internalData = DataSet(data);
    
    std::vector<Node> variables = internalData.getVariables();

    for (const Node& var : variables) {
        std::vector<Node> vars = expandVariable(internalData, var); // See expandVariable function below
        variablesPerNode.insert(std::pair<Node, std::vector<Node>>(var, vars));
    }

    int N = data.getNumRows();

    this->foldid = foldid;
    this->nfolds = foldid.max();

    this->coxRegression = CoxRegression(internalData);
    this->logisticRegression = LogisticRegression(internalData);
    this->regression = LinearRegression(internalData);
    this->verbose = false;

    for (Node n : scoreNodes) {
	if (n.isCensored()) {
	    std::vector<Node> emptySet;
	    CoxRegressionResult result = coxRegression.regress(n, emptySet);
	    arma::vec WZ(result.getResid());
	    n.setWZ(WZ);

	    if (internalData.updateNode(n)) {
		this->coxRegression = CoxRegression(internalData);
		this->logisticRegression = LogisticRegression(internalData);
		this->regression = LinearRegression(internalData);
	    }
	}
    }
}

bool SearchCV::checkFoldID(const arma::uvec& foldid) {
    uint nfolds = foldid.max();
    arma::uvec foldSampSize(nfolds, arma::fill::zeros);

    for (uint k = 1; k <= nfolds; k++) {
	arma::uvec kIdx = arma::find(foldid==k);
	if (kIdx.n_elem==0) {
	    return false;
	}
	foldSampSize(k-1) = kIdx.n_elem;
    }
    if (arma::accu(foldSampSize) != foldid.n_elem) {
	return false;
    }
    if (foldSampSize.max()-foldSampSize.min() > 0.05 * foldid.n_elem) {
	return false;
    }
    return true;
}

std::vector<Node> SearchCV::expandVariable(DataSet& dataSet, const Node& var) {
    if (var.isContinuous()) {
	std::vector<Node> contList;
	contList.push_back(var);
	return contList;
    }

    if (var.isCensored()) {
	std::vector<Node> censList;	
	censList.push_back(var);
	return censList;
    }

    if (var.isDiscrete() && var.getNumCategories() < 3) {
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


double SearchCV::multiTestLL(arma::mat& coeffs, const Node& dep,
			     std::vector<Node>& indep, arma::uvec testRows) {

    // RcppThread::Rcout << dep.getName() << " multiTestLL call\n";
    
    arma::uvec depIdx(1);
    depIdx[0] = internalData.getColumn(dep);
    
    arma::uvec indepIdx(indep.size());
    for (int i = 0; i < indep.size(); i++) {
        indepIdx[i] = internalData.getColumn(indep.at(i));
    }

    // RcppThread::Rcout << "depIdx:\n" << depIdx << std::endl;
    // RcppThread::Rcout << "indepIdx:\n" << indepIdx << std::endl;

    // RcppThread::Rcout << "testRows:\n" << testRows << std::endl;

    arma::vec depData = internalData.getSubsetData(testRows, depIdx).as_col();
    
    int N = depData.n_elem; // returns number of rows

    arma::mat indepData;
    if (indep.size() == 0)
        indepData = arma::mat(N, 1, arma::fill::ones); // filling it with ones
    else {
        indepData = internalData.getSubsetData(testRows, indepIdx);
        indepData.insert_cols(0, arma::mat(N, 1, arma::fill::ones));
    }

    // RcppThread::Rcout << "coeffs:\n" << coeffs << std::endl;

    arma::mat probs = indepData * coeffs;

    probs.insert_cols(0, arma::mat(indepData.n_rows, 1, arma::fill::zeros)); // reference class

    // RcppThread::Rcout << "log(probs):\n" << probs << std::endl;
    
    double ll = 0;
    for (int i = 0; i < N; i++) {
        arma::rowvec curRow = probs.row(i);
	curRow -= curRow.max();
        double logsumexp = std::log(arma::sum(arma::exp(curRow)));
	curRow -= logsumexp;
	ll += curRow((int) depData.at(i));	
    }
    
    if (std::isnan(ll)) {
        ll = -std::numeric_limits<double>::infinity();
    }

    // RcppThread::Rcout << dep.getName() << " multiTestLL finished\n";
    
    return ll;
}

double SearchCV::scoreTestLLTask(const Node& dep, std::vector<Node>& indep, int k) {

    // RcppThread::Rcout << dep.getName() << " scoreTestLLTask call\n";

    arma::uvec trainRows = arma::find(foldid != k);
    arma::uvec testRows = arma::find(foldid == k);

    double ll = 0;
    int N = testRows.n_elem;

    // RcppThread::Rcout << "Scoring " << dep << " fold " << k << ":\n  indep: ";

    // for (int i = 0; i < indep.size(); i++) {
    // 	RcppThread::Rcout << indep[i] << " ";
    // }
    
    // RcppThread::Rcout << std::endl;
    
    // 	if (indep[i].isCensored()) {
    // 	    std::vector<Node> emptySet;
    // 	    result = coxRegression.regress(n, emptySet, testRows);
    // 	    arma::vec WZ(result.getResid());
    // 	    indep[i].setWZ(WZ);
    // 	}
    // }

    if (dep.isContinuous()) {

	RegressionResult result = regression.regress(dep, indep, trainRows);

	arma::uvec depIdx(1);
	depIdx[0] = internalData.getColumn(dep);

	arma::vec depData = internalData.getSubsetData(testRows, depIdx).as_col();

	N = depData.n_elem; // returns number of rows
    
	arma::uvec indepIdx(indep.size());
	for (int i = 0; i < indep.size(); i++) {
	    indepIdx[i] = internalData.getColumn(indep.at(i));
	}

	arma::mat indepData;
	if (indep.size() == 0)
	    indepData = arma::mat(N, 1, arma::fill::ones); // filling it with ones
	else {
	    indepData = internalData.getSubsetData(testRows, indepIdx);
	    indepData.insert_cols(0, arma::mat(N, 1, arma::fill::ones));
	}

	// for (int i = 0; i < indep.size(); i++) {
	//     if (indep[i].isCensored()) {
		
	//     }
	// }

	arma::vec beta = result.getCoef();

	double rss = regression.rss(indepData, depData, beta);

	ll = -N / 2.0 * std::log(rss / N);
	
    } else if (dep.isDiscrete()) {

	LogisticRegressionResult result;
	arma::mat coeffs = arma::mat();
	
	for (int i = 0; i < variablesPerNode.at(dep).size(); i++) {
	    const Node& varDep = variablesPerNode.at(dep).at(i);

	    result = logisticRegression.regress(varDep, indep, trainRows);

	    // RcppThread::Rcout << varDep << " coefficients:\n" << result.getCoefs() << std::endl;

	    coeffs.insert_cols(i, result.getCoefs());
	}

	// RcppThread::Rcout << "All coefficients:\n" << coeffs << std::endl;

	ll = multiTestLL(coeffs, dep, indep, testRows);

    } else if (dep.isCensored()) {

	CoxRegressionResult result = coxRegression.regress(dep, indep, trainRows);

	// RcppThread::Rcout << "Cox Regression Result:\n\n" << result << std::endl;

        arma::uvec depIdx(1);
	depIdx[0] = internalData.getColumn(dep);

	arma::vec depData = internalData.getSubsetData(testRows, depIdx).as_col();

	Node depCopy(dep);

	arma::uvec censor = depCopy.getCensorVec();
	censor = censor(testRows);
	arma::uvec strata = depCopy.getStrata();
	strata = strata(testRows);
        depCopy.setCensor(depData, censor, strata);
	
	N = depData.n_elem; // returns number of rows
	
	arma::uvec indepIdx(indep.size());
	for (int i = 0; i < indep.size(); i++) {
	    indepIdx[i] = internalData.getColumn(indep.at(i));
	}

	arma::mat indepData;
	if (indep.size() == 0)
	    indepData = arma::mat(N, 1, arma::fill::ones); // filling it with ones
	else {
	    indepData = internalData.getSubsetData(testRows, indepIdx);
	    // indepData.insert_cols(0, arma::mat(N, 1, arma::fill::ones));
	}

	arma::vec beta = result.getCoef();

	N = depCopy.getNEvents();
	
	ll = coxRegression.loss(beta, indepData, depCopy);

    } else {
	throw std::invalid_argument("Node " + dep.getName() + " is an unrecognized type.");
    }

    // RcppThread::Rcout << dep.getName() << " scoreTestLLTask finished\n";

    return ll / N;
    
}


std::vector<double> SearchCV::scoreGraphTestLL(EdgeListGraph graph, int k) {

    // RcppThread::Rcout << "scoreGraphTestLL call\n" << graph << std::endl;

    int N = scoreNodes.size();
    
    RcppThread::ThreadPool pool(std::max(1, std::min(parallelism, N)));

    std::vector<std::future<std::vector<double>>> futures(N);
    arma::vec scores(N);
    arma::vec mbSizes(N);

    auto scoreTask = [&] (int i) {
	// RcppThread::Rcout << "  Scoring Node " << scoreNodes.at(i) << std::endl;
	std::vector<Node> mb = graph.getMarkovBlanket(scoreNodes.at(i));
	std::vector<Node> regressors;
	// RcppThread::Rcout << "    MB : { ";
	for (Node n : mb) {
	    // RcppThread::Rcout << n << " ";
	    std::vector<Node> temp = variablesPerNode.at(n);
	    regressors.insert(regressors.end(), temp.begin(), temp.end());
	}
	// RcppThread::Rcout << "}\n";
	std::vector<double> result = { scoreTestLLTask(scoreNodes.at(i), regressors, k), (double) mb.size() };
	return result;
    };

    for (uint i = 0; i < N; i++) {
	// std::vector<Node> mb = graph.getMarkovBlanket(scoreNodes.at(i));
	futures[i] = pool.pushReturn(scoreTask, i);
    }
    for (uint i = 0; i < N; i++) {
	std::vector<double> scoreOutput = futures[i].get();
	scores(i) = scoreOutput[0];
	mbSizes(i) = scoreOutput[1];
	// RcppThread::Rcout << scoreNodes.at(i) << " = " << scores(i) << std::endl;
    }

    pool.join();

    // RcppThread::Rcout << "scoreGraphTestLL finished\n";

    return { -arma::accu(scores), arma::mean(mbSizes) };
}


std::vector<EdgeListGraph> SearchCV::causalCV() {
    arma::vec lambdaVec;
    if (initialGraph != NULL) {
	std::string alg = initialGraph->getAlgorithm();
	if (alg == "MGM") {
	    lambdaVec = initialGraph->getHyperParam("lambda");
	} else if (alg == "CoxMGM") {
	    lambdaVec = initialGraph->getHyperParam("lambda");
	} else {
	    throw std::invalid_argument("Unsupported initial graph type for cross-validation. If a method other than MGM is being used, cross-validation must be done externally in R.");
	}
    }

    // RcppThread::Rcout << "causalCV call\n";

    std::vector<double> lambda(lambdaVec.begin(), lambdaVec.end());

    std::vector<std::vector<arma::vec>> loglik;
    std::vector<std::vector<arma::vec>> mbSize;
    for (int orIdx = 0; orIdx < orientRules.size(); orIdx++) {
	loglik.push_back(std::vector<arma::vec>());
	mbSize.push_back(std::vector<arma::vec>());
	for (int aIdx = 0; aIdx < alphas.n_elem; aIdx++) {
	    loglik.at(orIdx).push_back(arma::vec(nfolds));
	    mbSize.at(orIdx).push_back(arma::vec(nfolds));
	}
    }

    // Rcpp::Rcout << "result vectors created\n";
    
    // arma::mat loglik(nfolds, alphas.n_elem);

    // OrientRule orientRule;

    bool censFlag = originalData.isCensored();

    for (uint k = 1; k <= nfolds; k++) {
	arma::urowvec trainIdx = arma::conv_to<arma::urowvec>::from(arma::find(foldid!=k));
	DataSet train(originalData, trainIdx);

	EdgeListGraph g, ig;

	if (verbose) RcppThread::Rcout << "  Running Fold " << k << "...\n";
	
	if (initialGraph != NULL) {
	    // MGM mgm(train, lambda);
	    // ig = mgm.search();

	    if (!censFlag) {
	      MGM mgm(train, lambda);
	      ig = mgm.search();	       
	    } else {
	      CoxMGM coxmgm(train, lambda);
	      ig = coxmgm.search();
	    }
	}

	for (int aIdx = 0; aIdx < alphas.n_elem; aIdx++) {
	    IndTestMultiCox itm(train, alphas(aIdx));

	    if (verbose) RcppThread::Rcout << "\r    Alpha = " << alphas[aIdx];

	    if (censFlag) {
		std::vector<Node> emptySet;
		for (Node n : scoreNodes) {
		    if (n.isCensored()) {
		        itm.resetWZ(n, emptySet);
		    }
		}
	    }

	    if (alg == "pc") {
		PcStable causalAlg((IndependenceTest*) &itm);

		if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		causalAlg.setKnowledge(knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(fdr);
		causalAlg.setOrientRule(orientRules.at(0));
	    
		if (initialGraph != NULL) {
		    causalAlg.setInitialGraph(&ig);
		}
	    
		g = causalAlg.search();
	    
		std::vector<double> scoreOutput = scoreGraphTestLL(g, k);

		loglik[0][aIdx][k-1] = scoreOutput[0];
		mbSize[0][aIdx][k-1] = scoreOutput[1];
	      
		for (int orIdx = 1; orIdx < orientRules.size(); orIdx++) {
		    g = causalAlg.reorientWithRule(orientRules.at(orIdx));

		    std::vector<double> scoreOutput = scoreGraphTestLL(g, k);
		
		    loglik[orIdx][aIdx][k-1] = scoreOutput[0];
		    mbSize[orIdx][aIdx][k-1] = scoreOutput[1];

		}
	    } else if (alg=="fci") {
		Fci causalAlg((IndependenceTest*) &itm);

		if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		causalAlg.setKnowledge(knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(fdr);
		causalAlg.setOrientRule(orientRules.at(0));
	    
		if (initialGraph != NULL) {
		    causalAlg.setInitialGraph(&ig);
		}
	    
		g = causalAlg.search();
	    
		std::vector<double> scoreOutput = scoreGraphTestLL(g, k);

		loglik[0][aIdx][k-1] = scoreOutput[0];
		mbSize[0][aIdx][k-1] = scoreOutput[1];
	      
		for (int orIdx = 1; orIdx < orientRules.size(); orIdx++) {
		    g = causalAlg.reorientWithRule(orientRules.at(orIdx));

		    std::vector<double> scoreOutput = scoreGraphTestLL(g, k);
		
		    loglik[orIdx][aIdx][k-1] = scoreOutput[0];
		    mbSize[orIdx][aIdx][k-1] = scoreOutput[1];

		}
	    } else {
		throw std::invalid_argument("Unknown causal algorithm " + alg + " provided. Must be pc or fci");
	    }
	}
	if (verbose) RcppThread::Rcout << std::endl;
    }

    // std::set<CvResult> results;

    for (int orIdx = 0; orIdx < orientRules.size(); orIdx++) {
	for (int aIdx = 0; aIdx < alphas.n_elem; aIdx++) {
	    results.insert(CvResult(arma::mean(loglik[orIdx][aIdx]),
				    arma::stddev(loglik[orIdx][aIdx]),
				    arma::mean(mbSize[orIdx][aIdx]),
				    alphas(aIdx),
				    orientRules.at(orIdx)));
	}
    }

    CvResult minResult(1e20, 1e20, 0);
    uint minIdx = 0, idx = 0;
    for (CvResult res : results) {
	// Rcpp::Rcout << "  MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
	if (res.mean < minResult.mean) {
	    minIdx = idx;
	    minResult = res;
	}
	idx++;
    }

    CvResult seResult(1e20, 1e20, 0);
    uint seIdx = 0;
    for (CvResult res : results) {
	if (res.mean < minResult.mean + minResult.se) {
	    seResult = res;
	    break;
	}
	seIdx++;
    }

    // idx = 0;
    // for (CvResult res : results) {
    // 	if (idx == minIdx) {
    // 	    Rcpp::Rcout << "**MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	} else if (idx == seIdx) {
    // 	    Rcpp::Rcout << "* MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	} else {
    // 	    Rcpp::Rcout << "  MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	}
    // 	idx++;
    // }

    EdgeListGraph ig, gMin, g1se;

    if (initialGraph != NULL) {
	// MGM mgm(originalData, lambda);
	// ig = mgm.search();

	if (!censFlag) {
	    MGM mgm(originalData, lambda);
	    ig = mgm.search();	       
	} else {
	    CoxMGM coxmgm(originalData, lambda);
	    ig = coxmgm.search();
	}
    }

    if (verbose) RcppThread::Rcout << "\n  Min:    MB Size:  " << minResult.mbSize << "    Mean:  " << minResult.mean << "    SE:  " << minResult.se << "    Alpha:  " << minResult.alpha << "    Orient Rule:  " << SepsetProducer::rule2str(minResult.rule) << std::endl;
    
    if (verbose) RcppThread::Rcout << "\n  1 SE:    MB Size:  " << seResult.mbSize << "    Mean:  " << seResult.mean << "    SE:  " << seResult.se << "    Alpha:  " << seResult.alpha << "    Orient Rule:  " << SepsetProducer::rule2str(seResult.rule) << std::endl;

    IndTestMultiCox itmMin(originalData, minResult.alpha);
    IndTestMultiCox itm1se(originalData, seResult.alpha);

    if (alg == "pc") {
	PcStable causalMin((IndependenceTest*) &itmMin);

	if (parallelism > 0) causalMin.setThreads(parallelism);

	causalMin.setKnowledge(knowledge);
	causalMin.setVerbose(verbose);
	causalMin.setFDR(fdr);
	causalMin.setOrientRule(minResult.rule);
    
	if (initialGraph != NULL) {
	    causalMin.setInitialGraph(&ig);
	}

	gMin = causalMin.search();

	PcStable causal1se((IndependenceTest*) &itm1se);

	if (parallelism > 0) causal1se.setThreads(parallelism);

	causal1se.setKnowledge(knowledge);
	causal1se.setVerbose(verbose);
	causal1se.setFDR(fdr);
	causal1se.setOrientRule(seResult.rule);
    
	if (initialGraph != NULL) {
	    causal1se.setInitialGraph(&ig);
	}

	g1se = causal1se.search();
    } else if (alg == "fci") {
	Fci causalMin((IndependenceTest*) &itmMin);

	if (parallelism > 0) causalMin.setThreads(parallelism);

	causalMin.setKnowledge(knowledge);
	causalMin.setVerbose(verbose);
	causalMin.setFDR(fdr);
	causalMin.setOrientRule(minResult.rule);
    
	if (initialGraph != NULL) {
	    causalMin.setInitialGraph(&ig);
	}

	gMin = causalMin.search();

	Fci causal1se((IndependenceTest*) &itm1se);

	if (parallelism > 0) causal1se.setThreads(parallelism);

	causal1se.setKnowledge(knowledge);
	causal1se.setVerbose(verbose);
	causal1se.setFDR(fdr);
	causal1se.setOrientRule(seResult.rule);
    
	if (initialGraph != NULL) {
	    causal1se.setInitialGraph(&ig);
	}

	g1se = causal1se.search();
    } else {
	throw std::invalid_argument("Unknown causal algorithm " + alg + " provided. Must be pc or fci");
    }
    
    return { gMin, g1se };
}


std::vector<EdgeListGraph> SearchCV::causalMGMGridCV() {

    // std::sort(lambdas.begin(), lambdas.end(), std::greater<double>());
    lambdas = arma::sort(lambdas, "descend");
    std::vector<double> lambda;
    
    std::vector<std::vector<std::vector<arma::vec>>> loglik;
    std::vector<std::vector<std::vector<arma::vec>>> mbSize;
    for (int lIdx = 0; lIdx < lambdas.n_elem; lIdx++) {
	loglik.push_back(std::vector<std::vector<arma::vec>>());
	mbSize.push_back(std::vector<std::vector<arma::vec>>());
	for (int orIdx = 0; orIdx < orientRules.size(); orIdx++) {
	    loglik.at(lIdx).push_back(std::vector<arma::vec>());
	    mbSize.at(lIdx).push_back(std::vector<arma::vec>());
	    for (int aIdx = 0; aIdx < alphas.n_elem; aIdx++) {
		loglik.at(lIdx).at(orIdx).push_back(arma::vec(nfolds));
		mbSize.at(lIdx).at(orIdx).push_back(arma::vec(nfolds));
	    }
	}
    }

    // RcppThread::Rcout << "causalMGMGridCV call\n";

    bool censFlag = originalData.isCensored();
    MGM mgm;
    CoxMGM coxmgm;

    for (uint k = 1; k <= nfolds; k++) {
	arma::urowvec trainIdx = arma::conv_to<arma::urowvec>::from(arma::find(foldid!=k));
	DataSet train(originalData, trainIdx);

	EdgeListGraph g, ig;

	if (verbose) RcppThread::Rcout << "  Running Fold " << k << "...\n";
	
	if (!censFlag) {
	    lambda = { lambdas[0], lambdas[0], lambdas[0] };
	    mgm = MGM(train, lambda);
	} else {
	    lambda = { lambdas[0], lambdas[0], lambdas[0], lambdas[0], lambdas[0] };
	    coxmgm = CoxMGM(train, lambda);
	    // if (verbose) RcppThread::Rcout << "  CoxMGM initialized " << k << "...\n";
	}

	for (int lIdx = 0; lIdx < lambdas.n_elem; lIdx++) {

	    if (!censFlag) {
		lambda = { lambdas[lIdx], lambdas[lIdx], lambdas[lIdx] };
	    
		mgm.setLambda(lambda);
		ig = mgm.search();
	    } else {
		lambda = { lambdas[lIdx], lambdas[lIdx], lambdas[lIdx], lambdas[lIdx], lambdas[lIdx] };
	    
		coxmgm.setLambda(lambda);
		// if (verbose) RcppThread::Rcout << "  Learning CoxMGM " << k << "...\n";
		ig = coxmgm.search();
		// if (verbose) RcppThread::Rcout << "  Finished CoxMGM " << k << "\n";
		// if (verbose) RcppThread::Rcout << ig << "\n";
	    }
	
	    for (int aIdx = 0; aIdx < alphas.n_elem; aIdx++) {
		IndTestMultiCox itm(train, alphas(aIdx));
		if (verbose) RcppThread::Rcout << "\r    Lambda = " << lambdas[lIdx] << ", Alpha = " << alphas[aIdx];
		if (alg=="pc") {
		    PcStable causalAlg((IndependenceTest*) &itm);

		    if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		    causalAlg.setKnowledge(knowledge);
		    causalAlg.setVerbose(false);
		    causalAlg.setFDR(fdr);
		    causalAlg.setOrientRule(orientRules.at(0));
	    
		    causalAlg.setInitialGraph(&ig);
		
		    g = causalAlg.search();
	    
		    std::vector<double> scoreOutput = scoreGraphTestLL(g, k);

		    loglik[lIdx][0][aIdx][k-1] = scoreOutput[0];
		    mbSize[lIdx][0][aIdx][k-1] = scoreOutput[1];

		    for (int orIdx = 1; orIdx < orientRules.size(); orIdx++) {
			g = causalAlg.reorientWithRule(orientRules.at(orIdx));

			std::vector<double> scoreOutput = scoreGraphTestLL(g, k);
		
			loglik[lIdx][orIdx][aIdx][k-1] = scoreOutput[0];
			mbSize[lIdx][orIdx][aIdx][k-1] = scoreOutput[1];
		    }
		} else if (alg=="fci") {
		    Fci causalAlg((IndependenceTest*) &itm);

		    if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		    causalAlg.setKnowledge(knowledge);
		    causalAlg.setVerbose(false);
		    causalAlg.setFDR(fdr);
		    causalAlg.setOrientRule(orientRules.at(0));
	    
		    causalAlg.setInitialGraph(&ig);

		    // if (verbose) RcppThread::Rcout << "  Learning FCI " << k << "...\n";
		
		    g = causalAlg.search();

		    // if (verbose) RcppThread::Rcout << "  Finished FCI " << k << "\n";
		    // if (verbose) RcppThread::Rcout << g << "\n";

		    // if (verbose) RcppThread::Rcout << "  Scoring graph " << k << "...\n";
		    std::vector<double> scoreOutput = scoreGraphTestLL(g, k);
		    // if (verbose) RcppThread::Rcout << "  Finished scoring graph " << k << "\n";
		    

		    loglik[lIdx][0][aIdx][k-1] = scoreOutput[0];
		    mbSize[lIdx][0][aIdx][k-1] = scoreOutput[1];

		    // if (verbose) RcppThread::Rcout << "  Reorienting graph " << k << "...\n";

		    for (int orIdx = 1; orIdx < orientRules.size(); orIdx++) {
			g = causalAlg.reorientWithRule(orientRules.at(orIdx));

			std::vector<double> scoreOutput = scoreGraphTestLL(g, k);
		
			loglik[lIdx][orIdx][aIdx][k-1] = scoreOutput[0];
			mbSize[lIdx][orIdx][aIdx][k-1] = scoreOutput[1];
		    }

		    // if (verbose) RcppThread::Rcout << "  Finished reorienting graph " << k << "\n";
		} else {
		    throw std::invalid_argument("Unknown causal algorithm " + alg + " provided. Must be pc or fci");
		}
	    }
	}
	if (verbose) RcppThread::Rcout << std::endl;
    }

    for (int lIdx = 0; lIdx < lambdas.n_elem; lIdx++) {
	lambda = { lambdas[lIdx], lambdas[lIdx], lambdas[lIdx] };
	for (int orIdx = 0; orIdx < orientRules.size(); orIdx++) {
	    for (int aIdx = 0; aIdx < alphas.n_elem; aIdx++) {
		results.insert(CvResult(arma::mean(loglik[lIdx][orIdx][aIdx]),
					arma::stddev(loglik[lIdx][orIdx][aIdx]),
					arma::mean(mbSize[lIdx][orIdx][aIdx]),
					alphas(aIdx),
					lambda,
					orientRules.at(orIdx)));
	    }
	}
    }

    CvResult minResult(1e20, 1e20, 0);
    uint minIdx = 0, idx = 0;
    for (CvResult res : results) {
	if (res.mean < minResult.mean) {
	    minIdx = idx;
	    minResult = res;
	}
	idx++;
    }

    CvResult seResult(1e20, 1e20, 0);
    uint seIdx = 0;
    for (CvResult res : results) {
	if (res.mean < minResult.mean + minResult.se) {
	    seResult = res;
	    break;
	}
	seIdx++;
    }

    // Rcpp::Rcout << std::endl;
    // idx = 0;
    // for (CvResult res : results) {
    // 	if (idx == minIdx) {
    // 	    Rcpp::Rcout << " ** MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Lambda:  " << res.lambda.at(0) << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	} else if (idx == seIdx) {
    // 	    Rcpp::Rcout << " *  MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Lambda:  " << res.lambda.at(0) << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	} else {
    // 	    Rcpp::Rcout << "    MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Lambda:  " << res.lambda.at(0) << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	}
    // 	idx++;
    // }
    // Rcpp::Rcout << std::endl;

    EdgeListGraph igMin, ig1se, gMin, g1se;

    if (verbose) RcppThread::Rcout << "\n  Min:    MB Size:  " << minResult.mbSize << "    Mean:  " << minResult.mean << "    SE:  " << minResult.se << "    Lambda:  " << minResult.lambda.at(0) << "    Alpha:  " << minResult.alpha << "    Orient Rule:  " << SepsetProducer::rule2str(minResult.rule) << std::endl;
    
    if (verbose) RcppThread::Rcout << "\n  1 SE:    MB Size:  " << seResult.mbSize << "    Mean:  " << seResult.mean << "    SE:  " << seResult.se << "    Lambda:  " << seResult.lambda.at(0) << "    Alpha:  " << seResult.alpha << "    Orient Rule:  " << SepsetProducer::rule2str(seResult.rule) << std::endl;

    if (!censFlag) {
	
	mgm = MGM(originalData, seResult.lambda);
	ig1se = mgm.search();

	mgm.setLambda(minResult.lambda);
	igMin = mgm.search();
	
    } else {
	
	coxmgm = CoxMGM(originalData, seResult.lambda);
	ig1se = coxmgm.search();

	coxmgm.setLambda(minResult.lambda);
	igMin = coxmgm.search();
    }
    
    IndTestMultiCox itmMin(originalData, minResult.alpha);
    IndTestMultiCox itm1se(originalData, seResult.alpha);

    if (alg == "pc") {
	PcStable causalMin((IndependenceTest*) &itmMin);
	
	if (parallelism > 0) causalMin.setThreads(parallelism);
	
	causalMin.setKnowledge(knowledge);
	causalMin.setVerbose(verbose);
	causalMin.setFDR(fdr);
	causalMin.setOrientRule(minResult.rule);
    
	causalMin.setInitialGraph(&igMin);

	gMin = causalMin.search();

	PcStable causal1se((IndependenceTest*) &itm1se);

	if (parallelism > 0) causal1se.setThreads(parallelism);

	causal1se.setKnowledge(knowledge);
	causal1se.setVerbose(verbose);
	causal1se.setFDR(fdr);
	causal1se.setOrientRule(seResult.rule);
    
	causal1se.setInitialGraph(&ig1se);
	
	g1se = causal1se.search();
    } else if (alg == "fci") {
	Fci causalMin((IndependenceTest*) &itmMin);
	
	if (parallelism > 0) causalMin.setThreads(parallelism);
	
	causalMin.setKnowledge(knowledge);
	causalMin.setVerbose(verbose);
	causalMin.setFDR(fdr);
	causalMin.setOrientRule(minResult.rule);
    
	causalMin.setInitialGraph(&igMin);

	gMin = causalMin.search();

	Fci causal1se((IndependenceTest*) &itm1se);

	if (parallelism > 0) causal1se.setThreads(parallelism);

	causal1se.setKnowledge(knowledge);
	causal1se.setVerbose(verbose);
	causal1se.setFDR(fdr);
	causal1se.setOrientRule(seResult.rule);
    
	causal1se.setInitialGraph(&ig1se);
	
	g1se = causal1se.search();
    } else {
	throw std::invalid_argument("Unknown causal algorithm " + alg + " provided. Must be pc or fci");
    }
    
    
    return { gMin, g1se };
}


std::vector<EdgeListGraph> SearchCV::causalMGMRandCV() {

    if (lambdas.n_elem != alphas.n_elem) {
	throw std::runtime_error("The number of lambdas and alphas must be the same for Random CV.");
    }

    // std::sort(lambdas.begin(), lambdas.end(), std::greater<double>());
    lambdas = arma::sort(lambdas, "descend");
    std::vector<double> lambda;

    int trials = lambdas.n_elem;
    
    std::vector<std::vector<arma::vec>> loglik;
    std::vector<std::vector<arma::vec>> mbSize;
    for (int orIdx = 0; orIdx < orientRules.size(); orIdx++) {
	loglik.push_back(std::vector<arma::vec>());
	mbSize.push_back(std::vector<arma::vec>());
	for (int tIdx = 0; tIdx < trials; tIdx++) {
	    loglik.at(orIdx).push_back(arma::vec(nfolds));
	    mbSize.at(orIdx).push_back(arma::vec(nfolds));
	}
    }

    // RcppThread::Rcout << "causalMGMRandCV call\n";

    bool censFlag = originalData.isCensored();
    MGM mgm;
    CoxMGM coxmgm;
    
    for (uint k = 1; k <= nfolds; k++) {
	arma::urowvec trainIdx = arma::conv_to<arma::urowvec>::from(arma::find(foldid!=k));
	DataSet train(originalData, trainIdx);

	EdgeListGraph g, ig;

	if (verbose) RcppThread::Rcout << "  Running Fold " << k << "...\n";
	
	// MGM mgm(train, lambda);

	if (!censFlag) {
	    lambda = { lambdas[0], lambdas[0], lambdas[0] };
	    mgm = MGM(train, lambda);
	} else {
	    lambda = { lambdas[0], lambdas[0], lambdas[0], lambdas[0], lambdas[0] };
	    coxmgm = CoxMGM(train, lambda);
	}


	for (int tIdx = 0; tIdx < trials; tIdx++) {

	    // lambda = { lambdas[tIdx], lambdas[tIdx], lambdas[tIdx] };
	    
	    // mgm.setLambda(lambda);
	    // ig = mgm.search();

	    if (!censFlag) {
		lambda = { lambdas[tIdx], lambdas[tIdx], lambdas[tIdx] };
	    
		mgm.setLambda(lambda);
		ig = mgm.search();
	    } else {
		lambda = { lambdas[tIdx], lambdas[tIdx], lambdas[tIdx], lambdas[tIdx], lambdas[tIdx] };
	    
		coxmgm.setLambda(lambda);
		ig = coxmgm.search();
	    }
	
	    IndTestMultiCox itm(train, alphas(tIdx));

	    if (verbose) RcppThread::Rcout << "\r    Lambda = " << lambdas[tIdx] << ", Alpha = " << alphas[tIdx];

	    if (alg=="pc") {
		PcStable causalAlg((IndependenceTest*) &itm);

		if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		causalAlg.setKnowledge(knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(fdr);
		causalAlg.setOrientRule(orientRules.at(0));
	    
		causalAlg.setInitialGraph(&ig);
		
		g = causalAlg.search();
	    
		std::vector<double> scoreOutput = scoreGraphTestLL(g, k);

		loglik[0][tIdx][k-1] = scoreOutput[0];
		mbSize[0][tIdx][k-1] = scoreOutput[1];

		for (int orIdx = 1; orIdx < orientRules.size(); orIdx++) {
		    g = causalAlg.reorientWithRule(orientRules.at(orIdx));

		    std::vector<double> scoreOutput = scoreGraphTestLL(g, k);
		
		    loglik[orIdx][tIdx][k-1] = scoreOutput[0];
		    mbSize[orIdx][tIdx][k-1] = scoreOutput[1];
		}
	    } else if (alg=="fci") {
		Fci causalAlg((IndependenceTest*) &itm);

		if (parallelism > 0) causalAlg.setThreads(parallelism);
	
		causalAlg.setKnowledge(knowledge);
		causalAlg.setVerbose(false);
		causalAlg.setFDR(fdr);
		causalAlg.setOrientRule(orientRules.at(0));
	    
		causalAlg.setInitialGraph(&ig);
		
		g = causalAlg.search();
	    
		std::vector<double> scoreOutput = scoreGraphTestLL(g, k);

		loglik[0][tIdx][k-1] = scoreOutput[0];
		mbSize[0][tIdx][k-1] = scoreOutput[1];

		for (int orIdx = 1; orIdx < orientRules.size(); orIdx++) {
		    g = causalAlg.reorientWithRule(orientRules.at(orIdx));

		    std::vector<double> scoreOutput = scoreGraphTestLL(g, k);
		
		    loglik[orIdx][tIdx][k-1] = scoreOutput[0];
		    mbSize[orIdx][tIdx][k-1] = scoreOutput[1];
		}
	    } else {
		throw std::invalid_argument("Unknown causal algorithm " + alg + " provided. Must be pc or fci");
	    }
	}

	if (verbose) RcppThread::Rcout << std::endl;
    }

    int count = 0;
    for (int tIdx = 0; tIdx < trials; tIdx++) {
	lambda = { lambdas[tIdx], lambdas[tIdx], lambdas[tIdx] };
	for (int orIdx = 0; orIdx < orientRules.size(); orIdx++) {
	    count++;
	    results.insert(CvResult(arma::mean(loglik[orIdx][tIdx]),
				    arma::stddev(loglik[orIdx][tIdx]),
				    arma::mean(mbSize[orIdx][tIdx]),
				    alphas(tIdx),
				    lambda,
				    orientRules.at(orIdx)));
	}
    }

    CvResult minResult(1e20, 1e20, 0);
    uint minIdx = 0, idx = 0;
    for (CvResult res : results) {
	if (res.mean < minResult.mean) {
	    minIdx = idx;
	    minResult = res;
	}
	idx++;
    }

    CvResult seResult(1e20, 1e20, 0);
    uint seIdx = 0;
    for (CvResult res : results) {
	if (res.mean < minResult.mean + minResult.se) {
	    seResult = res;
	    break;
	}
	seIdx++;
    }

    // Rcpp::Rcout << std::endl;
    // idx = 0;
    // for (CvResult res : results) {
    // 	if (idx == minIdx) {
    // 	    Rcpp::Rcout << " ** MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Lambda:  " << res.lambda.at(0) << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	} else if (idx == seIdx) {
    // 	    Rcpp::Rcout << " *  MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Lambda:  " << res.lambda.at(0) << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	} else {
    // 	    Rcpp::Rcout << "    MB Size:  " << res.mbSize << "    Mean:  " << res.mean << "    SE:  " << res.se << "    Lambda:  " << res.lambda.at(0) << "    Alpha:  " << res.alpha << "    Orient Rule:  " << orientStrings.at((int) res.rule) << std::endl;
    // 	}
    // 	idx++;
    // }
    // Rcpp::Rcout << std::endl;

    // Rcpp::Rcout << "    Expected total hyperparameter sets:  " << orientRules.size() * trials << std::endl;
    // Rcpp::Rcout << "    Inserted total hyperparameter sets:  " << count << std::endl;
    // Rcpp::Rcout << "    Actual total hyperparameter sets:    " << results.size() << std::endl << std::endl;

    EdgeListGraph igMin, ig1se, gMin, g1se;

    if (verbose) RcppThread::Rcout << "\n  Min:    MB Size:  " << minResult.mbSize << "    Mean:  " << minResult.mean << "    SE:  " << minResult.se << "    Lambda:  " << minResult.lambda.at(0) << "    Alpha:  " << minResult.alpha << "    Orient Rule:  " << SepsetProducer::rule2str(minResult.rule) << std::endl;
    
    if (verbose) RcppThread::Rcout << "\n  1 SE:    MB Size:  " << seResult.mbSize << "    Mean:  " << seResult.mean << "    SE:  " << seResult.se << "    Lambda:  " << seResult.lambda.at(0) << "    Alpha:  " << seResult.alpha << "    Orient Rule:  " << SepsetProducer::rule2str(seResult.rule) << std::endl;

    // MGM mgm1se(originalData, seResult.lambda);
    // ig1se = mgm1se.search();

    // MGM mgmMin(originalData, minResult.lambda);
    // igMin = mgmMin.search();

    if (!censFlag) {
	
	mgm = MGM(originalData, seResult.lambda);
	ig1se = mgm.search();

	mgm.setLambda(minResult.lambda);
	igMin = mgm.search();
	
    } else {
	
	coxmgm = CoxMGM(originalData, seResult.lambda);
	ig1se = coxmgm.search();

	coxmgm.setLambda(minResult.lambda);
	igMin = coxmgm.search();
    }

    IndTestMultiCox itmMin(originalData, minResult.alpha);
    IndTestMultiCox itm1se(originalData, seResult.alpha);

    if (alg == "pc") {
	PcStable causalMin((IndependenceTest*) &itmMin);
	
	if (parallelism > 0) causalMin.setThreads(parallelism);
	
	causalMin.setKnowledge(knowledge);
	causalMin.setVerbose(verbose);
	causalMin.setFDR(fdr);
	causalMin.setOrientRule(minResult.rule);
    
	causalMin.setInitialGraph(&igMin);

	gMin = causalMin.search();

	PcStable causal1se((IndependenceTest*) &itm1se);

	if (parallelism > 0) causal1se.setThreads(parallelism);

	causal1se.setKnowledge(knowledge);
	causal1se.setVerbose(verbose);
	causal1se.setFDR(fdr);
	causal1se.setOrientRule(seResult.rule);
    
	causal1se.setInitialGraph(&ig1se);
	
	g1se = causal1se.search();
    } else if (alg == "fci") {
	Fci causalMin((IndependenceTest*) &itmMin);
	
	if (parallelism > 0) causalMin.setThreads(parallelism);
	
	causalMin.setKnowledge(knowledge);
	causalMin.setVerbose(verbose);
	causalMin.setFDR(fdr);
	causalMin.setOrientRule(minResult.rule);
    
	causalMin.setInitialGraph(&igMin);

	gMin = causalMin.search();

	Fci causal1se((IndependenceTest*) &itm1se);

	if (parallelism > 0) causal1se.setThreads(parallelism);

	causal1se.setKnowledge(knowledge);
	causal1se.setVerbose(verbose);
	causal1se.setFDR(fdr);
	causal1se.setOrientRule(seResult.rule);
    
	causal1se.setInitialGraph(&ig1se);
	
	g1se = causal1se.search();
    } else {
	throw std::invalid_argument("Unknown causal algorithm " + alg + " provided. Must be pc or fci");
    }

    // EdgeListGraph igMin, ig1se, gMin, g1se;

    
    // if (verbose) RcppThread::Rcout << "Min:    MB Size:  " << minResult.mbSize << "    Mean:  " << minResult.mean << "    SE:  " << minResult.se << "    Lambda:  " << minResult.lambda.at(0) << "    Alpha:  " << minResult.alpha << "    Orient Rule:  " << orientStrings.at((int) minResult.rule) << std::endl;

    // MGM mgmMin(originalData, minResult.lambda);
    // igMin = mgmMin.search();

    // IndTestMultiCox itmMin(originalData, minResult.alpha);

    // PcStable pcMin((IndependenceTest*) &itmMin);

    // if (parallelism > 0) pcMin.setThreads(parallelism);

    // pcMin.setKnowledge(knowledge);
    // pcMin.setVerbose(verbose);
    // pcMin.setFDR(false);
    // pcMin.setOrientRule(minResult.rule);
    
    // pcMin.setInitialGraph(&igMin);
    
    // gMin = pcMin.search();

    // if (verbose) RcppThread::Rcout << "1 SE:    MB Size:  " << seResult.mbSize << "    Mean:  " << seResult.mean << "    SE:  " << seResult.se << "    Lambda:  " << seResult.lambda.at(0) << "    Alpha:  " << seResult.alpha << "    Orient Rule:  " << orientStrings.at((int) seResult.rule) << std::endl;

    // MGM mgm1se(originalData, seResult.lambda);
    // ig1se = mgm1se.search();

    // IndTestMultiCox itm1se(originalData, seResult.alpha);

    // PcStable pc1se((IndependenceTest*) &itm1se);

    // if (parallelism > 0) pc1se.setThreads(parallelism);

    // pc1se.setKnowledge(knowledge);
    // pc1se.setVerbose(verbose);
    // pc1se.setFDR(false);
    // pc1se.setOrientRule(seResult.rule);
    
    // pc1se.setInitialGraph(&ig1se);
    
    // g1se = pc1se.search();
    
    return { gMin, g1se };
}


// // [[Rcpp::export]]
// arma::vec scoreGraphTest1(const Rcpp::DataFrame& data, Rcpp::List graph, arma::uvec foldid) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     EdgeListGraph g = EdgeListGraph(graph, ds);

//     SearchCV cv(ds, foldid);

//     arma::vec scores(cv.nfolds);

//     for (uint k = 1; k <= cv.nfolds; k++) {
// 	std::vector<double> scoreOutput = cv.scoreGraphTestLL(g, k);
// 	scores(k-1) = scoreOutput[0];
//     }

//     return scores;
// }

// // [[Rcpp::export]]
// arma::vec scoreGraphTest2(const Rcpp::DataFrame& data, Rcpp::List graph, uint nfolds) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     EdgeListGraph g = EdgeListGraph(graph, ds);

//     SearchCV cv(ds, nfolds);

//     arma::vec scores(cv.nfolds);

//     for (uint k = 1; k <= cv.nfolds; k++) {
//         std::vector<double> scoreOutput = cv.scoreGraphTestLL(g, k);
// 	scores(k-1) = scoreOutput[0];
//     }

//     return scores;
// }

// // [[Rcpp::export]]
// Rcpp::List pcCvTest1(const Rcpp::DataFrame& data,
// 		     arma::vec alphas, arma::uvec foldid,
// 		     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     SearchCV cv(ds, foldid);
//     cv.verbose = true;
//     cv.alphas = alphas;
//     cv.orientRules = { ORIENT_MAJORITY, ORIENT_MAXP, ORIENT_CONSERVATIVE };

//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         cv.setInitialGraph(&ig);
//     }

//     std::vector<EdgeListGraph> cvOut = cv.pcStableCV();

//     int nValues = cv.results.size();
//     arma::vec mbSize(nValues);
//     arma::vec mean(nValues);
//     arma::vec se(nValues);

//     uint idx = 0;
//     for (auto res : cv.results) {
// 	mbSize(idx) = res.mbSize;
// 	mean(idx) = res.mean;
// 	se(idx) = res.se;
// 	idx++;
//     }

//     return Rcpp::List::create(Rcpp::_["graph.min"] = cvOut[0].toList(),
// 			      Rcpp::_["graph.1se"] = cvOut[1].toList(),
// 			      Rcpp::_["mbSize"] = mbSize,
// 			      Rcpp::_["mean"] = mean,
// 			      Rcpp::_["se"] = se);
// }


// // [[Rcpp::export]]
// Rcpp::List pcCvTest2(const Rcpp::DataFrame& data,
// 		     arma::vec alphas, uint nfolds,
// 		     Rcpp::Nullable<Rcpp::List> initialGraph = R_NilValue) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     SearchCV cv(ds, nfolds);
//     cv.verbose = true;
//     cv.alphas = alphas;
//     cv.orientRules = { ORIENT_MAJORITY, ORIENT_MAXP, ORIENT_CONSERVATIVE };

//     EdgeListGraph ig;
//     if (!initialGraph.isNull()) {
//         Rcpp::List _initialGraph(initialGraph);
//         ig = EdgeListGraph(_initialGraph, ds);
//         cv.setInitialGraph(&ig);
//     }

//     std::vector<EdgeListGraph> cvOut = cv.pcStableCV();

//     int nValues = cv.results.size();
//     arma::vec mbSize(nValues);
//     arma::vec mean(nValues);
//     arma::vec se(nValues);

//     uint idx = 0;
//     for (auto res : cv.results) {
// 	mbSize(idx) = res.mbSize;
// 	mean(idx) = res.mean;
// 	se(idx) = res.se;
// 	idx++;
//     }

//     return Rcpp::List::create(Rcpp::_["graph.min"] = cvOut[0].toList(),
// 			      Rcpp::_["graph.1se"] = cvOut[1].toList(),
// 			      Rcpp::_["mbSize"] = mbSize,
// 			      Rcpp::_["mean"] = mean,
// 			      Rcpp::_["se"] = se);
// }


// // [[Rcpp::export]]
// Rcpp::List pcGridCvTest2(const Rcpp::DataFrame& data, arma::vec lambdas,
// 			 arma::vec alphas, uint nfolds) {
//     DataSet ds = DataSet(data);
//     ds.dropMissing();

//     SearchCV cv(ds, nfolds);
//     cv.verbose = true;
//     cv.alphas = alphas;
//     cv.lambdas = lambdas;
//     cv.orientRules = { ORIENT_MAJORITY, ORIENT_MAXP, ORIENT_CONSERVATIVE };

//     std::vector<EdgeListGraph> cvOut = cv.pcStableGridCV();

//     int nValues = cv.results.size();
//     arma::vec mbSize(nValues);
//     arma::vec mean(nValues);
//     arma::vec se(nValues);

//     uint idx = 0;
//     for (auto res : cv.results) {
// 	mbSize(idx) = res.mbSize;
// 	mean(idx) = res.mean;
// 	se(idx) = res.se;
// 	idx++;
//     }

//     return Rcpp::List::create(Rcpp::_["graph.min"] = cvOut[0].toList(),
// 			      Rcpp::_["graph.1se"] = cvOut[1].toList(),
// 			      Rcpp::_["mbSize"] = mbSize,
// 			      Rcpp::_["mean"] = mean,
// 			      Rcpp::_["se"] = se,
// 			      Rcpp::_["foldid"] = cv.foldid);
// }


// no export // [[Rcpp::export]]
Rcpp::List causalMGMRandCvTest(const Rcpp::DataFrame& data,
			       std::string alg,
			       arma::vec lambdas,
			       arma::vec alphas,
			       uint nfolds,
			       bool fdr) {
    DataSet ds = DataSet(data);
    ds.dropMissing();

    SearchCV cv(ds, alg, nfolds);
    cv.setVerbose(true);
    cv.setAlphas(alphas);
    cv.setLambdas(lambdas);
    cv.setOrientRules({ ORIENT_MAJORITY, ORIENT_MAXP, ORIENT_CONSERVATIVE });

    std::vector<EdgeListGraph> cvOut = cv.causalMGMRandCV();

    std::set<CvResult> results = cv.getCVResults();

    int nValues = results.size();
    arma::vec mbSize(nValues);
    arma::vec mean(nValues);
    arma::vec se(nValues);

    uint idx = 0;
    for (auto res : results) {
	mbSize(idx) = res.mbSize;
	mean(idx) = res.mean;
	se(idx) = res.se;
	idx++;
    }

    return Rcpp::List::create(Rcpp::_["graph.min"] = cvOut[0].toList(),
			      Rcpp::_["graph.1se"] = cvOut[1].toList(),
			      Rcpp::_["mbSize"] = mbSize,
			      Rcpp::_["mean"] = mean,
			      Rcpp::_["se"] = se,
			      Rcpp::_["foldid"] = cv.foldid);
}
