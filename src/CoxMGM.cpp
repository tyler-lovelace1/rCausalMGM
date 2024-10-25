#include "CoxMGM.hpp"


// CoxMGM::CoxMGM(arma::mat& x, arma::mat& y, std::vector<Node>& variables, std::vector<int>& l, std::vector<double>& lambda) {
    
//     if (l.size() != y.n_cols)
//         throw std::invalid_argument("length of l doesn't match number of variables in Y");

//     if (y.n_rows != x.n_rows)
//         throw std::invalid_argument("different number of samples for x and y");

//     //lambda should have 3 values corresponding to cc, cd, and dd
//     if (lambda.size() != 3)
//         throw std::invalid_argument("Lambda should have three values for cc, cd, and dd edges respectively");
    
//     this->xDat = x;
//     this->yDat = y;
//     this->l = l;
//     this->p = x.n_cols;
//     this->q = y.n_cols;
//     this->n = x.n_rows;
//     this->variables = variables;
//     this->initVariables = variables;

//     this->lambda = arma::vec(lambda);
//     fixData();
//     initParameters();
//     calcWeights();
//     makeDummy();

// }

CoxMGM::CoxMGM(DataSet& ds) {

    bool mixed = true;

    if (!ds.isCensored()) {
	throw std::runtime_error("Cannot run CoxMGM: The dataset does not contain censored variables.");
    }

    if (ds.isContinuous()) {
	dummyVar = Node(new DiscreteVariable("dummy.gLpkx1Hs6x", 2));
	ds.addVariable(dummyVar);
	arma::uword j = ds.getColumn(dummyVar);
	for (arma::uword i = 0; i < ds.getNumRows(); i++) {
	    ds.set(i, j, std::floor(R::runif(0,2)));
	}
	mixed = false;
	qDummy = 1;
    }

    if (ds.isDiscrete()) {
	dummyVar = Node(new ContinuousVariable("dummy.qCm6jaC1VK"));
	ds.addVariable(dummyVar);
	arma::uword j = ds.getColumn(dummyVar);
	for (arma::uword i = 0; i < ds.getNumRows(); i++) {
	    ds.set(i, j, std::floor(R::rnorm(0,1)));
	}
	mixed = false;
	pDummy = 1;
    }
    

    this->xDat = ds.getContinuousData();
    this->yDat = ds.getDiscreteData();
    this->cDat = ds.getCensoredData();
    this->l = ds.getDiscLevels();
    this->p = xDat.n_cols;
    this->q = yDat.n_cols;
    this->r = cDat.n_cols;
    this->n = xDat.n_rows;

    this->zDat = arma::mat(n, r, arma::fill::zeros);
    this->wzDat = arma::mat(n, r, arma::fill::zeros);
    this->coxgrad = arma::mat(n, r, arma::fill::zeros);
    this->diagHess = arma::mat(n, r, arma::fill::zeros);
    // this->orderMat = arma::umat(n, r, arma::fill::zeros);
    this->censMat = arma::umat(n, r, arma::fill::zeros);
    this->numStrata = arma::uvec(r, arma::fill::ones);

    this->fitWeight = arma::rowvec(r, arma::fill::zeros);

    this->orderList = std::vector<std::vector<arma::uvec>>(r);
    this->censList = std::vector<std::vector<arma::uvec>>(r);
    this->idxList = std::vector<std::vector<arma::uvec>>(r);
    this->HList = std::vector<std::vector<arma::uvec>>(r);

    //the variables are now ordered continuous first then discrete
    std::vector<Node> cVar = ds.getContinuousVariables();
    std::vector<Node> dVar = ds.getDiscreteVariables();
    std::vector<Node> sVar = ds.getCensoredVariables();
    this->variables = std::vector<Node>();
    this->variables.reserve(p+q+r);
    this->variables.insert(this->variables.end(), cVar.begin(), cVar.end());
    this->variables.insert(this->variables.end(), dVar.begin(), dVar.end());
    this->variables.insert(this->variables.end(), sVar.begin(), sVar.end());

    for (int m = 0; m < r; m++) {
	this->censMat.col(m) = sVar[m].getCensorVec();
	this->orderList[m] = sVar[m].getOrder();
	this->censList[m] = sVar[m].getCensor();
	this->HList[m] = sVar[m].getH();
	this->idxList[m] = sVar[m].getIndex();
	this->numStrata(m) = sVar[m].getNumStrata();
    }
    
    this->initVariables = ds.getVariables();
    this->lambda = 5 * std::sqrt(std::log10(p+q+r)/((double)n)) * arma::vec(5, arma::fill::ones);

    // this->lambda = arma::pow(4/3.0 * this->lambda, 1.5);

    //Data is checked for 0 or 1 indexing and for missing levels and N(0,1) Standardizes continuous data
    fixData();

    // this->resid = arma::mat(xDat);
    // this->catResid = arma::mat(yDat);
    // arma::rowvec pHat = arma::mean(yDat, 0);
    // this->catResid.each_row([&pHat] (arma::rowvec& r) { r = pHat-r; });

    //Initialize all parameters to zeros
    initParameters();

    //Sets continuous variable weights to standard deviation and discrete variable weights to p*(1-p) for each category
    calcWeights();

    //Creates dummy variables for each category of discrete variables (stored in dDat)
    makeDummy();

    arma::vec dDatSum = arma::sum(dDat, 0).t();

    params.alpha2 = arma::log(dDatSum + 1) - arma::log(arma::vec(lsum, arma::fill::value(n)) - dDatSum + 1);
}


CoxMGM::CoxMGM(DataSet& ds, std::vector<double>& lambda) {

    bool mixed = true;

    if (!ds.isCensored()) {
	throw std::runtime_error("Cannot run CoxMGM: The dataset does not contain censored variables.");
    }

    if (ds.isContinuous()) {
	dummyVar = Node(new DiscreteVariable("dummy.gLpkx1Hs6x", 2));
	ds.addVariable(dummyVar);
	arma::uword j = ds.getColumn(dummyVar);
	for (arma::uword i = 0; i < ds.getNumRows(); i++) {
	    ds.set(i, j, std::floor(R::runif(0,2)));
	}
	mixed = false;
	qDummy = 1;
    }

    if (ds.isDiscrete()) {
	dummyVar = Node(new ContinuousVariable("dummy.qCm6jaC1VK"));
	ds.addVariable(dummyVar);
	arma::uword j = ds.getColumn(dummyVar);
	for (arma::uword i = 0; i < ds.getNumRows(); i++) {
	    ds.set(i, j, std::floor(R::rnorm(0,1)));
	}
	mixed = false;
	pDummy = 1;
    }
    

    this->xDat = ds.getContinuousData();
    this->yDat = ds.getDiscreteData();
    this->cDat = ds.getCensoredData();
    this->l = ds.getDiscLevels();
    this->p = xDat.n_cols;
    this->q = yDat.n_cols;
    this->r = cDat.n_cols;
    this->n = xDat.n_rows;

    this->zDat = arma::mat(n, r, arma::fill::zeros);
    this->wzDat = arma::mat(n, r, arma::fill::zeros);
    this->coxgrad = arma::mat(n, r, arma::fill::zeros);
    this->diagHess = arma::mat(n, r, arma::fill::zeros);
    // this->orderMat = arma::umat(n, r, arma::fill::zeros);
    this->censMat = arma::umat(n, r, arma::fill::zeros);
    this->numStrata = arma::uvec(r, arma::fill::ones);

    this->fitWeight = arma::rowvec(r, arma::fill::zeros);

    // this->HList = std::vector<arma::uvec>(r);
    this->orderList = std::vector<std::vector<arma::uvec>>(r);
    this->censList = std::vector<std::vector<arma::uvec>>(r);
    this->idxList = std::vector<std::vector<arma::uvec>>(r);
    this->HList = std::vector<std::vector<arma::uvec>>(r);

    //the variables are now ordered continuous first then discrete
    std::vector<Node> cVar = ds.getContinuousVariables();
    std::vector<Node> dVar = ds.getDiscreteVariables();
    std::vector<Node> sVar = ds.getCensoredVariables();
    this->variables = std::vector<Node>();
    this->variables.reserve(p+q+r);
    this->variables.insert(this->variables.end(), cVar.begin(), cVar.end());
    this->variables.insert(this->variables.end(), dVar.begin(), dVar.end());
    this->variables.insert(this->variables.end(), sVar.begin(), sVar.end());

    // for (int m = 0; m < r; m++) {
    // 	this->orderMat.col(m) = sVar[m].getOrder();
    // 	this->censMat.col(m) = sVar[m].getCensor();
    // 	this->HList[m] = sVar[m].getH();
    // }

    for (int m = 0; m < r; m++) {
	this->censMat.col(m) = sVar[m].getCensorVec();
	this->orderList[m] = sVar[m].getOrder();
	this->censList[m] = sVar[m].getCensor();
	this->HList[m] = sVar[m].getH();
	this->idxList[m] = sVar[m].getIndex();
	this->numStrata(m) = sVar[m].getNumStrata();
    }
    
    this->initVariables = ds.getVariables();
    this->lambda = arma::vec(lambda);
    // this->lambda = arma::pow(4/3.0 * this->lambda, 1.5);
    
    //Data is checked for 0 or 1 indexing and for missing levels and N(0,1) Standardizes continuous data
    fixData();

    // this->resid = arma::mat(xDat);
    // this->catResid = arma::mat(yDat);
    // arma::rowvec pHat = arma::mean(yDat, 0);
    // this->catResid.each_row([&pHat] (arma::rowvec& r) { r = pHat-r; });

    //Initialize all parameters to zeros
    initParameters();

    //Sets continuous variable weights to standard deviation and discrete variable weights to p*(1-p) for each category
    calcWeights();

    //Creates dummy variables for each category of discrete variables (stored in dDat)
    makeDummy();

    arma::vec dDatSum = arma::sum(dDat, 0).t();
    
    params.alpha2 = arma::log(dDatSum + 1) - arma::log(arma::vec(lsum, arma::fill::value(n)) - dDatSum + 1);
}

// init all parameters to zeros except for betad which is set to 1s
void CoxMGM::initParameters() {
    lcumsum = std::vector<int>(l.size()+1);
    lcumsum[0] = 0;
    for(int i = 0; i < l.size(); i++){
        lcumsum[i+1] = lcumsum[i] + l[i];
    }
    lsum = lcumsum[l.size()];

    // arma::mat beta((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros);

    params = CoxMGMParams(
        arma::mat((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros),  // beta
        arma::vec((int) xDat.n_cols,                    arma::fill::ones),   // betad
        arma::mat(lsum,              (int) xDat.n_cols, arma::fill::zeros),  // theta
        arma::mat(lsum,              lsum,              arma::fill::zeros),  // phi
	arma::mat((int) xDat.n_cols, (int) zDat.n_cols, arma::fill::zeros),  // gamma
	arma::mat(lsum,              (int) zDat.n_cols, arma::fill::zeros),  // psi
        arma::vec((int) xDat.n_cols,                    arma::fill::zeros),  // alpha1
        arma::vec(lsum,                                 arma::fill::zeros),  // alpha2
	arma::vec((int) zDat.n_cols,                    arma::fill::zeros)   // alpha3
    );
}

// avoid underflow in log(sum(exp(x))) calculation
double CoxMGM::logsumexp(arma::vec x) {
    // arma::vec myX = arma::vec(x);
    double maxX = x.max();
    x -= maxX;
    return std::log(arma::sum(arma::exp(x))) + maxX;
}

// avoid underflow in log(sum(exp(x))) calculation
arma::vec CoxMGM::logsumexp(arma::mat x) {
    arma::vec maxX = arma::max(x, 1);
    x.each_col() -= maxX;
    return arma::log(arma::sum(arma::exp(x), 1)) + maxX;
}


//calculate parameter weights as in Lee and Hastie
void CoxMGM::calcWeights() {
    // Rcpp::Rcout << "calcWeights called\n";
    weights = arma::vec(p+q+r, arma::fill::ones);

    //Continuous variable weights are standard deviations
    for (arma::uword i = 0; i < p; i++) {
        weights(i) = arma::stddev(xDat.col(i));
    }

    //Discrete variable weights for each variable-category pair are p(1-p) where p is the percentage of times that category appears
    for (arma::uword j = 0; j < q; j++) {
        double curWeight = 0;
        for (int k = 0; k < l[j]; k++) {
            arma::vec equalityVec = arma::vec(yDat.col(j)).transform( [k](double val) { return val == k+1 ? 1 : 0; } );
            double curp = arma::sum(equalityVec) / (double) n;
            curWeight += curp * (1-curp);
        }
        weights(p+j) = std::sqrt(curWeight);
    }

    arma::mat eta = arma::zeros(n,r);
    // eta.each_col( [](arma::vec& c) { c -= arma::mean(c); } );
    
    arma::mat coxgrad(arma::size(eta), arma::fill::zeros);
    arma::mat diagHess(arma::size(eta), arma::fill::zeros);
    coxgrad.zeros();
    diagHess.zeros();

    oldCoxloss = coxGradHess(eta, coxgrad, diagHess);
    nullCoxloss = oldCoxloss;

    coxWeights = -diagHess;
    diagHess.replace(0, -1e-10);
    zDat = eta - coxgrad / diagHess;

    arma::rowvec eventRate = arma::mean(arma::conv_to<arma::mat>::from(censMat), 0);
    weights.subvec(p+q, p+q+r-1) = 0.5 * arma::sqrt(eventRate.t());

    // Rcpp::Rcout << "calcWeights finished\n";

    // weightMat = weights * weights.t();

    // arma::mat wZDat = coxWeights % zDat;
    // arma::rowvec wZSd = arma::stddev(wZDat, 0, 0);

    // for (arma::uword i = 0; i < p; i++) {
    // 	double Ex2 = arma::mean(arma::square(xDat.col(i)));
    // 	for (arma::uword m = 0; m < r; m++) {
    // 	    weightMat(i, p+q+m) += 0.5 * wZSd(m) * std::sqrt(Ex2);
    // 	    weightMat(p+q+m, i) = weightMat(i, p+q+m); 
    // 	}
    // }

    // for (arma::uword j = 0; j < q; j++) {
    //     for (arma::uword m = 0; m < r; m++) {
    // 	    weightMat(p+j, p+q+m) += 0.5 * wZSd(m);
    // 	    weightMat(p+q+m, p+j) = weightMat(p+j, p+q+m); 
    // 	}
    // }

    // Rcpp::Rcout << "Censored Variable Weights:\n" << weightMat.submat(0, p+q, p+q-1, p+q+r-1) << std::endl;
    
    // arma::mat wZDat = coxWeights % zDat;
    // arma::rowvec weightSums = arma::sum(coxWeights,0);
    // arma::rowvec zDatBar = arma::sum(wZDat) / weightSums;
    // arma::rowvec zDat2Bar = arma::sum(wZDat % zDat) / weightSums;
    // arma::rowvec zDatWSd = arma::sqrt(weightSums / (weightSums-1) * (zDat2Bar - arma::square(zDatBar)));
    
    // weights.subvec(p+q, p+q+r-1) = weightSums / ((double) n) * zDatWSd;
      // (2 * arma::sqrt(arma::mean(arma::conv_to<arma::mat>::from(censMat), 0)));

    // arma::vec w2;
    // double wsum, wZmu, wZvar, normConstant, eventRate;
    // for (arma::uword i = 0; i < r; i++) {

    // 	// Rcpp::Rcout << "Z(" << i << ")\n"
    // 	// 	    << "weight" << coxWeights.col(i).t()
    // 	// 	    << "Z" << zDat.col(i).t() << std::endl;
    // 	eventRate = 1; // arma::mean(arma::conv_to<arma::vec>::from(censMat.col(i)));
    // 	w2 = arma::square(coxWeights.col(i));
    // 	wsum = arma::sum(w2);
    // 	wZmu = arma::sum(coxWeights.col(i) % zDat.col(i)) / arma::sum(coxWeights.col(i));
    // 	// wZmu = arma::mean(zDat.col(i));
    // 	normConstant = (wsum / (wsum*wsum - arma::sum(arma::square(w2))));
    // 	RcppThread::Rcout << normConstant << std::endl << (1/(double)(n-1)) / normConstant << std::endl;
    // 	wZvar = arma::sum(w2 % arma::square(zDat.col(i) - wZmu)) * normConstant;
    // 	weights(p+q+i) = std::sqrt(wZvar/eventRate);
    // }
    
    // weightMat = weights * weights.t();

    // RcppThread::Rcout << std::endl << "Z Weights:\n" << weights.subvec(p+q, p+q+r-1).t() << std::endl;

    // RcppThread::Rcout << std::endl << "Z Weights 2:\n" << arma::stddev(coxWeights % zDat, 0, 0) << std::endl;

    // RcppThread::Rcout << std::endl << "Z Weights 3:\n" << arma::stddev(arma::square(coxWeights) % zDat, 0, 0) << std::endl;

    // Rcpp::Rcout << std::endl << "GradWeights Psi:\n" << weightMat.submat(p, p+q, p+q-1, p+q+r-1) << std::endl;

    // Cox variable weights will be recalculated with each new Z and coxWeight
    // calcZWeights();
}

// void CoxMGM::calcZWeights() {

//     arma::mat eta = xDat * par.gamma + dDat * par.psi;
//     eta.each_col( [](arma::vec& c) { c -= arma::mean(c); } );
    
//     arma::mat coxgrad(arma::size(eta), arma::fill::zeros);
//     arma::mat diagHess(arma::size(eta), arma::fill::zeros);

//     arma::vec oldCoxloss = coxGradHess(eta, coxgrad, diagHess);

//     coxWeights = -diagHess;
//     diagHess.replace(0, -1e-10);
//     zDat = eta - coxgrad / diagHess;

//     weights.subvec(p+q, p+q+r-1) = arma::stddev(coxWeights % zDat, 0, 0);


//     for (int m = 0; m < r; m++) {
// 	censWeights = arma::vec(p+q+r, arma::fill::ones);

// 	//Continuous variable weights are standard deviations
// 	for (arma::uword i = 0; i < p; i++) {
// 	    censWeights(i) = arma::stddev(xDat.col(i).elem(arma::find(censMat.col(m) == 1)));
// 	}

// 	//Discrete variable weights for each variable-category pair are p(1-p) where p is the percentage of times that category appears
// 	for (arma::uword j = 0; j < q; j++) {
// 	    double curWeight = 0;
// 	    for (int k = 0; k < l[j]; k++) {
// 		arma::vec equalityVec = arma::vec(yDat.col(j)).transform( [k](double val) { return val == k+1 ? 1 : 0; } );
// 		double curp = arma::mean(equalityVec.elem(arma::find(censMat.col(m) == 1)));
// 		curWeight += curp * (1-curp);
// 	    }
// 	    censWeights(p+j) = std::sqrt(curWeight);
// 	}

	

//     }

//     // arma::vec w2;
//     // double wsum, wZmu, wZvar, normConstant, eventRate;
//     // for (arma::uword i = 0; i < r; i++) {

//     // 	// Rcpp::Rcout << "Z(" << i << ")\n"
//     // 	// 	    << "weight" << coxWeights.col(i).t()
//     // 	// 	    << "Z" << zDat.col(i).t() << std::endl;

//     // 	eventRate = 0.5; // arma::mean(arma::conv_to<arma::vec>::from(censMat.col(i)));
//     // 	w2 = arma::square(coxWeights.col(i));
//     // 	wsum = arma::sum(w2);
//     // 	wZmu = arma::sum(coxWeights.col(i) % zDat.col(i)) / arma::sum(coxWeights.col(i));
//     // 	normConstant = (wsum / (wsum*wsum - arma::sum(arma::square(w2))));
//     // 	wZvar = arma::sum(w2 % arma::square(zDat.col(i) - wZmu)) * normConstant;
//     // 	weights(p+q+i) = std::sqrt(wZvar/eventRate);
//     // }
    
//     // weightMat = weights * weights.t();

//     // Rcpp::Rcout << std::endl << "Z Weights:\n" << weights.subvec(p+q, p+q+r-1).t() << std::endl;
    
//     // arma::vec wz;
//     // double wXmu, wXvar, wZmu, wZvar, wsum, normConstant;
//     // for (arma::uword i = 0; i < r; i++) {
//     // 	// wz = coxWeights.col(i) % zDat.col(i);
//     // 	wsum = arma::sum(coxWeights.col(i));
//     // 	if (wsum==0) continue;
//     // 	normConstant = (wsum / (wsum*wsum - arma::sum(arma::square(coxWeights.col(i)))));
//     // 	// normConstant = 1 / wsum;
//     // 	wZmu = arma::sum(coxWeights.col(i) % zDat.col(i)) / wsum;
//     // 	// arma::vec w2 = arma::square(coxWeights.col(i));
//     // 	// Zmu = arma::mean(zDat.col(i));
//     // 	wZvar = arma::sum(coxWeights.col(i) % arma::square((zDat.col(i) - wZmu))) * normConstant;
//     // 	// wZvar *= normConstant;
//     // 	// wvar /= (n-1);
//     //     // weightMat.col(i+p+q) = std::max(std::sqrt(wZvar), 1e-5);
//     // 	// weightMat.row(i+p+q) = std::max(std::sqrt(wZvar), 1e-5);

//     // 	Rcpp::Rcout << "curWeight Z(" << i << "):   " << std::sqrt(wZvar) << std::endl;

//     // 	for (arma::uword j = 0; j < p; j++) {
//     // 	    wXmu = arma::sum(coxWeights.col(i) % xDat.col(j)) / wsum;
//     // 	    // Xmu = arma::mean(xDat.col(i));
//     // 	    wXvar = arma::sum(coxWeights.col(i) % arma::square((xDat.col(j)-wXmu))) * normConstant;
//     // 	    // wXvar *= normConstant;
//     // 	    Rcpp::Rcout << "  curWeight X(" << j << "):   " << std::sqrt(wXvar) << std::endl;
	    
//     // 	    weightMat(j,i+p+q) = std::sqrt(wZvar) * std::sqrt(wXvar);
//     // 	    // weightMat(j,i+p+q) = std::sqrt(arma::sum(arma::square(xDat.col(j) % (wZmu - coxWeights.col(i) % zDat.col(i)) + zDat.col(i) % (wXmu - coxWeights.col(i) % xDat.col(j)))) / n);
//     // 	    weightMat(i+p+q,j) = weightMat(j,i+p+q);
//     // 	}

//     // 	//Discrete variable weights for each variable-category pair are p(1-p) where p is the percentage of times that category appears
//     // 	for (arma::uword j = 0; j < q; j++) {
//     // 	    double curWeight = 0;
//     // 	    double varSum = 0;
//     // 	    for (int k = 0; k < l[j]; k++) {
//     // 		arma::vec equalityVec = arma::vec(yDat.col(j)).transform( [k](double val) { return val == k+1 ? 1 : 0; } );
//     // 		double curp = arma::sum(coxWeights.col(i) % equalityVec) / wsum;
//     // 		// double curp = arma::mean(equalityVec);
//     // 		// curWeight += arma::sum(arma::square(equalityVec % (wZmu - coxWeights.col(i) % zDat.col(i)) + zDat.col(i) % (curp - coxWeights.col(i) % equalityVec))) / n;
//     // 		curWeight += arma::sum(coxWeights.col(i) % arma::square((equalityVec - curp))) * normConstant;
//     // 		// varSum += curp * (1-curp);
//     // 	    }
//     // 	    Rcpp::Rcout << "  curWeight Y(" << j << "):   " << std::sqrt(curWeight) << std::endl;
//     // 	    // Rcpp::Rcout << "  wp(1-wp) Y(" << j << "):   " << varSum << std::endl;
//     // 	    // curWeight *= normConstant;
//     // 	    weightMat(p+j, p+q+i) = std::sqrt(wZvar) * std::sqrt(curWeight);
//     // 	    // weightMat(p+j, p+q+i) = std::sqrt(curWeight);
//     // 	    weightMat(p+q+i, p+j) = weightMat(p+j, p+q+i);
//     // 	}
//     // }
    
//     // Rcpp::Rcout << std::endl << "GradWeights Gamma:\n" << weightMat.submat(0, p+q, p-1, p+q+r-1) << std::endl;

//     // Rcpp::Rcout << std::endl << "GradWeights Psi:\n" << weightMat.submat(p, p+q, p+q-1, p+q+r-1) << std::endl;

//     // for (arma::uword i = 0; i < r; i++) {
//     // 	wz = coxWeights.col(i) % zDat.col(i);
//     // 	// wmean = arma::mean(coxWeights.col(i));
//     // 	// if (wmean==0) continue;
//     // 	weights(i+p+q) = std::max(std::sqrt(2) * arma::stddev(wz), 0.01);
//     // }
//     // Rcpp::Rcout << "GradWeights: " << weights.t();
// }

// convert discrete data (in yDat) to a matrix of dummy variables (stored in dDat)
void CoxMGM::makeDummy() {
    dDat = arma::mat(n, lsum, arma::fill::zeros);
    for(int i = 0; i < q; i++) {
        for(int j = 0; j < l[i]; j++) {
            arma::vec curCol = arma::vec(yDat.col(i)).transform( [j](double val) { return val == j+1 ? 1 : 0; } );
            if (arma::sum(curCol) == 0)
                throw std::invalid_argument("Discrete data is missing a level: variable " + variables[p+i].getName() + " level " + std::to_string(j));
            dDat.col(lcumsum[i]+j) = curCol;
        }
    }
}

// checks if yDat is zero indexed and converts to 1 index. zscores x
void CoxMGM::fixData() {
    double ymin = yDat.min();
    if(ymin < 0 || ymin > 1)
        throw std::invalid_argument("Discrete data must be either zero or one indexed. Found min index: " + std::to_string(ymin));
    
    if (ymin == 0) {
        yDat += 1;
    }

    // z-score of columns of x
    xDat.each_col( [](arma::vec& c) {c = (c - arma::mean(c)) / arma::stddev(c); } );
}

arma::vec CoxMGM::coxGradHess(arma::mat& eta, arma::mat& grad, arma::mat& diagHess) {
    arma::mat theta = arma::exp(eta);
    arma::mat theta_weight = arma::zeros(n, r);
    arma::mat theta_weight2 = arma::zeros(n, r);
    arma::vec loss = arma::zeros(r);
    arma::uvec censor, order, index;

    // Rcpp::Rcout << "coxGradHess called\n";
    
    // arma::uvec order = target->getOrder();
    // arma::uvec H = target->getH();
    // arma::uvec censor = target->getCensor();
    double HsumTheta, nEvents, sub, d, phi, theta_weight_sum, theta_weight2_sum, dSum, dSum2, rs_sum;
    
    // arma::vec rs_sum = arma::sum(theta,0).t();

    grad += arma::conv_to<arma::mat>::from(censMat);

    for (int m = 0; m < r; m++) {
	arma::uvec col = { (uint) m };
	for (int strat = 0; strat < numStrata(m); strat++) {
	    // censor = censList[m][strat];
	    // order = orderList[m][strat];
	    // index = idxList[m][strat];
	    rs_sum = arma::accu(theta(idxList[m][strat], col));
	    theta_weight_sum = 0;
	    theta_weight2_sum = 0;
	    int i = 0;
	    for (int j = 0; j < HList[m][strat].n_elem; j++) {
		HsumTheta = 0;
		nEvents = 0;
		sub = 0;
		dSum = 0;
		dSum2 = 0;

		for (int k = 0; k < HList[m][strat][j]; k++) {
		    if (censList[m][strat](orderList[m][strat](i+k))) {
			nEvents++;
			HsumTheta += theta(idxList[m][strat](orderList[m][strat](i+k)), m);
			loss[m] += eta(idxList[m][strat](orderList[m][strat](i+k)), m);
		    }
		    sub += theta(idxList[m][strat](orderList[m][strat](i+k)), m);
		}

		if (nEvents > 0) {
		    if (sub - rs_sum > 1e-5) {
			if ((HList[m][strat][j] + i) != n) {
			    arma::uvec indices = orderList[m][strat](arma::span(i,n-1));
			    // arma::uvec col = { (uint) m };
			    rs_sum = arma::accu(theta(indices, col));
			    if (sub - rs_sum > 1e-5) {
				throw std::runtime_error("Error in CoxMGM coxGradHess, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
			    }
			} else {
			    sub = rs_sum;
			}
		    }

		    for (int l = 0; l < nEvents; l++) {
			d = l / ((double) nEvents);
			phi = rs_sum - d * HsumTheta;
			loss[m] -= std::log(phi);
			theta_weight_sum += 1 / phi;
			dSum += d / phi;
			theta_weight2_sum += 1 / std::pow(phi, 2);
			dSum2 += (2*d - std::pow(d, 2)) / std::pow(phi, 2);
		    }
		}

		for (int k = 0; k < HList[m][strat][j]; k++) {
		    theta_weight(idxList[m][strat](orderList[m][strat](i+k)), m) = theta_weight_sum - censList[m][strat](orderList[m][strat](i+k)) * dSum;
		    theta_weight2(idxList[m][strat](orderList[m][strat](i+k)), m) = theta_weight2_sum - censList[m][strat](orderList[m][strat](i+k)) * dSum2;
		}
	
		i += HList[m][strat][j];
		rs_sum -= sub;
	    }
	}
    }

    grad -= theta % theta_weight;

    diagHess = arma::square(theta) % theta_weight2 - theta % theta_weight;

    // diagHess.replace(0, -1e-10);

    // Rcpp::Rcout << "coxGradHess finished.\n";

    return loss;
}


double CoxMGM::calcLambdaMax() {
    // Rcpp::Rcout << "Running calcLambdaMax...\n";
    lcumsum = std::vector<int>(l.size()+1);
    lcumsum[0] = 0;
    for(int i = 0; i < l.size(); i++){
        lcumsum[i+1] = lcumsum[i] + l[i];
    }
    lsum = lcumsum[l.size()];

    arma::mat beta((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros);

    CoxMGMParams nullParams(
        arma::mat((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros),  // beta
        arma::vec((int) xDat.n_cols,                    arma::fill::ones),   // betad
        arma::mat(lsum,              (int) xDat.n_cols, arma::fill::zeros),  // theta
        arma::mat(lsum,              lsum,              arma::fill::zeros),  // phi
	arma::mat((int) xDat.n_cols, (int) zDat.n_cols, arma::fill::zeros),  // gamma
	arma::mat(lsum,              (int) zDat.n_cols, arma::fill::zeros),  // psi
        arma::vec((int) xDat.n_cols,                    arma::fill::zeros),  // alpha1
        arma::vec(lsum,                                 arma::fill::zeros),  // alpha2
	arma::vec((int) zDat.n_cols,                    arma::fill::zeros)   // alpha3
    );

    arma::vec nullParams1D = nullParams.toMatrix1D();
    arma::vec nullGrad1D = smoothGradient(nullParams1D);

    CoxMGMParams nullGrad(nullGrad1D, p, lsum, r);

    // double lambdaMax = 0.0;

    arma::mat weightMat = weights * weights.t();

    //weight beta
    //betaw = (wv(1:p)'*wv(1:p)).*abs(beta);
    //betanorms=sum(betaw(:));
    arma::mat betaNorms = arma::abs(nullGrad.beta) / arma::mat(weightMat.submat(0, 0, p-1, p-1));
    double lambdaMax = betaNorms.max();

    // Rcpp::Rcout << "lambdaMax after beta = " << lambdaMax << std::endl;

    // Rcpp::Rcout << "NSV betaNorms = " << betaNorms << std::endl;

    /*
    thetanorms=0;
    for s=1:p
        for j=1:q
            tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
            thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
        end
    end
    */
    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            arma::vec tempVec = nullGrad.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
	    // Rcpp::Rcout << "mean theta(" << i << ", " << j << ") = " << arma::mean(tempVec) << std::endl;
            lambdaMax = std::max(lambdaMax, arma::norm(tempVec, 2) / weightMat(i, p+j));
        }
    }
    // Rcpp::Rcout << "NSV thetaNorms = " << thetaNorms << std::endl;

    // Rcpp::Rcout << "lambdaMax after theta = " << lambdaMax << std::endl;

    /*
    for r=1:q
        for j=1:q
            if r<j
                tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
                tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
                phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
                phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
            end
        end
    end
    */
    for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
        for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
            arma::mat tempMat = nullGrad.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
	    // Rcpp::Rcout << "mean phi(" << i << ", " << j << ") = " << arma::mean(arma::mean(tempMat)) << std::endl;
            lambdaMax = std::max(lambdaMax, arma::norm(tempMat, "fro") / weightMat(p+i, p+j));
        }
    }

    // Rcpp::Rcout << "Independent Gamma^2 for Survival0:\n" << arma::square(nullGrad.gamma.col(0).t());

    // Rcpp::Rcout << "lambdaMax after phi = " << lambdaMax << std::endl;

    double gammaMax = (abs(nullGrad.gamma) / weightMat.submat(0, p+q, p-1, p+q+r-1)).max();

    lambdaMax = std::max(lambdaMax, gammaMax);

    // Rcpp::Rcout << "lambdaMax after gamma = " << lambdaMax << std::endl;


    // Rcpp::Rcout << "Independent Psi^2 for Survival0:\n";
    
    for (arma::uword i = 0; i < r; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            arma::subview_col<double> tempVec = nullGrad.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
	    // if (i == 0) Rcpp::Rcout << arma::sum(tempVec % tempVec) << "    ";
	    lambdaMax = std::max(lambdaMax, arma::norm(tempVec, 2) / weightMat(p+j, p+q+i));
        }
    }
    // Rcpp::Rcout << "\n";

    // Rcpp::Rcout << "lambdaMax after psi = " << lambdaMax << std::endl;

    // Rcpp::Rcout << "calcLambdaMax Finished\n";

    return lambdaMax;

}


/**
 * Calculate value of g(X) and gradient of g(X) at the same time for efficiency reasons.
 *
 * @param X input Vector
 * @param Xout gradient of g(X)
 * @return value of g(X)
 */
double CoxMGM::smooth(arma::vec& parIn, arma::vec& gradOutVec) {
    CoxMGMParams par(parIn, p, lsum, r);

    // Rcpp::Rcout << "par: \n" << par << std::endl;

    // for (arma::uword j = 0; j < q; j++) {
    //  	par.psi.row(lcumsum[j]).fill(0);
    // }

    // Rcpp::Rcout << "Running smooth...\n";

    arma::mat eta = xDat * par.gamma + dDat * par.psi;
    eta.each_row() += par.alpha3.t();

    // arma::mat curEta(eta);
    // curEta.each_col( [](arma::vec& c) { c -= arma::mean(c); } );

    
    // arma::mat curCoxgrad(arma::size(eta), arma::fill::zeros);
    // arma::mat curDiagHess(arma::size(eta), arma::fill::zeros);

    // arma::vec curCoxloss = coxGradHess(curEta, curCoxgrad, curDiagHess);

    // arma::mat curW(-curDiagHess);
    // curDiagHess.replace(0, -1e-10);
    // arma::mat curWZ = curW % (curEta - curCoxgrad / curDiagHess;

    // Rcpp::Rcout << "eta: " << eta << std::endl;
    // Rcpp::Rcout << "Z: " << zDat << std::endl;
    // Rcpp::Rcout << "coxWeights: " << coxWeights << std::endl;

    CoxMGMParams gradOut;

    // gradOut.gamma = arma::zeros(p, r);
    // gradOut.psi = arma::zeros(lsum, r);
    // gradOut.xid = arma::zeros(r);
    // gradOut.alpha3 = arma::zeros(r);

    for(arma::uword i = 0; i < par.betad.size(); i++) {
        if(par.betad(i) <= 0)
            return arma::datum::inf;
    }

    //beta=beta-diag(diag(beta));
    //for r=1:q
    //  phi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    //beta=triu(beta); phi=triu(phi);
    //beta=beta+beta';
    //phi=phi+phi';
    par.beta = arma::symmatu(par.beta);
    par.beta.diag(0).zeros();

    for (arma::uword i = 0; i < q; i++) {
        par.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    par.phi = arma::symmatu(par.phi);

    //Xbeta=X*beta*diag(1./betad);
    arma::mat divBetaD = arma::diagmat(arma::vec(p, arma::fill::ones) / par.betad);
    arma::mat xBeta = xDat * par.beta * divBetaD;

    //Dtheta=D*theta*diag(1./betad);
    arma::mat dTheta = dDat * par.theta * divBetaD;

    arma::mat zGamma = wzDat * par.gamma.t() * divBetaD;

    // Rcpp::Rcout << "zGamma: \n" << zGamma << std::endl;

    // zGamma.fill(0);

    // for (arma::uword m = 0; m < r; m++) {
    // 	zGamma += arma::diagmat(coxWeights.col(m)) * zDat.col(m) * par.gamma.col(m).t() * divBetaD;
    // }
    // // zGamma = zGamma * divBetaD;
    
    // Rcpp::Rcout << "zGamma2: \n" << zGamma << std::endl;

    // Squared loss
    //tempLoss =  (X-e*alpha1'-Xbeta-Dtheta) = -res (in gradient code)
    arma::mat tempLoss(n, xDat.n_cols);

    //wxprod=X*(theta')+D*phi+e*alpha2';
    arma::mat wxProd = xDat * par.theta.t() + dDat * par.phi + wzDat * par.psi.t();

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < xDat.n_cols; j++) {
            tempLoss(i, j) = xDat(i,j) - par.alpha1(j) - xBeta(i,j) - dTheta(i,j) - zGamma(i,j);
        }
        for (arma::uword j = 0; j < dDat.n_cols; j++) {
            wxProd(i,j) = wxProd(i,j) + par.alpha2(j);
        }
    }

    //sqloss=-n/2*sum(log(betad))+...
    //.5*norm((X-e*alpha1'-Xbeta-Dtheta)*diag(sqrt(betad)),'fro')^2;
    double sqloss = -n/2.0*arma::sum(arma::log(arma::vec(par.betad))) +
                    0.5 * std::pow(arma::norm(tempLoss * arma::diagmat(arma::sqrt(arma::vec(par.betad))), "fro"), 2);

    //ok now tempLoss = res
    tempLoss *= -1;

    //gradbeta=X'*(res);
    gradOut.beta = xDat.t() * tempLoss;

    //gradbeta=gradbeta-diag(diag(gradbeta)); % zero out diag
    //gradbeta=tril(gradbeta)'+triu(gradbeta);
    gradOut.beta = arma::trimatl(gradOut.beta, -1).t() + arma::trimatu(gradOut.beta, 1);

    //gradalpha1=diag(betad)*sum(res,0)';
    gradOut.alpha1 = arma::diagmat(par.betad) * arma::sum(tempLoss, 0).t();

    //gradtheta=D'*(res);
    gradOut.theta = dDat.t() * tempLoss;

    // gradOut.gamma = (wzDat.t() * tempLoss).t();
    // gradOut.gamma = (5/2.0) * ((coxWeights % zDat).t() * (tempLoss / (4 + arma::square(resid)))).t();
    
    // gradOut.gamma.each_row( [this] (arma::rowvec& r) { r %= (fitWeight / (1 + fitWeight)); } );
    // gradOut.gamma.each_row( [this] (arma::rowvec& r) { r %= fitWeight; } );

    // gradOut.gamma = arma::mat(p, r, arma::fill::zeros);

    // Rcpp::Rcout << "Gamma grad 1: \n" << (wzDat.t() * tempLoss).t();

    // categorical loss
    /*catloss=0;
    wxprod=X*(theta')+D*phi+e*alpha2'; %this is n by Ltot
    for r=1:q
        wxtemp=wxprod(:,Lsum(r)+1:Lsum(r)+L(r));
        denom= logsumexp(wxtemp,2); %this is n by 1
        catloss=catloss-sum(wxtemp(sub2ind([n L(r)],(1:n)',Y(:,r))));
        catloss=catloss+sum(denom);
    end
    */
    double catloss = 0;
    for (arma::uword i = 0; i < yDat.n_cols; i++) {
        arma::subview<double> wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));

        //need to copy init values for calculating nll
        arma::mat wxTemp0(wxTemp);

        // does this need to be done in log space??
        wxTemp = arma::exp(wxTemp.each_col() - logsumexp(wxTemp0));
        // arma::vec invDenom = arma::vec(n, arma::fill::ones) / arma::sum(wxTemp, 1);
        // wxTemp = arma::diagmat(invDenom) * wxTemp;

        for (arma::uword k = 0; k < n; k++) {
            arma::vec curRow0(arma::conv_to<arma::vec>::from(wxTemp0.row(k)));

            catloss -= curRow0((arma::uword) yDat(k,i) - 1);
            catloss += logsumexp(curRow0);

            //wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))=wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))-1;
            wxTemp(k, (arma::uword) yDat(k,i)-1) -= 1;
        }
    }

    //gradalpha2=sum(wxprod,0)';
    gradOut.alpha2 = arma::sum(wxProd, 0).t();

    //gradw=X'*wxprod;
    arma::mat gradW = xDat.t() * wxProd;

    //gradtheta=gradtheta+gradw';
    gradOut.theta += gradW.t();

    //gradphi=D'*wxprod;
    gradOut.phi = dDat.t() * wxProd;

    // gradOut.psi = (wzDat.t() * wxProd).t();
    // gradOut.psi = (5/2.0) * ((coxWeights % zDat).t() * (wxProd / (4 + arma::square(catResid)))).t();
    // gradOut.psi = arma::mat(lsum, r, arma::fill::zeros);

    // gradOut.psi.each_row( [this] (arma::rowvec& r) { r /= wzSd; } );
    // gradOut.psi.each_row( [this] (arma::rowvec& r) { r %= (fitWeight / (1 + fitWeight)); } );
    // gradOut.psi.each_row( [this] (arma::rowvec& r) { r %= fitWeight; } );

    // Rcpp::Rcout << "Psi grad 1: \n" << (wzDat.t() * wxProd).t();

    //zero out gradphi diagonal
    //for r=1:q
    //gradphi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    for (arma::uword i = 0; i < q; i++) {
        gradOut.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    //gradphi=tril(gradphi)'+triu(gradphi);
    gradOut.phi = arma::trimatl(gradOut.phi, 0).t() + arma::trimatu(gradOut.phi, 0);

    // Cox loss
    arma::mat coxresid = eta - zDat;

    gradOut.gamma = xDat.t() * (coxWeights % coxresid);
    gradOut.psi = dDat.t() * (coxWeights % coxresid);

    // Rcpp::Rcout << "Gamma grad 2: \n" << xDat.t() * (coxWeights % coxresid);
    // Rcpp::Rcout << "Gamma grad p1/p2: \n" << gradOut.gamma / (xDat.t() * (coxWeights % coxresid));

    // Rcpp::Rcout << "Average gamma grad p1/p2:" << arma::exp(arma::mean(arma::log(arma::abs(gradOut.gamma)) - arma::log(arma::abs(xDat.t() * (coxWeights % coxresid))), 0)) << std::endl;

    // Rcpp::Rcout << "Psi grad 2: \n" << dDat.t() * (coxWeights % coxresid);
    // Rcpp::Rcout << "Psi grad p1/p2: \n" << gradOut.psi / (dDat.t() * (coxWeights % coxresid));

    // Rcpp::Rcout << "Average psi grad p1/p2:" << arma::exp(arma::mean(arma::log(arma::abs(gradOut.psi)) - arma::log(arma::abs(dDat.t() * (coxWeights % coxresid))), 0)) << std::endl;
    
    double coxloss = 0;
    for (arma::uword i = 0; i < zDat.n_cols; i++) {
	// coxresid.col(i) += par.alpha3(i);
	coxloss -= arma::as_scalar(-0.5 * coxresid.col(i).t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i));
	coxloss += 0.5 * arma::sum(arma::square(coxgrad.col(i)) / diagHess.col(i));
	coxloss -= oldCoxloss(i);
        // gradOut.gamma.col(i) += (1 / (1 + fitWeight(i))) * xDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);
	// gradOut.psi.col(i) += (1 / (1 + fitWeight(i))) * dDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);
	// gradOut.gamma.col(i) += xDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);

	// arma::mat psiGradTemp = arma::diagmat(coxWeights.col(i) % coxresid.col(i)) * dDat;

	// for (arma::uword j = 0; j < q; j++) {
	//     arma::subview<double> wxTemp = psiGradTemp.submat(0, lcumsum[j], n-1, lcumsum[j+1]-1);
	//     wxTemp.each_row( [this] (arma::rowvec& r) { r -= arma::mean(r); } );
	// }
	
        // gradOut.psi.col(i) += arma::mean(psiGradTemp,0);
	
	// gradOut.psi.col(i) += dDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);


	// gradOut.gamma.col(i) /= arma::sum(coxWeights.col(i));
	// gradOut.psi.col(i) /= arma::sum(coxWeights.col(i));
    }


    // for (arma::uword j = 0; j < q; j++) {
    // 	gradOut.psi.row(lcumsum[j]).fill(0);
    // }

    gradOut.alpha3 = arma::sum(coxWeights % coxresid, 0).t();
    
    /*
    for s=1:p
        gradbetad(s)=-n/(2*betad(s))+1/2*norm(res(:,s))^2-res(:,s)'*(Xbeta(:,s)+Dtheta(:,s));
    end
    */
    gradOut.betad = arma::vec(xDat.n_cols);
    for(arma::uword i = 0; i < p; i++){
        gradOut.betad(i) = -n / (2.0 * par.betad(i)) + std::pow(arma::norm(tempLoss.col(i), 2), 2) / 2.0 -
	    arma::as_scalar(tempLoss.col(i).t() * (xBeta.col(i) + dTheta.col(i) + zGamma.col(i) + par.alpha1(i) / par.betad(i)));
    }

    // Rcpp::Rcout << "beta grad:\n" << gradOut.beta;

    // Rcpp::Rcout << "Sum of coxWeights:\n" << arma::sum(coxWeights, 0);

    gradOut.alpha1 /= (double) n;
    gradOut.alpha2 /= (double) n;
    gradOut.alpha3 /= (double) n; // arma::sum(coxWeights,0).t();
    gradOut.betad /= (double) n;
    gradOut.beta /= (double) n;
    gradOut.theta /= (double) n;
    gradOut.phi /= (double) n;
    gradOut.gamma /= (double) n;
    gradOut.psi /= (double) n;

    // Rcpp::Rcout << "gradOut: \n" << gradOut << std::endl;

    if (pDummy != 0) {
	gradOut.alpha1.fill(0);
	gradOut.betad.fill(0);
	gradOut.beta.fill(0);
	gradOut.theta.fill(0);
	gradOut.gamma.fill(0);
    }

    if (qDummy != 0) {
	gradOut.alpha2.fill(0);
	gradOut.theta.fill(0);
	gradOut.phi.fill(0);
	gradOut.psi.fill(0);
    }

    gradOutVec = gradOut.toMatrix1D();

    // Rcpp::Rcout << "Finished\n";
    
    return (sqloss + catloss + coxloss)/((double) n);
}

double CoxMGM::smoothValue(arma::vec& parIn) {
    CoxMGMParams par(parIn, p, lsum, r);

    // Rcpp::Rcout << "smoothValue called\n";

    // for (arma::uword j = 0; j < q; j++) {
    //  	par.psi.row(lcumsum[j]).fill(0);
    // }

    // Rcpp::Rcout << "par: \n" << par << std::endl;

    arma::mat eta = xDat * par.gamma + dDat * par.psi;
    eta.each_row() += par.alpha3.t();
    // eta.each_col( [](arma::vec& c) { c -= arma::mean(c); } );
    
    // arma::mat coxgrad(arma::size(eta), arma::fill::zeros);
    // arma::mat diagHess(arma::size(eta), arma::fill::zeros);

    // arma::vec oldCoxloss = coxGradHess(eta, coxgrad, diagHess);

    // coxWeights = -diagHess;
    // diagHess.replace(0, -1e-10);
    // zDat = eta - coxgrad / diagHess;
    
    // Rcpp::Rcout << "eta: " << eta << std::endl;
    // Rcpp::Rcout << "Z: " << zDat << std::endl;
    // Rcpp::Rcout << "coxWeights: " << coxWeights << std::endl;

    for(arma::uword i = 0; i < par.betad.size(); i++) {
        if(par.betad(i) <= 0)
            return std::numeric_limits<double>::infinity();
    }

    // for(arma::uword i = 0; i < par.xid.size(); i++) {
    //     if(par.xid(i) <= 0)
    //         return std::numeric_limits<double>::infinity();
    // }

    //double nll = 0;
    //int n = xDat.rows();
    //beta=beta+beta';
    //phi=phi+phi';
    par.beta = arma::symmatu(par.beta);
    par.beta.diag(0).zeros();

    for (arma::uword i = 0; i < q; i++) {
        par.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    par.phi = arma::symmatu(par.phi);

    arma::mat divBetaD = arma::diagmat(arma::vec(p, arma::fill::ones) / par.betad);

    //Xbeta=X*beta*diag(1./betad);
    //Dtheta=D*theta*diag(1./betad);
    arma::mat xBeta = xDat * par.beta * divBetaD;
    arma::mat dTheta = dDat * par.theta * divBetaD;
    arma::mat zGamma = wzDat * par.gamma.t() * divBetaD;

    // Rcpp::Rcout << "xBeta: \n" << xBeta << std::endl;
    // Rcpp::Rcout << "dTheta: \n" << dTheta << std::endl;
    // Rcpp::Rcout << "zGamma: \n" << dTheta << std::endl;

    // Squared loss
    //sqloss=-n/2*sum(log(betad))+...
    //.5*norm((X-e*alpha1'-Xbeta-Dtheta)*diag(sqrt(betad)),'fro')^2;
    arma::mat tempLoss(n, xDat.n_cols);

    //wxprod=X*(theta')+D*phi+e*alpha2';
    arma::mat wxProd = xDat * par.theta.t() + dDat * par.phi + wzDat * par.psi.t();

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < xDat.n_cols; j++) {
            tempLoss(i, j) = xDat(i,j) - par.alpha1(j) - xBeta(i,j) - dTheta(i,j) - zGamma(i,j);
        }
        for (arma::uword j = 0; j < dDat.n_cols; j++) {
            wxProd(i,j) = wxProd(i,j) + par.alpha2(j);
        }
    }

    // Rcpp::Rcout << "tempLoss: \n" << tempLoss << std::endl;

    double sqloss = -n/2.0*arma::sum(arma::log(arma::vec(par.betad))) +
                    0.5 * std::pow(arma::norm(tempLoss * arma::diagmat(arma::sqrt(arma::vec(par.betad))), "fro"), 2);
    
    // Rcpp::Rcout << "sqloss: " << sqloss << std::endl;

    // categorical loss
    /*catloss=0;
    wxprod=X*(theta')+D*phi+e*alpha2'; %this is n by Ltot
    for r=1:q
        wxtemp=wxprod(:,Lsum(r)+1:Lsum(r)+L(r));
        denom= logsumexp(wxtemp,2); %this is n by 1
        catloss=catloss-sum(wxtemp(sub2ind([n L(r)],(1:n)',Y(:,r))));
        catloss=catloss+sum(denom);
    end
    */
    double catloss = 0;
    for (arma::uword i = 0; i < yDat.n_cols; i++) {
        arma::subview<double> wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));
        for (arma::uword k = 0; k < n; k++) {
            arma::vec curRow = arma::conv_to<arma::vec>::from(wxTemp.row(k));

            catloss -= curRow((arma::uword) yDat(k,i) - 1);
            catloss += logsumexp(curRow);
        }
    }

    // Rcpp::Rcout << "catloss: " << catloss << std::endl;

    // Cox loss
    arma::mat coxresid = eta - zDat;
    double coxloss = 0;
    for (arma::uword i = 0; i < zDat.n_cols; i++) {
	// coxresid.col(i) += par.alpha3(i);
	coxloss -= arma::as_scalar(-0.5 * coxresid.col(i).t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i));
	coxloss += 0.5 * arma::sum(arma::square(coxgrad.col(i)) / diagHess.col(i));
	coxloss -= oldCoxloss(i);
    }

    // Rcpp::Rcout << "Old Cox Loss: " << arma::sum(oldCoxloss) << std::endl;
    // Rcpp::Rcout << "coxloss: " << coxloss << std::endl;

    return (sqloss + catloss + coxloss) / (double) n;
}

/**
 * Calculate value of h(X) and proxOperator of h(X) at the same time for efficiency reasons.
 *
 * @param t positive parameter for prox operator
 * @param X input vector
 * @param Xout vector solution to prox_t(X)
 * @return value of h(X)
 */
double CoxMGM::nonSmooth(double t, arma::vec& X, arma::vec& pX) {
    double nonSmooth = 0;

    // Rcpp::Rcout << "nonSmooth called\n";

    // Rcpp::Rcout << "Running nonSmooth...\n";

    arma::vec tlam = lambda * t;

    //Constructor copies and checks dimension
    //par is a copy so we can update it
    CoxMGMParams par(X, p, lsum, r);

    // calcZWeights();

    // Rcpp::Rcout << "NS par: \n" << par << std::endl;

    //penbeta = t(1).*(wv(1:p)'*wv(1:p));
    //betascale=zeros(size(beta));
    //betascale=max(0,1-penbeta./abs(beta));
    arma::mat weightMat = weights * weights.t();

    // Rcpp::Rcout << "weightMat = \n" << weightMat << std::endl;

    const arma::mat& betaWeight = weightMat.submat(0, 0, p-1, p-1);
    arma::mat betaScale = betaWeight * -tlam(0);
    arma::mat absBeta = arma::abs(par.beta); 
    
    betaScale /= absBeta;
    betaScale += 1;
    betaScale.transform( [](double val) {return std::max(val, 0.0); } );

    // Rcpp::Rcout << "betaScale = \n" << betaScale << std::endl;

    double betaNorms = 0;

    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < p; j++) {
            double curVal = par.beta(i,j);
            if (curVal != 0) {
                curVal *= betaScale(i,j);
                // Rcpp::Rcout << "curVal = " << curVal << std::endl;
                par.beta(i,j) = curVal;
                betaNorms += std::abs(betaWeight(i,j)*curVal);
                // Rcpp::Rcout << "curBetaNorms = " << betaNorms << std::endl;

            }
        }
    }
    // Rcpp::Rcout << "NS betaNorms = \n" << betaNorms << std::endl;

    //weight beta
    //betaw = (wv(1:p)'*wv(1:p)).*beta;
    //betanorms=sum(abs(betaw(:)));
    //double betaNorm = betaWeight.copy().assign(par.beta, Functions.mult).assign(Functions.abs).zSum();

    /*
    thetanorms=0;
    for s=1:p
        for j=1:q
            tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
            tempvec=max(0,1-t(2)*(wv(s)*wv(p+j))/norm(tempvec))*tempvec;
            thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
            theta(Lsums(j)+1:Lsums(j+1),s)=tempvec(1:L(j));
        end
    end
    */
    double thetaNorms = 0;
    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            arma::subview_col<double> tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
            double thetaScale = std::max(0.0, 1 - tlam(1)*weightMat(i,p+j)/arma::norm(tempVec, 2));
            par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * thetaScale;
            thetaNorms += weightMat(i, p+j) * arma::norm(tempVec, 2);
        }
    }
    // Rcpp::Rcout << "NS thetaNorms = \n" << thetaNorms << std::endl;

    /*
    for r=1:q
        for j=1:q
            if r<j
                tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
                tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
                phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
                phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
            end
        end
    end
    */
    double phiNorms = 0;
    for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
        for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
            arma::subview<double> tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
            double phiScale = std::max(0.0, 1 - tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, "fro"));
            // Use the tempMat subview again to set the values (doesn't work with const)
            par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) = tempMat * phiScale;
            phiNorms += weightMat(p+i, p+j) * arma::norm(tempMat, "fro");
        }
    }

    // Rcpp::Rcout << "NS phiNorms = \n" << phiNorms << std::endl;

    arma::rowvec zWeight = weights.subvec(p+q, p+q+r-1).t();

    arma::mat gammaWeight = arma::mat(weightMat.submat(0, p+q, p-1, p+q+r-1));
    // gammaWeight.each_row( [this,&zWeight] (arma::rowvec& r) {
    // 			    r = zWeight + fitWeight % r;
    // 			  } );
    arma::mat gammaScale = gammaWeight * -tlam(3);
    arma::mat absGamma = arma::abs(par.gamma); 
    
    gammaScale /= absGamma;
    gammaScale += 1;
    gammaScale.transform( [](double val) {return std::max(val, 0.0); } );

    // Rcpp::Rcout << "betaScale = \n" << betaScale << std::endl;

    double gammaNorms = 0;

    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < r; j++) {
            double curVal = par.gamma(i,j);
            if (curVal != 0) {
                curVal *= gammaScale(i,j);
                // Rcpp::Rcout << "curVal = " << curVal << std::endl;
                par.gamma(i,j) = curVal;
                gammaNorms += std::abs(gammaWeight(i,j)*curVal);
                // Rcpp::Rcout << "curGammaNorms = " << gammaNorms << std::endl;

            }
        }
    }

    double psiNorms = 0;
    for (arma::uword i = 0; i < r; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
            double psiScale = std::max(0.0, 1 - tlam(4)*weightMat(p+j, p+q+i)/arma::norm(tempVec, 2));
            par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * psiScale;
            psiNorms += weightMat(p+j, p+q+i) * arma::norm(tempVec, 2);
        }
    }

    // Rcpp::Rcout << "Finished\n";

    pX = par.toMatrix1D();
    return lambda(0)*betaNorms + lambda(1)*thetaNorms + lambda(2)*phiNorms + lambda(3)*gammaNorms + lambda(4)*psiNorms;
}

// /**
//  * Calculate value of h(X) and proxOperator of h(X) at the same time for efficiency reasons.
//  *
//  * @param t positive parameter for prox operator
//  * @param X input vector
//  * @param Xout vector solution to prox_t(X)
//  * @return value of h(X)
//  */
// double CoxMGM::nonSmooth(double t, arma::vec& X, arma::vec& pX) {
//     double nonSmooth = 0;

//     arma::vec tlam = lambda * t;

//     double a = 4;

//     //Constructor copies and checks dimension
//     //par is a copy so we can update it
//     CoxMGMParams par(X, p, lsum, r);

//     // calcZWeights();

//     // Rcpp::Rcout << "NS par: \n" << par << std::endl;

//     //penbeta = t(1).*(wv(1:p)'*wv(1:p));
//     //betascale=zeros(size(beta));
//     //betascale=max(0,1-penbeta./abs(beta));
//     arma::mat weightMat = weights * weights.t();

//     double betaNorms = 0;
//     double thetaNorms = 0;
//     double phiNorms = 0;
//     double gammaNorms = 0;
//     double psiNorms = 0;
    
//     arma::mat betaNorm = arma::abs(par.beta);
//     arma::mat thetaNorm = arma::mat(q, p, arma::fill::zeros);
//     arma::mat phiNorm = arma::mat(q, q, arma::fill::zeros);
//     arma::mat gammaNorm = arma::abs(par.gamma);
//     arma::mat psiNorm = arma::mat(q, r, arma::fill::zeros);

//     for (arma::uword i = 0; i < p; i++) {
//         for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//             arma::subview_col<double> tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
//             thetaNorm(j,i) = arma::norm(tempVec, 2);
//         }
//     }

//     for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
//         for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
//             const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);

//     	    phiNorm(i,j) = arma::norm(tempMat, "fro");
// 	}
//     }
    
//     for (arma::uword i = 0; i < r; i++) {
// 	for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
// 	    arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
//             psiNorm(j,i) = arma::norm(tempVec, 2);
//         }
//     }

//     // L_1/2 thresholding

//     for (arma::uword i = 0; i < p; i++) {
//         for (arma::uword j = i+1; j < p; j++) {
// 	    double betaScale = 0.0;
// 	    if (betaNorm(i,j) > std::pow(0.75 * tlam(0) * weightMat(i,j), 2 / 3.0)) {
// 		double temp = std::acos(tlam(0) * weightMat(i,j) / 8.0 *
// 					std::pow(betaNorm(i,j) / 3, -1.5));
// 		betaScale = 2 / 3.0 * (1 + std::cos(2 * M_PI / 3.0 - 2 * temp / 3.0));
// 	    }
// 	    par.beta(i,j) *= betaScale;
// 	    betaNorms += weightMat(i,j) * std::sqrt(std::abs(par.beta(i,j)));
// 	}
//     }

    
//     for (arma::uword i = 0; i < p; i++) {
//         for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//             arma::subview_col<double> tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
// 	    double thetaScale = 0.0;
// 	    if (thetaNorm(j,i) > std::pow(0.75 * tlam(1) * weightMat(i,p+j), 2 / 3.0)) {
// 		double temp = std::acos(tlam(1) * weightMat(i,p+j) / 8.0 *
// 					std::pow(thetaNorm(j,i) / 3, -1.5));
// 		thetaScale = 2 / 3.0 * (1 + std::cos(2 * M_PI / 3.0 - 2 * temp / 3.0));
// 	    }
// 	    tempVec *= thetaScale;
// 	    thetaNorms += weightMat(i,p+j) * std::sqrt(arma::norm(tempVec, 2));
// 	}
//     }

//     for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
//         for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
//             arma::subview<double> tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
// 	    double phiScale = 0.0;
// 	    if (phiNorm(i,j) > std::pow(0.75 * tlam(2) * weightMat(p+i,p+j), 2 / 3.0)) {
// 		double temp = std::acos(tlam(2) * weightMat(p+i,p+j) / 8.0 *
// 					std::pow(phiNorm(i,j) / 3, -1.5));
// 		phiScale = 2 / 3.0 * (1 + std::cos(2 * M_PI / 3.0 - 2 * temp / 3.0));
// 	    }
// 	    tempMat *= phiScale;
// 	    phiNorms += weightMat(p+i,p+j) * std::sqrt(arma::norm(tempMat, 2));
// 	}
//     }


//     for (arma::uword i = 0; i < p; i++) {
//         for (arma::uword j = 0; j < r; j++) {
// 	    double gammaScale = 0.0;
// 	    if (gammaNorm(i,j) > std::pow(0.75 * tlam(3) * weightMat(i,p+q+j), 2 / 3.0)) {
// 		double temp = std::acos(tlam(3) * weightMat(i,p+q+j) / 8.0 *
// 					std::pow(gammaNorm(i,j) / 3, -1.5));
// 		gammaScale = 2 / 3.0 * (1 + std::cos(2 * M_PI / 3.0 - 2 * temp / 3.0));
// 	    }
// 	    par.gamma(i,j) *= gammaScale;
// 	    gammaNorms += weightMat(i,p+q+j) * std::sqrt(std::abs(par.gamma(i,j)));
// 	}
//     }


//     for (arma::uword i = 0; i < r; i++) {
//         for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//             arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
// 	    double psiScale = 0.0;
// 	    if (psiNorm(j,i) > std::pow(0.75 * tlam(4) * weightMat(p+j,p+q+i), 2 / 3.0)) {
// 		double temp = std::acos(tlam(4) * weightMat(p+j,p+q+i) / 8.0 *
// 					std::pow(psiNorm(j,i) / 3, -1.5));
// 		psiScale = 2 / 3.0 * (1 + std::cos(2 * M_PI / 3.0 - 2 * temp / 3.0));
// 	    }
// 	    tempVec *= psiScale;
// 	    psiNorms += weightMat(p+j,p+q+i) * std::sqrt(arma::norm(tempVec, 2));
// 	}
//     }


//     // Rcpp::Rcout << "old gamma: " << par.gamma << std::endl;

//     // Rcpp::Rcout << "Gamma norm: " << gammaNorm << std::endl;

//     //// MCP Penalty

//     // arma::mat betaX1 = arma::min(a * lambda(0) * weightMat.submat(0, 0, p-1, p-1) / betaNorm,
//     // 		       arma::max(a * (1 - tlam(0) * weightMat.submat(0, 0, p-1, p-1)
//     // 				      / betaNorm) / (a-1),
//     // 				 arma::zeros(arma::size(betaNorm))));
//     // arma::mat betaX2 = arma::max(a * lambda(0) * weightMat.submat(0, 0, p-1, p-1) / betaNorm,
//     // 		       arma::ones(arma::size(betaNorm)));

//     // arma::mat thetaX1 = arma::min(a * lambda(1) * weightMat.submat(p, 0, p+q-1, p-1)
//     // 				  / thetaNorm,
//     // 			arma::max(a * (1 - tlam(1) * weightMat.submat(p, 0, p+q-1, p-1)
//     // 				       / thetaNorm) / (a-1),
//     // 				  arma::zeros(arma::size(thetaNorm))));
//     // arma::mat thetaX2 = arma::max(a * lambda(1) * weightMat.submat(p, 0, p+q-1, p-1)
//     // 				  / thetaNorm,
//     // 			arma::ones(arma::size(thetaNorm)));

//     // arma::mat phiX1 = arma::min(a * lambda(2) * weightMat.submat(p, p, p+q-1, p+q-1)
//     // 				/ phiNorm,
//     // 		      arma::max(a * (1 - tlam(2) * weightMat.submat(p, p, p+q-1, p+q-1)
//     // 				     / phiNorm) / (a-1),
//     // 				arma::zeros(arma::size(phiNorm))));
//     // arma::mat phiX2 = arma::max(a * lambda(2) * weightMat.submat(p, p, p+q-1, p+q-1)
//     // 				/ phiNorm,
//     // 		      arma::ones(arma::size(phiNorm)));

//     // arma::mat gammaX1 = arma::min(a * lambda(3) * weightMat.submat(0, p+q, p-1, p+q+r-1)
//     // 				  / gammaNorm,
//     // 			arma::max(a * (1 - tlam(3) * weightMat.submat(0, p+q, p-1, p+q+r-1)
//     // 				       / gammaNorm) / (a-1),
//     // 				  arma::zeros(arma::size(gammaNorm))));
//     // arma::mat gammaX2 = arma::max(a * lambda(3) * weightMat.submat(0, p+q, p-1, p+q+r-1)
//     // 				  / gammaNorm,
//     // 			arma::ones(arma::size(gammaNorm)));

//     // arma::mat psiX1 = arma::min(a * lambda(4) * weightMat.submat(p, p+q, p+q-1, p+q+r-1)
//     // 				/ psiNorm,
//     // 		      arma::max(a * (1 - tlam(4) * weightMat.submat(p, p+q, p+q-1, p+q+r-1)
//     // 				     / psiNorm) / (a-1),
//     // 				arma::zeros(arma::size(psiNorm))));
//     // arma::mat psiX2 = arma::max(a * lambda(4) * weightMat.submat(p, p+q, p+q-1, p+q+r-1)
//     // 				/ psiNorm,
//     // 		      arma::ones(arma::size(psiNorm)));

//     // // Rcpp::Rcout << "GammaX1 solution: " << gammaX1 << std::endl;
//     // // Rcpp::Rcout << "GammaX2 solution: " << gammaX2 << std::endl;
    
//     // double betaNorms = 0;
//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = i+1; j < p; j++) {
//     // 	    // arma::vec xOpts = {0.0,
//     // 	    // 		       a * lambda(0) * weightMat(i,j),
//     // 	    // 		       betaX1(i,j) * betaNorm(i,j)};
//     // 	    // arma::vec zDists(3);
//     //         // zDists(0) = 0.5 * std::pow(betaNorm(i,j), 2);
//     // 	    // zDists(1) = 0.5 * std::pow(xOpts(1) - betaNorm(i,j), 2)
//     // 	    // 	+ tlam(0) * weightMat(i,j) * xOpts(1)
//     // 	    // 	- 0.5 * std::pow(xOpts(1), 2) / a;
//     // 	    // zDists(2) = 0.5 * std::pow(xOpts(2) - betaNorm(i,j), 2)
//     // 	    // 	+ tlam(0) * weightMat(i,j) * xOpts(2)
//     // 	    // 	- 0.5 * std::pow(xOpts(2), 2) / a;
//     // 	    // betaX1(i,j) = xOpts(zDists.index_min()) / betaNorm(i,j);

//     // 	    double hx1 = 0.5 * std::pow(betaX1(i,j) * betaNorm(i,j) - betaNorm(i,j), 2)
//     // 		+ tlam(0) * weightMat(i,j) * betaX1(i,j) * betaNorm(i,j)
//     // 		- 0.5 * std::pow(betaX1(i,j) * betaNorm(i,j), 2) / a;
	    
//     // 	    double hx2 = 0.5 * std::pow(betaX2(i,j) * betaNorm(i,j) - betaNorm(i,j), 2)
//     // 		+ 0.5 * a * std::pow(tlam(0) * weightMat(i,j), 2);

//     // 	    if (hx1 <= hx2) {
//     // 		par.beta(i,j) *= betaX1(i,j);
//     // 		betaNorms += lambda(0) * weightMat(i,j) * std::abs(par.beta(i,j))
//     // 		    - 0.5 * std::pow(par.beta(i,j), 2) / a;
//     // 	    } else {
//     // 		Rcpp::Rcout << "beta h1(x) > h2(x): (" << i << "," << j << ")\n";
//     // 		par.beta(i,j) *= betaX2(i,j);
//     // 		betaNorms += 0.5 * a * std::pow(lambda(0) * weightMat(i,j), 2);
//     // 	    }
//     //     }
//     // }

//     // double thetaNorms = 0;
//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     //         arma::subview_col<double> tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);

//     // 	    // arma::vec xOpts = {0.0,
//     // 	    // 		       a * lambda(1) * weightMat(i,p+j),
//     // 	    // 		       thetaX1(j,i) * thetaNorm(j,i)};
//     // 	    // arma::vec zDists(3);
//     //         // zDists(0) = 0.5 * std::pow(thetaNorm(j,i), 2);
//     // 	    // zDists(1) = 0.5 * std::pow(xOpts(1) - thetaNorm(j,i), 2)
//     // 	    // 	+ tlam(1) * weightMat(i,p+j) * xOpts(1)
//     // 	    // 	- 0.5 * std::pow(xOpts(1), 2) / a;
//     // 	    // zDists(2) = 0.5 * std::pow(xOpts(2) - thetaNorm(j,i), 2)
//     // 	    // 	+ tlam(1) * weightMat(i,p+j) * xOpts(2)
//     // 	    // 	- 0.5 * std::pow(xOpts(2), 2) / a;
//     // 	    // thetaX1(j,i) = xOpts(zDists.index_min()) / thetaNorm(j,i);

//     // 	    double hx1 = 0.5 * std::pow(thetaX1(j,i) * thetaNorm(j,i) - thetaNorm(j,i), 2)
//     // 		+ tlam(1) * weightMat(i,p+j) * thetaX1(j,i) * thetaNorm(j,i)
//     // 		- 0.5 * std::pow(thetaX1(j,i) * thetaNorm(j,i), 2) / a;
	    
//     // 	    double hx2 = 0.5 * std::pow(thetaX2(j,i) * thetaNorm(j,i) - thetaNorm(j,i), 2)
//     // 		+ 0.5 * a * std::pow(tlam(1) * weightMat(i,p+j), 2);

//     // 	    if (hx1 <= hx2) {
//     // 	        tempVec *= thetaX1(j,i);
//     // 		thetaNorms += lambda(1) * weightMat(i,p+j) * arma::norm(tempVec, 2)
//     // 		    - 0.5 * std::pow(arma::norm(tempVec, 2), 2) / a;
//     // 	    } else {
//     // 		tempVec *= thetaX2(j,i);
//     // 		thetaNorms += 0.5 * a * std::pow(lambda(1) * weightMat(i,p+j), 2);
//     // 	    }
//     // 	}
//     // }

//     // double phiNorms = 0;
//     // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
//     //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
//     //         arma::subview<double> tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);

//     // 	    // arma::vec xOpts = {0.0,
//     // 	    // 		       a * lambda(2) * weightMat(p+i,p+j),
//     // 	    // 		       phiX1(i,j) * phiNorm(i,j)};
//     // 	    // arma::vec zDists(3);
//     //         // zDists(0) = 0.5 * std::pow(phiNorm(i,j), 2);
//     // 	    // zDists(1) = 0.5 * std::pow(xOpts(1) - phiNorm(i,j), 2)
//     // 	    // 	+ tlam(2) * weightMat(p+i,p+j) * xOpts(1)
//     // 	    // 	- 0.5 * std::pow(xOpts(1), 2) / a;
//     // 	    // zDists(2) = 0.5 * std::pow(xOpts(2) - phiNorm(i,j), 2)
//     // 	    // 	+ tlam(2) * weightMat(p+i,p+j) * xOpts(2)
//     // 	    // 	- 0.5 * std::pow(xOpts(2), 2) / a;
//     // 	    // phiX1(i,j) = xOpts(zDists.index_min()) / phiNorm(i,j);

//     // 	    double hx1 = 0.5 * std::pow(phiX1(i,j) * phiNorm(i,j) - phiNorm(i,j), 2)
//     // 		+ tlam(2) * weightMat(p+i,p+j) * phiX1(i,j) * phiNorm(i,j)
//     // 		- 0.5 * std::pow(phiX1(i,j) * phiNorm(i,j), 2) / a;
	    
//     // 	    double hx2 = 0.5 * std::pow(phiX2(i,j) * phiNorm(i,j) - phiNorm(i,j), 2)
//     // 		+ 0.5 * a * std::pow(tlam(2) * weightMat(p+i,p+j), 2);

//     // 	    if (hx1 <= hx2) {
//     // 	        tempMat *= phiX1(i,j);
//     // 		phiNorms += lambda(2) * weightMat(p+i,p+j) * arma::norm(tempMat, "fro")
//     // 		    - 0.5 * std::pow(arma::norm(tempMat, "fro"), 2) / a;
//     // 	    } else {
//     // 		tempMat *= phiX2(i,j);
//     // 		phiNorms += 0.5 * a * std::pow(lambda(2) * weightMat(p+i,p+j), 2);
//     // 	    }
	    
//     // 	    phiNorm(i,j) = weightMat(p+i,p+j) * arma::norm(tempMat, "fro");
//     // 	}
//     // }

//     // double gammaNorms = 0;
//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < r; j++) {
//     // 	    // arma::vec xOpts = {0.0,
//     // 	    // 		       a * lambda(3) * weightMat(i,p+q+j),
//     // 	    // 		       gammaX1(i,j) * gammaNorm(i,j)};
//     // 	    // arma::vec zDists(3);
//     //         // zDists(0) = 0.5 * std::pow(gammaNorm(i,j), 2);
//     // 	    // zDists(1) = 0.5 * std::pow(xOpts(1) - gammaNorm(i,j), 2)
//     // 	    // 	+ tlam(3) * weightMat(i,p+q+j) * xOpts(1)
//     // 	    // 	- 0.5 * std::pow(xOpts(1), 2) / a;
//     // 	    // zDists(2) = 0.5 * std::pow(xOpts(2) - gammaNorm(i,j), 2)
//     // 	    // 	+ tlam(3) * weightMat(i,p+q+j) * xOpts(2)
//     // 	    // 	- 0.5 * std::pow(xOpts(2), 2) / a;
//     // 	    // gammaX1(i,j) = xOpts(zDists.index_min()) / gammaNorm(i,j);

//     // 	    double hx1 = 0.5 * std::pow(gammaX1(i,j) * gammaNorm(i,j) - gammaNorm(i,j), 2)
//     // 		+ tlam(3) * weightMat(i,p+q+j) * gammaX1(i,j) * gammaNorm(i,j)
//     // 		- 0.5 * std::pow(gammaX1(i,j) * gammaNorm(i,j), 2) / a;
	    
//     // 	    double hx2 = 0.5 * std::pow(gammaX2(i,j) * gammaNorm(i,j) - gammaNorm(i,j), 2)
//     // 		+ 0.5 * a * std::pow(tlam(3) * weightMat(i,p+q+j), 2);

//     // 	    if (hx1 <= hx2) {
//     // 		par.gamma(i,j) *= gammaX1(i,j);
//     // 		gammaNorms += lambda(3) * weightMat(i,p+q+j) * std::abs(par.gamma(i,j))
//     // 		    - 0.5 * std::pow(par.gamma(i,j), 2) / a;
//     // 	    } else {
//     // 		Rcpp::Rcout << "gamma h1(x) > h2(x): (" << i << "," << j << ")\n";
//     // 		par.gamma(i,j) *= gammaX2(i,j);
//     // 		gammaNorms += 0.5 * a * std::pow(lambda(3) * weightMat(i,p+q+j), 2);
//     // 	    }
//     //     }
//     // }

//     //  double psiNorms = 0;
//     // for (arma::uword i = 0; i < r; i++) {
//     // 	for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     // 	    arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
//     // 	    // arma::vec xOpts = {0.0,
//     // 	    // 		       a * lambda(4) * weightMat(p+j,p+q+i),
//     // 	    // 		       psiX1(j,i) * psiNorm(j,i)};
//     // 	    // arma::vec zDists(3);
//     //         // zDists(0) = 0.5 * std::pow(psiNorm(j,i), 2);
//     // 	    // zDists(1) = 0.5 * std::pow(xOpts(1) - psiNorm(j,i), 2)
//     // 	    // 	+ tlam(4) * weightMat(p+j,p+q+i) * xOpts(1)
//     // 	    // 	- 0.5 * std::pow(xOpts(1), 2) / a;
//     // 	    // zDists(2) = 0.5 * std::pow(xOpts(2) - psiNorm(j,i), 2)
//     // 	    // 	+ tlam(4) * weightMat(p+j,p+q+i) * xOpts(2)
//     // 	    // 	- 0.5 * std::pow(xOpts(2), 2) / a;
//     // 	    // psiX1(j,i) = xOpts(zDists.index_min()) / psiNorm(j,i);

//     // 	    double hx1 = 0.5 * std::pow(psiX1(j,i) * psiNorm(j,i) - psiNorm(j,i), 2)
//     // 		+ tlam(4) * weightMat(p+j,p+q+i) * psiX1(j,i) * psiNorm(j,i)
//     // 		- 0.5 * std::pow(psiX1(j,i) * psiNorm(j,i), 2) / a;
	    
//     // 	    double hx2 = 0.5 * std::pow(psiX2(j,i) * psiNorm(j,i) - psiNorm(j,i), 2)
//     // 		+ 0.5 * a * std::pow(tlam(4) * weightMat(p+j,p+q+i), 2);

//     // 	    if (hx1 <= hx2) {
//     // 	        tempVec *= psiX1(j,i);
//     // 		psiNorms += lambda(4) * weightMat(p+j,p+q+i) * arma::norm(tempVec, 2)
//     // 		    - 0.5 * std::pow(arma::norm(tempVec, 2), 2) / a;
//     // 	    } else {
//     // 		tempVec *= psiX2(j,i);
//     // 		psiNorms += 0.5 * a * std::pow(lambda(4) * weightMat(p+j,p+q+i), 2);
//     // 	    }
//     // 	}
//     // }

//     // Rcpp::Rcout << "new beta: \n" << par.beta << std::endl;
//     // Rcpp::Rcout << "nonSmooth loss for beta: " << betaNorms << std::endl;

//     // Rcpp::Rcout << "new theta: \n" << par.theta << std::endl;
//     // Rcpp::Rcout << "nonSmooth loss for theta: " << thetaNorms << std::endl;

//     // Rcpp::Rcout << "new phi: \n" << par.phi << std::endl;
//     // Rcpp::Rcout << "nonSmooth loss for phi: " << phiNorms << std::endl;

//     // Rcpp::Rcout << "new gamma: \n" << par.gamma << std::endl;
//     // Rcpp::Rcout << "nonSmooth loss for gamma: " << gammaNorms << std::endl;

//     // Rcpp::Rcout << "new psi: \n" << par.psi << std::endl;
//     // Rcpp::Rcout << "nonSmooth loss for psi: " << psiNorms << std::endl;


//     // pX = par.toMatrix1D();
//     // return betaNorms + thetaNorms + phiNorms + gammaNorms + psiNorms;

    
    
//     // Rcpp::Rcout << "weightMat = \n" << weightMat << std::endl;

//     // const arma::mat& betaWeight = weightMat.submat(0, 0, p-1, p-1);
//     // arma::mat betaScale = betaWeight * tlam(0);
    
//     // betaScale /= arma::abs(par.beta);
//     // // betaScale += 1;
//     // betaScale.transform( [&a](double val) {
//     // 			     if (val >= 0.5) {
//     // 				 return std::max(1 - val, 0.0);
//     // 			     } else if (val >= (1/a)) {
//     // 				 return (a-1) / (a-2) * std::max(1 - a * val / (a-1), 0.0);
//     // 			     } else {
//     // 				 return 1.0;
//     // 			     }
//     // 			 } );

//     // par.beta = par.beta % betaScale;

//     // double betaNorms = 0;
//     // arma::mat betaLambda = tlam(0) * weightMat.submat(0,0,p-1,p-1);
//     // arma::mat absBeta = arma::abs(par.beta);
//     // arma::mat betaPenalty = arma::mat(p,p,arma::fill::zeros);

//     // for (arma::uword i = 0; i < p; i++) {
//     // 	for (arma::uword j = i+1; j < p; j++) {
//     // 	    if (absBeta(i,j) < betaLambda(i,j)) {
//     // 		betaPenalty(i,j) = betaLambda(i,j) * absBeta(i,j);
//     // 	    } else if (absBeta(i,j) < a * betaLambda(i,j)) {
//     // 		betaPenalty(i,j) = (2 * a * betaLambda(i,j) * absBeta(i,j) - std::pow(absBeta(i,j), 2) - std::pow(betaLambda(i,j), 2)) / (2 * (a-1));
//     // 	    } else {
//     // 		betaPenalty(i,j) = std::pow(betaLambda(i,j), 2) * (a + 1) / 2;
//     // 	    }
//     // 	}
//     // }

    
//     // // Rcpp::Rcout << "Beta penalties: \n" << betaPenalty;

//     // betaNorms = arma::accu(betaPenalty);

//     // // Rcpp::Rcout << "NSV betaNorms = " << betaNorms << std::endl;

//     // /*
//     // thetanorms=0;
//     // for s=1:p
//     //     for j=1:q
//     //         tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
//     //         thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
//     //     end
//     // end
//     // */
//     // double thetaNorms = 0;
//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     //         const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);

//     // 	    double thetaScale = tlam(1)*weightMat(i,p+j)/arma::norm(tempVec, 2);
//     // 	    if (thetaScale > 0.5) {
//     // 		thetaScale = std::max(1-thetaScale, 0.0);
//     // 	    } else if (thetaScale > (1/a)) {
//     // 		thetaScale = (a-1) / (a-2) * std::max(1 - a * thetaScale / (a-1), 0.0);
//     // 	    } else {
//     // 		thetaScale = 1;
//     // 	    }
//     //         par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * thetaScale;
	    
	    
//     // 	    double wLam = weightMat(i, p+j) * tlam(1);
//     // 	    double normTheta = arma::norm(tempVec, 2);
//     // 	    if (normTheta < wLam) {
//     // 		thetaNorms += wLam * normTheta;
//     // 	    } else if (normTheta < a * wLam) {
//     // 	        thetaNorms += (2 * a * wLam * normTheta - std::pow(normTheta, 2) - std::pow(wLam, 2)) / (2 * (a-1));
//     // 	    } else {
//     // 	        thetaNorms += std::pow(wLam, 2) * (a + 1) / 2;
//     // 	    }
//     //     }
//     // }
//     // // Rcpp::Rcout << "NSV thetaNorms = " << thetaNorms << std::endl;

//     // /*
//     // for r=1:q
//     //     for j=1:q
//     //         if r<j
//     //             tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
//     //             tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
//     //             phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
//     //             phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
//     //         end
//     //     end
//     // end
//     // */
//     // double phiNorms = 0;
//     // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
//     //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
//     //         const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);

//     // 	    double phiScale = tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, "fro");
//     // 	    if (phiScale > 0.5) {
//     // 		phiScale = std::max(1-phiScale, 0.0);
//     // 	    } else if (phiScale > (1/a)) {
//     // 		phiScale = (a-1) / (a-2) * std::max(1 - a * phiScale / (a-1), 0.0);
//     // 	    } else {
//     // 		phiScale = 1;
//     // 	    }
//     //         // Use the tempMat subview again to set the values (doesn't work with const)
//     //         par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) = tempMat * phiScale;
	    
//     // 	    double wLam = weightMat(p+i, p+j) * tlam(2);
//     // 	    double normPhi = arma::norm(tempMat, "fro");
//     // 	    if (normPhi < wLam) {
//     // 		phiNorms += wLam * normPhi;
//     // 	    } else if (normPhi < a * wLam) {
//     // 	        phiNorms += (2 * a * wLam * normPhi - std::pow(normPhi, 2) - std::pow(wLam, 2)) / (2 * (a-1));
//     // 	    } else {
//     // 	        phiNorms += std::pow(wLam, 2) * (a + 1) / 2;
//     // 	    }
//     //     }
//     // }

    
//     // // double gammaNorms = arma::accu(weightMat.submat(0, p+q, p-1, p+q+r-1) % abs(par.gamma));

//     // const arma::mat& gammaWeight = weightMat.submat(0, p+q, p-1, p+q+r-1);

//     // // Rcpp::Rcout << "GammaWeight:\n" << gammaWeight;
    
//     // arma::mat gammaScale = gammaWeight * tlam(3);
//     // arma::mat absGamma = arma::abs(par.gamma); 
    
//     // gammaScale /= absGamma;
//     // // gammaScale += 1;

//     // // Rcpp::Rcout << "GammaScale:\n" << gammaScale;
    
//     // gammaScale.transform( [&a](double val) {
//     // 			      if (val > 0.5) {
//     // 				  return std::max(1 - val, 0.0);
//     // 			      } else if (val > (1/a)) {
//     // 				  return (a-1) / (a-2) * std::max(1 - a * val / (a-1), 0.0);
//     // 			      } else {
//     // 				  return 1.0;
//     // 			      }
//     // 			  } );

//     // // Rcpp::Rcout << "Thresholded GammaScale:\n" << gammaScale;

//     // par.gamma = par.gamma % gammaScale;

//     // double gammaNorms = 0;
//     // arma::mat gammaLambda = tlam(3) * weightMat.submat(0, p+q, p-1, p+q+r-1);
//     // absGamma = arma::abs(par.gamma);
//     // arma::mat gammaPenalty = arma::mat(p,r,arma::fill::zeros);

//     // for (arma::uword i = 0; i < p; i++) {
//     // 	for (arma::uword j = 0; j < r; j++) {
//     // 	    if (absGamma(i,j) < gammaLambda(i,j)) {
//     // 		gammaPenalty(i,j) = gammaLambda(i,j) * absGamma(i,j);
//     // 	    } else if (absGamma(i,j) < a * gammaLambda(i,j)) {
//     // 		gammaPenalty(i,j) = (2 * a * gammaLambda(i,j) * absGamma(i,j) - std::pow(absGamma(i,j), 2) - std::pow(gammaLambda(i,j), 2)) / (2 * (a-1));
//     // 	    } else {
//     // 		gammaPenalty(i,j) = std::pow(gammaLambda(i,j), 2) * (a + 1) / 2;
//     // 	    }
//     // 	}
//     // }

    
//     // // Rcpp::Rcout << "Gamma penalties: \n" << gammaPenalty;

//     // gammaNorms = arma::accu(gammaPenalty);
    
//     // double psiNorms = 0;
//     // for (arma::uword i = 0; i < r; i++) {
//     //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     //         arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);

//     // 	    double psiScale = tlam(4)*weightMat(p+j, p+q+i)/arma::norm(tempVec, 2);
//     // 	    if (psiScale > 0.5) {
//     // 		psiScale = std::max(1-psiScale, 0.0);
//     // 	    } else if (psiScale > (1/a)) {
//     // 		psiScale = (a-1) / (a-2) * std::max(1 - a * psiScale / (a-1), 0.0);
//     // 	    } else {
//     // 		psiScale = 1;
//     // 	    }
//     //         par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * psiScale;
	    
//     // 	    double wLam = weightMat(p+j, p+q+i) * tlam(4);
//     // 	    double normPsi = arma::norm(tempVec, 2);
//     // 	    if (normPsi < wLam) {
//     // 		psiNorms += wLam * normPsi;
//     // 	    } else if (normPsi < a * wLam) {
//     // 	        psiNorms += (2 * a * wLam * normPsi - std::pow(normPsi, 2) - std::pow(wLam, 2)) / (2 * (a-1));
//     // 	    } else {
//     // 	        psiNorms += std::pow(wLam, 2) * (a + 1) / 2;
//     // 	    }
	    
//     //     }
//     // }

//     // pX = par.toMatrix1D();
//     // return betaNorms + thetaNorms + phiNorms + gammaNorms + psiNorms;

//     // const arma::mat& betaWeight = weightMat.submat(0, 0, p-1, p-1);
//     // arma::mat betaScale = betaWeight * -tlam(0);
//     // arma::mat absBeta = arma::abs(par.beta); 
    
//     // betaScale /= absBeta;
//     // betaScale += 1;
//     // betaScale.transform( [](double val) {return std::max(val, 0.0); } );

//     // // Rcpp::Rcout << "betaScale = \n" << betaScale << std::endl;

//     // double betaNorms = 0;

//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < p; j++) {
//     //         double curVal = par.beta(i,j);
//     //         if (curVal != 0) {
//     //             curVal *= betaScale(i,j);
//     //             // Rcpp::Rcout << "curVal = " << curVal << std::endl;
//     //             par.beta(i,j) = curVal;
//     //             betaNorms += std::abs(betaWeight(i,j)*curVal);
//     //             // Rcpp::Rcout << "curBetaNorms = " << betaNorms << std::endl;

//     //         }
//     //     }
//     // }
//     // // Rcpp::Rcout << "NS betaNorms = \n" << betaNorms << std::endl;

//     // //weight beta
//     // //betaw = (wv(1:p)'*wv(1:p)).*beta;
//     // //betanorms=sum(abs(betaw(:)));
//     // //double betaNorm = betaWeight.copy().assign(par.beta, Functions.mult).assign(Functions.abs).zSum();

//     // /*
//     // thetanorms=0;
//     // for s=1:p
//     //     for j=1:q
//     //         tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
//     //         tempvec=max(0,1-t(2)*(wv(s)*wv(p+j))/norm(tempvec))*tempvec;
//     //         thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
//     //         theta(Lsums(j)+1:Lsums(j+1),s)=tempvec(1:L(j));
//     //     end
//     // end
//     // */
//     // double thetaNorms = 0;
//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     //         arma::subview_col<double> tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
//     //         double thetaScale = std::max(0.0, 1 - tlam(1)*weightMat(i,p+j)/arma::norm(tempVec, 2));
//     //         par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * thetaScale;
//     //         thetaNorms += weightMat(i, p+j) * arma::norm(tempVec, 2);
//     //     }
//     // }
//     // // Rcpp::Rcout << "NS thetaNorms = \n" << thetaNorms << std::endl;

//     // /*
//     // for r=1:q
//     //     for j=1:q
//     //         if r<j
//     //             tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
//     //             tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
//     //             phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
//     //             phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
//     //         end
//     //     end
//     // end
//     // */
//     // double phiNorms = 0;
//     // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
//     //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
//     //         arma::subview<double> tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
//     //         double phiScale = std::max(0.0, 1 - tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, "fro"));
//     //         // Use the tempMat subview again to set the values (doesn't work with const)
//     //         par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) = tempMat * phiScale;
//     //         phiNorms += weightMat(p+i, p+j) * arma::norm(tempMat, "fro");
//     //     }
//     // }

//     // // Rcpp::Rcout << "NS phiNorms = \n" << phiNorms << std::endl;

//     // const arma::mat& gammaWeight = weightMat.submat(0, p+q, p-1, p+q+r-1);
//     // arma::mat gammaScale = gammaWeight * -tlam(3);
//     // arma::mat absGamma = arma::abs(par.gamma); 
    
//     // gammaScale /= absGamma;
//     // gammaScale += 1;
//     // gammaScale.transform( [](double val) {return std::max(val, 0.0); } );

//     // // Rcpp::Rcout << "betaScale = \n" << betaScale << std::endl;

//     // double gammaNorms = 0;

//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < r; j++) {
//     //         double curVal = par.gamma(i,j);
//     //         if (curVal != 0) {
//     //             curVal *= gammaScale(i,j);
//     //             // Rcpp::Rcout << "curVal = " << curVal << std::endl;
//     //             par.gamma(i,j) = curVal;
//     //             gammaNorms += std::abs(gammaWeight(i,j)*curVal);
//     //             // Rcpp::Rcout << "curGammaNorms = " << gammaNorms << std::endl;

//     //         }
//     //     }
//     // }

//     // double psiNorms = 0;
//     // for (arma::uword i = 0; i < r; i++) {
//     //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     //         arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
//     //         double psiScale = std::max(0.0, 1 - tlam(4)*weightMat(p+j, p+q+i)/arma::norm(tempVec, 2));
//     //         par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * psiScale;
//     //         psiNorms += weightMat(p+j, p+q+i) * arma::norm(tempVec, 2);
//     //     }
//     // }

//     // const arma::mat& betaWeight = weightMat.submat(0, 0, p-1, p-1);
//     // arma::mat betaScale = betaWeight * -tlam(0);
//     // arma::mat absBeta = arma::abs(par.beta); 
    
//     // betaScale /= absBeta;
//     // betaScale += 1;
//     // betaScale.transform( [](double val) { return (val >= 0.0) ? 1.0 : 0.0; } );

//     // // Rcpp::Rcout << "betaScale = \n" << betaScale << std::endl;

//     // double betaNorms = 0;

//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < p; j++) {
//     //         double curVal = par.beta(i,j);
//     //         if (curVal != 0) {
//     //             curVal *= betaScale(i,j);
//     //             // Rcpp::Rcout << "curVal = " << curVal << std::endl;
//     //             par.beta(i,j) = curVal;
//     //             betaNorms += std::abs(betaWeight(i,j) * (curVal!=0) );
//     //             // Rcpp::Rcout << "curBetaNorms = " << betaNorms << std::endl;

//     //         }
//     //     }
//     // }
//     // // Rcpp::Rcout << "NS betaNorms = \n" << betaNorms << std::endl;

//     // //weight beta
//     // //betaw = (wv(1:p)'*wv(1:p)).*beta;
//     // //betanorms=sum(abs(betaw(:)));
//     // //double betaNorm = betaWeight.copy().assign(par.beta, Functions.mult).assign(Functions.abs).zSum();

//     // /*
//     // thetanorms=0;
//     // for s=1:p
//     //     for j=1:q
//     //         tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
//     //         tempvec=max(0,1-t(2)*(wv(s)*wv(p+j))/norm(tempvec))*tempvec;
//     //         thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
//     //         theta(Lsums(j)+1:Lsums(j+1),s)=tempvec(1:L(j));
//     //     end
//     // end
//     // */
//     // double thetaNorms = 0;
//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     //         arma::subview_col<double> tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
//     //         double thetaScale = (1 - tlam(1)*weightMat(i,p+j)/arma::norm(tempVec, 2) >= 0) ? 1.0 : 0.0;
//     //         par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * thetaScale;
//     //         thetaNorms += weightMat(i, p+j) * (arma::norm(tempVec, 2) != 0);
//     //     }
//     // }
//     // // Rcpp::Rcout << "NS thetaNorms = \n" << thetaNorms << std::endl;

//     // /*
//     // for r=1:q
//     //     for j=1:q
//     //         if r<j
//     //             tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
//     //             tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
//     //             phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
//     //             phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
//     //         end
//     //     end
//     // end
//     // */
//     // double phiNorms = 0;
//     // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
//     //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
//     //         arma::subview<double> tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
//     //         double phiScale = (1 - tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, "fro") >= 0) ? 1.0 : 0.0;
//     //         // Use the tempMat subview again to set the values (doesn't work with const)
//     //         par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) = tempMat * phiScale;
//     //         phiNorms += weightMat(p+i, p+j) * (arma::norm(tempMat, "fro") != 0);
//     //     }
//     // }

//     // // Rcpp::Rcout << "NS phiNorms = \n" << phiNorms << std::endl;

//     // const arma::mat& gammaWeight = weightMat.submat(0, p+q, p-1, p+q+r-1);
//     // arma::mat gammaScale = gammaWeight * -tlam(3);
//     // arma::mat absGamma = arma::abs(par.gamma); 
    
//     // gammaScale /= absGamma;
//     // gammaScale += 1;
//     // gammaScale.transform( [](double val) {return (val >= 0) ? 1.0 : 0.0; } );

//     // // Rcpp::Rcout << "betaScale = \n" << betaScale << std::endl;

//     // double gammaNorms = 0;

//     // for (arma::uword i = 0; i < p; i++) {
//     //     for (arma::uword j = 0; j < r; j++) {
//     //         double curVal = par.gamma(i,j);
//     //         if (curVal != 0) {
//     //             curVal *= gammaScale(i,j);
//     //             // Rcpp::Rcout << "curVal = " << curVal << std::endl;
//     //             par.gamma(i,j) = curVal;
//     //             gammaNorms += std::abs(gammaWeight(i,j)*(curVal!=0));
//     //             // Rcpp::Rcout << "curGammaNorms = " << gammaNorms << std::endl;

//     //         }
//     //     }
//     // }

//     // double psiNorms = 0;
//     // for (arma::uword i = 0; i < r; i++) {
//     //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
//     //         arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
//     //         double psiScale = (1 - tlam(4)*weightMat(p+j, p+q+i)/arma::norm(tempVec, 2) >= 0) ? 1.0 : 0.0;
//     //         par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * psiScale;
//     //         psiNorms += weightMat(p+j, p+q+i) * (arma::norm(tempVec, 2) != 0);
//     //     }
//     // }

//     pX = par.toMatrix1D();
//     return lambda(0)*betaNorms + lambda(1)*thetaNorms + lambda(2)*phiNorms + lambda(3)*gammaNorms + lambda(4)*psiNorms;
// }

/**
 * Calculate value of h(X)
 *
 * @param X input vector
 * @return value of h(X)
 */
double CoxMGM::nonSmoothValue(arma::vec& parIn) {
    //Dimension checked in constructor
    //par is a copy so we can update it
    CoxMGMParams par(parIn, p, lsum, r);

    double a = 5;

    // Rcpp::Rcout << "nonSmoothValue called\n";

    // calcZWeights();

    //penbeta = t(1).*(wv(1:p)'*wv(1:p));
    //betascale=zeros(size(beta));
    //betascale=max(0,1-penbeta./abs(beta));
    arma::mat weightMat = weights * weights.t();

    double betaNorms = 0;
    arma::mat betaLambda = lambda(0) * weightMat.submat(0,0,p-1,p-1);
    arma::mat absBeta = arma::abs(par.beta);
    arma::mat betaPenalty = arma::mat(p,p,arma::fill::zeros);

    for (arma::uword i = 0; i < p; i++) {
	for (arma::uword j = i+1; j < p; j++) {
	    if (absBeta(i,j) < betaLambda(i,j)) {
		betaPenalty(i,j) = betaLambda(i,j) * absBeta(i,j);
	    } else if (absBeta(i,j) < a * betaLambda(i,j)) {
		betaPenalty(i,j) = (2 * a * betaLambda(i,j) * absBeta(i,j) - std::pow(absBeta(i,j), 2) - std::pow(betaLambda(i,j), 2)) / (2 * (a-1));
	    } else {
		betaPenalty(i,j) = std::pow(betaLambda(i,j), 2) * (a + 1) / 2;
	    }
	}
    }

    betaNorms = arma::accu(betaPenalty);

    // Rcpp::Rcout << "NSV betaNorms = " << betaNorms << std::endl;

    /*
    thetanorms=0;
    for s=1:p
        for j=1:q
            tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
            thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
        end
    end
    */
    double thetaNorms = 0;
    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);

	    double wLam = weightMat(i, p+j) * lambda(1);
	    double normTheta = arma::norm(tempVec, 2);
	    if (normTheta < wLam) {
		thetaNorms += wLam * normTheta;
	    } else if (normTheta < a * wLam) {
	        thetaNorms += (2 * a * wLam * normTheta - std::pow(normTheta, 2) - std::pow(wLam, 2)) / (2 * (a-1));
	    } else {
	        thetaNorms += std::pow(wLam, 2) * (a + 1) / 2;
	    }
        }
    }
    // Rcpp::Rcout << "NSV thetaNorms = " << thetaNorms << std::endl;

    /*
    for r=1:q
        for j=1:q
            if r<j
                tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
                tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
                phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
                phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
            end
        end
    end
    */
    double phiNorms = 0;
    for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
        for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
            const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
	    double wLam = weightMat(p+i, p+j) * lambda(2);
	    double normPhi = arma::norm(tempMat, "fro");
	    if (normPhi < wLam) {
		phiNorms += wLam * normPhi;
	    } else if (normPhi < a * wLam) {
	        phiNorms += (2 * a * wLam * normPhi - std::pow(normPhi, 2) - std::pow(wLam, 2)) / (2 * (a-1));
	    } else {
	        phiNorms += std::pow(wLam, 2) * (a + 1) / 2;
	    }
        }
    }

    
    // double gammaNorms = arma::accu(weightMat.submat(0, p+q, p-1, p+q+r-1) % abs(par.gamma));

    double gammaNorms = 0;
    arma::mat gammaLambda = lambda(3) * weightMat.submat(0, p+q, p-1, p+q+r-1);
    arma::mat absGamma = arma::abs(par.gamma);
    arma::mat gammaPenalty = arma::mat(p,r,arma::fill::zeros);

    for (arma::uword i = 0; i < p; i++) {
	for (arma::uword j = 0; j < r; j++) {
	    if (absGamma(i,j) < gammaLambda(i,j)) {
		gammaPenalty(i,j) = gammaLambda(i,j) * absGamma(i,j);
	    } else if (absGamma(i,j) < a * gammaLambda(i,j)) {
		gammaPenalty(i,j) = (2 * a * gammaLambda(i,j) * absGamma(i,j) - std::pow(absGamma(i,j), 2) - std::pow(gammaLambda(i,j), 2)) / (2 * (a-1));
	    } else {
		gammaPenalty(i,j) = std::pow(gammaLambda(i,j), 2) * (a + 1) / 2;
	    }
	}
    }

    gammaNorms = arma::accu(gammaPenalty);
    
    double psiNorms = 0;
    for (arma::uword i = 0; i < r; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
	    double wLam = weightMat(p+j, p+q+i) * lambda(4);
	    double normPsi = arma::norm(tempVec, 2);
	    if (normPsi < wLam) {
		psiNorms += wLam * normPsi;
	    } else if (normPsi < a * wLam) {
	        psiNorms += (2 * a * wLam * normPsi - std::pow(normPsi, 2) - std::pow(wLam, 2)) / (2 * (a-1));
	    } else {
	        psiNorms += std::pow(wLam, 2) * (a + 1) / 2;
	    }
	    
        }
    }

    // Rcpp::Rcout << "NSV phiNorms = " << phiNorms << std::endl;

    return betaNorms + thetaNorms + phiNorms + gammaNorms + psiNorms;

    // //weight beta
    // //betaw = (wv(1:p)'*wv(1:p)).*abs(beta);
    // //betanorms=sum(betaw(:));
    // double betaNorms = arma::accu(arma::mat(weightMat.submat(0, 0, p-1, p-1)) % arma::abs(par.beta));

    // // Rcpp::Rcout << "NSV betaNorms = " << betaNorms << std::endl;

    // /*
    // thetanorms=0;
    // for s=1:p
    //     for j=1:q
    //         tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
    //         thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
    //     end
    // end
    // */
    // double thetaNorms = 0;
    // for (arma::uword i = 0; i < p; i++) {
    //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
    //         const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
    //         thetaNorms += weightMat(i, p+j) * arma::norm(tempVec, 2);
    //     }
    // }
    // // Rcpp::Rcout << "NSV thetaNorms = " << thetaNorms << std::endl;

    // /*
    // for r=1:q
    //     for j=1:q
    //         if r<j
    //             tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
    //             tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
    //             phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
    //             phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
    //         end
    //     end
    // end
    // */
    // double phiNorms = 0;
    // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
    //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
    //         const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
    //         phiNorms += weightMat(p+i, p+j) * arma::norm(tempMat, "fro");
    //     }
    // }

    
    // double gammaNorms = arma::accu(weightMat.submat(0, p+q, p-1, p+q+r-1) % abs(par.gamma));
    
    // double psiNorms = 0;
    // for (arma::uword i = 0; i < r; i++) {
    //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
    //         arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
    //         psiNorms += weightMat(p+j, p+q+i) * arma::norm(tempVec, 2);
    //     }
    // }

    // // Rcpp::Rcout << "NSV phiNorms = " << phiNorms << std::endl;

    // return lambda(0)*betaNorms + lambda(1)*thetaNorms + lambda(2)*phiNorms + lambda(3)*gammaNorms + lambda(4)*psiNorms;
}

/**
 * Gradient of the pseudolikelihood
 *
 * @param parIn
 * @return
 */
arma::vec CoxMGM::smoothGradient(arma::vec& parIn) {
    int n = xDat.n_rows;
    CoxMGMParams grad;

    // Rcpp::Rcout << "Running smoothGradient...\n";

    // grad.gamma = arma::zeros(p, r);
    // grad.psi = arma::zeros(lsum, r);
    // grad.xid = arma::zeros(r);
    // grad.alpha3 = arma::zeros(r);

    CoxMGMParams par(parIn, p, lsum, r);

    // for (arma::uword j = 0; j < q; j++) {
    //  	par.psi.row(lcumsum[j]).fill(0);
    // }

    // Rcpp::Rcout << "par: \n" << par << std::endl;

    arma::mat eta = xDat * par.gamma + dDat * par.psi;
    eta.each_row() += par.alpha3.t();
    // eta.each_col( [](arma::vec& c) { c -= arma::mean(c); } );

    
    // arma::mat coxgrad(arma::size(eta), arma::fill::zeros);
    // arma::mat diagHess(arma::size(eta), arma::fill::zeros);

    // arma::vec oldCoxloss = coxGradHess(eta, coxgrad, diagHess);

    // coxWeights = -diagHess;
    // diagHess.replace(0, -1e-10);
    // zDat = eta - coxgrad / diagHess;

    // Rcpp::Rcout << "eta: " << eta << std::endl;
    // Rcpp::Rcout << "Z: " << zDat << std::endl;
    // Rcpp::Rcout << "coxWeights: " << coxWeights << std::endl;

    par.beta = arma::symmatu(par.beta);
    par.beta.diag(0).zeros();
    // Rcpp::Rcout << "par.beta: \n" << par.beta << std::endl;

    for (arma::uword i = 0; i < q; i++) {
        par.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    par.phi = arma::symmatu(par.phi);
    // Rcpp::Rcout << "par.phi: \n" << par.phi << std::endl;

    arma::mat divBetaD = arma::diagmat(arma::vec(p, arma::fill::ones) / par.betad);
    // Rcpp::Rcout << "divBetaD: \n" << divBetaD << std::endl;

    //Xbeta=X*beta*diag(1./betad);
    //Dtheta=D*theta*diag(1./betad);
    arma::mat xBeta = xDat * par.beta * divBetaD;
    arma::mat dTheta = dDat * par.theta * divBetaD;
    arma::mat zGamma = wzDat * par.gamma.t() * divBetaD;
    // Rcpp::Rcout << "xBeta: \n" << xBeta << std::endl;
    // Rcpp::Rcout << "dTheta: \n" << dTheta << std::endl;

    // Rcpp::Rcout << "zGamma: \n" << zGamma << std::endl;

    // zGamma.fill(0);

    // for (arma::uword m = 0; m < r; m++) {
    // 	zGamma += arma::diagmat(coxWeights.col(m)) * zDat.col(m) * par.gamma.col(m).t() * divBetaD;
    // }
    // // zGamma = zGamma * divBetaD;
    
    // Rcpp::Rcout << "zGamma2: \n" << zGamma << std::endl;

    //res=Xbeta-X+e*alpha1'+Dtheta;
    //wxprod=X*(theta')+D*phi+e*alpha2';
    arma::mat negLoss(n, xDat.n_cols);

    arma::mat wxProd = xDat*par.theta.t() + dDat*par.phi + wzDat * par.psi.t();
    // Rcpp::Rcout << "wxProd1: \n" << wxProd << std::endl;

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < p; j++) {
            negLoss(i,j) = xBeta(i,j) - xDat(i,j) + par.alpha1(j) + dTheta(i,j) + zGamma(i,j);
        }
        for (arma::uword j = 0; j < dDat.n_cols; j++) {
            wxProd(i, j) += par.alpha2(j);
        }
    }
    // Rcpp::Rcout << "negLoss: \n" << negLoss << std::endl;
    // Rcpp::Rcout << "wxProd2: \n" << wxProd << std::endl;

    //gradbeta=X'*(res);
    grad.beta = xDat.t() * negLoss;
    // Rcpp::Rcout << "grad.beta1: \n" << grad.beta << std::endl;


    //gradbeta=gradbeta-diag(diag(gradbeta)); % zero out diag
    //gradbeta=tril(gradbeta)'+triu(gradbeta);
    grad.beta = arma::trimatl(grad.beta, -1).t() + arma::trimatu(grad.beta, 1);
    // Rcpp::Rcout << "grad.beta2: \n" << grad.beta << std::endl;

    // Rcpp::Rcout << "negLoss: \n" << negLoss << std::endl;

    //gradalpha1=diag(betad)*sum(res,0)';
    grad.alpha1 = arma::diagmat(par.betad) * arma::sum(negLoss, 0).t();

    // Rcpp::Rcout << "grad.alpha1: \n" << grad.alpha1 << std::endl;

    //gradtheta=D'*(res);
    grad.theta = dDat.t() * negLoss;
    // Rcpp::Rcout << "grad.theta: \n" << grad.theta << std::endl;

    // grad.gamma = (wzDat.t() * negLoss).t();
    // grad.gamma = (5/2.0) * ((coxWeights % zDat).t() * (negLoss / (4 + arma::square(resid)))).t();
    // grad.gamma.each_row( [this] (arma::rowvec& r) { r %= (fitWeight / (1 + fitWeight)); } );
    // grad.gamma.each_row( [this] (arma::rowvec& r) { r %= fitWeight; } );
    // grad.gamma = arma::mat(p, r, arma::fill::zeros);

    for (arma::uword i = 0; i < yDat.n_cols; i++) {
        arma::subview<double> wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));

        // does this need to be done in log space??
	wxTemp = arma::exp(wxTemp.each_col() - logsumexp(arma::mat(wxTemp)));
        // wxTemp = arma::exp(wxTemp);
        // arma::vec invDenom = arma::vec(n, arma::fill::ones) / arma::sum(wxTemp, 1);
        // wxTemp = arma::diagmat(invDenom) * wxTemp;

        for (arma::uword k = 0; k < n; k++) {
            //wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))=wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))-1;
            wxTemp(k, (arma::uword) yDat(k,i)-1) -= 1;
        }
    }
    // Rcpp::Rcout << "wxProd: \n" << wxProd << std::endl;

    //gradalpha2=sum(wxprod,0)';
    grad.alpha2 = arma::sum(wxProd, 0).t();
    // Rcpp::Rcout << "grad.alpha2: \n" << grad.alpha2 << std::endl;

    //gradw=X'*wxprod;
    arma::mat gradW = xDat.t() * wxProd;
    // Rcpp::Rcout << "gradW: \n" << gradW << std::endl;

    //gradtheta=gradtheta+gradw';
    grad.theta += gradW.t();
    // Rcpp::Rcout << "grad.theta: \n" << grad.theta << std::endl;

    //gradphi=D'*wxprod;
    grad.phi = dDat.t() * wxProd;
    // Rcpp::Rcout << "grad.phi1: \n" << grad.phi << std::endl;

    // grad.psi = (wzDat.t() * wxProd).t();
    // grad.psi = (5/2.0) * ((coxWeights % zDat).t() * (wxProd / (4 + arma::square(catResid)))).t();
    // grad.psi.each_row( [this] (arma::rowvec& r) { r %= (fitWeight / (1 + fitWeight)); } );
    // grad.psi.each_row( [this] (arma::rowvec& r) { r %= fitWeight; } );
    // grad.psi = arma::mat(lsum, r, arma::fill::zeros);

    //zero out gradphi diagonal
    //for r=1:q
    //gradphi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    for (arma::uword i = 0; i < q; i++) {
        grad.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }
    //gradphi=tril(gradphi)'+triu(gradphi);
    grad.phi = arma::trimatl(grad.phi, 0).t() + arma::trimatu(grad.phi, 0);

    // arma::mat expGamma2(p,r,arma::fill::zeros);
    // arma::mat expPsi2(q,r,arma::fill::zeros);

    arma::mat coxresid = eta - zDat;
    grad.gamma = xDat.t() * (coxWeights % coxresid);
    grad.psi = dDat.t() * (coxWeights % coxresid);
    
    // for (arma::uword i = 0; i < zDat.n_cols; i++) {
    // 	// coxresid.col(i) += par.alpha3(i);
    // 	// grad.gamma.col(i) += (1 / (1+fitWeight(i))) * xDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);
    //     // grad.psi.col(i) += (1 / (1+fitWeight(i))) * dDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);

    // 	grad.gamma.col(i) += xDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);

    // 	// arma::mat psiGradTemp = arma::diagmat(coxWeights.col(i) % coxresid.col(i)) * dDat;

    // 	// for (arma::uword j = 0; j < q; j++) {
    // 	//     arma::subview<double> wxTemp = psiGradTemp.submat(0, lcumsum[j], n-1, lcumsum[j+1]-1);
    // 	//     wxTemp.each_row( [this] (arma::rowvec& r) { r -= arma::mean(r); } );
    // 	// }
	
    //     // grad.psi.col(i) += arma::mean(psiGradTemp,0);
    // 	grad.psi.col(i) += dDat.t() * arma::diagmat(coxWeights.col(i)) * coxresid.col(i);
	
    // 	// grad.gamma.col(i) /= arma::sum(coxWeights.col(i));
    // 	// grad.psi.col(i) /= arma::sum(coxWeights.col(i));

    // 	// if (i == 0) {
    // 	// arma::mat gamma2 = arma::square(arma::diagmat(coxWeights.col(i) % coxresid.col(i)) * xDat + fitWeight(i) * arma::diagmat(coxWeights.col(i) % zDat.col(i)) * negLoss);
    // 	//     // Rcpp::Rcout << "Gamma^2:\n" << gamma2 << std::endl;

    // 	// Rcpp::Rcout << "cont->surv("<< i << ") gamma grad:\n"
    // 	// 	    << arma::diagmat(coxWeights.col(i) % coxresid.col(i)) * xDat;
	
    // 	// Rcpp::Rcout << "surv("<< i << ")->cont gamma grad:\n"
    // 	// 	    << arma::diagmat(coxWeights.col(i) % zDat.col(i)) * negLoss;

    // 	// expGamma2.col(i) = arma::mean(gamma2,0);

    // 	// Rcpp::Rcout << "E[Gamma^2]:\n" << arma::mean(gamma2,0) << std::endl;

    // 	// Rcpp::Rcout << "sqrt(E[Gamma^2]):\n" << arma::sqrt(arma::mean(gamma2,0)) << std::endl;

    // 	// arma::mat temp = arma::square(psiGradTemp + fitWeight(i) * arma::diagmat(coxWeights.col(i) % zDat.col(i)) * wxProd);
    // 	// arma::mat psi2(n, q, arma::fill::zeros);

    // 	// Rcpp::Rcout << "disc->surv("<< i << ") psi grad:\n"
    // 	// 	    << 2 * (1-fitWeight(i)) * psiGradTemp;
	
    // 	// Rcpp::Rcout << "surv("<< i << ")->disc psi grad:\n"
    // 	// 	    << 2 * fitWeight(i) * arma::diagmat(coxWeights.col(i) % zDat.col(i)) * wxProd;


    // 	// for (arma::uword j = 0; j < q; j++) {
    // 	//     arma::mat wxTemp = temp.submat(0, lcumsum[j], n-1, lcumsum[j+1]-1);
    // 	//     psi2.col(j) = arma::sum(wxTemp, 1);
    // 	//     // if (j == 0) {
    // 	//     // 	Rcpp::Rcout << wxTemp << std::endl;
    // 	//     // 	Rcpp::Rcout << arma::sum(wxTemp, 1) << std::endl;
    // 	//     // }
    // 	// }
    // 	// // Rcpp::Rcout << "Psi^2:\n" << psi2 << std::endl;

    // 	// expPsi2.col(i) = arma::mean(psi2,0);
	
    // 	// Rcpp::Rcout << "E[Psi^2]:\n" << arma::mean(psi2,0) << std::endl;
	
    // 	// Rcpp::Rcout << "sqrt(E[Psi^2]):\n" << arma::sqrt(arma::mean(psi2,0)) << std::endl;
    // 	// }
    // }

    grad.alpha3 = arma::sum(coxWeights % coxresid, 0).t();

    // Rcpp::Rcout << "sqrt(E[Gamma^2]):\n" << arma::sqrt(expGamma2)/2 << std::endl;

    // Rcpp::Rcout << "sqrt(E[Psi^2]):\n" << arma::sqrt(expPsi2)/2 << std::endl;

    // for (arma::uword j = 0; j < q; j++) {
    // 	grad.psi.row(lcumsum[j]).fill(0);
    // }
    
    /*
    for s=1:p
        gradbetad(s)=-n/(2*betad(s))+1/2*norm(res(:,s))^2-res(:,s)'*(Xbeta(:,s)+Dtheta(:,s));
    end
        */
    grad.betad = arma::vec(xDat.n_cols);
    for(arma::uword i = 0; i < p; i++){
        // Rcpp::Rcout << "BETA NORM (i = " << i << "):\n" << std::pow(arma::norm(negLoss.col(i), 2), 2) << std::endl;
        // Rcpp::Rcout << "BETA SCALAR (i = " << i << "):\n" << arma::as_scalar(negLoss.col(i).t() * (xBeta.col(i) + dTheta.col(i))) << std::endl;
        grad.betad(i) = -n / (2.0 * par.betad(i)) + std::pow(arma::norm(negLoss.col(i), 2), 2) / 2.0 -
	    arma::as_scalar(negLoss.col(i).t() * (xBeta.col(i) + dTheta.col(i) + zGamma.col(i) + par.alpha1(i) / par.betad(i)));
    }

    grad.alpha1 /= (double) n;
    grad.alpha2 /= (double) n;
    grad.alpha3 /= (double) n; // arma::sum(coxWeights,0).t();
    grad.betad /= (double) n;
    grad.beta /= (double) n;
    grad.theta /= (double) n;
    grad.phi /= (double) n;
    grad.gamma /= (double) n;
    grad.psi /= (double) n;

    if (pDummy != 0) {
	grad.alpha1.fill(0);
	grad.betad.fill(0);
	grad.beta.fill(0);
	grad.theta.fill(0);
	grad.gamma.fill(0);
    }

    if (qDummy != 0) {
	grad.alpha2.fill(0);
	grad.theta.fill(0);
	grad.phi.fill(0);
	grad.psi.fill(0);
    }

    // Rcpp::Rcout << "Finished\n";

    return grad.toMatrix1D();
}

/**
 * A proximal operator is the solution to this optimization problem:
 *     prox_t(x) = argmin_z \frac{1}{2t} \|x-z\|^2_2 + h(x)
 *
 * @param t positive parameter for prox operator
 * @param X input vector
 * @return vector solution to prox_t(X)
 */
arma::vec CoxMGM::proximalOperator(double t, arma::vec& X) {
    if (t <= 0)
        throw std::invalid_argument("t must be positive: " + std::to_string(t));

    double a = 5;

    // Rcpp::Rcout << "proximalOperator called\n";
    // Rcpp::Rcout << "Lambda:\n" << lambda.t() << std::endl;
    arma::vec tlam = lambda * t;

    //Constructor copies and checks dimension
    //par is a copy so we can update it
    CoxMGMParams par(X, p, lsum, r);

    // calcZWeights();

    //penbeta = t(1).*(wv(1:p)'*wv(1:p));
    //betascale=zeros(size(beta));
    //betascale=max(0,1-penbeta./abs(beta));
    arma::mat weightMat = weights * weights.t();

    const arma::mat& betaWeight = weightMat.submat(0, 0, p-1, p-1);
    arma::mat betaScale = betaWeight * tlam(0);
    
    betaScale /= arma::abs(par.beta);
    // betaScale += 1;
    betaScale.transform( [&a](double val) {
			     if (val > 0.5) {
				 return std::max(1 - val, 0.0);
			     } else if (val > (1/a)) {
				 return (a-1) / (a-2) * std::max(1 - a * val / (a-1), 0.0);
			     } else {
				 return 1.0;
			     }
			 } );

    par.beta = par.beta % betaScale;

    //weight beta
    //betaw = (wv(1:p)'*wv(1:p)).*beta;
    //betanorms=sum(abs(betaw(:)));

    /*
    thetanorms=0;
    for s=1:p
        for j=1:q
            tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
            tempvec=max(0,1-t(2)*(wv(s)*wv(p+j))/norm(tempvec))*tempvec;
            thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
            theta(Lsums(j)+1:Lsums(j+1),s)=tempvec(1:L(j));
        end
    end
    */
    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
            double thetaScale = tlam(1)*weightMat(i,p+j)/arma::norm(tempVec, 2);
	    if (thetaScale > 0.5) {
		thetaScale = std::max(1-thetaScale, 0.0);
	    } else if (thetaScale > (1/a)) {
		thetaScale = (a-1) / (a-2) * std::max(1 - a * thetaScale / (a-1), 0.0);
	    } else {
		thetaScale = 1;
	    }
            par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * thetaScale;
        }
    }

    /*
    for r=1:q
        for j=1:q
            if r<j
                tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
                tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
                phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
                phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
            end
        end
    end
    */
    for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
        for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
            const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
            double phiScale = tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, "fro");
	    if (phiScale > 0.5) {
		phiScale = std::max(1-phiScale, 0.0);
	    } else if (phiScale > (1/a)) {
		phiScale = (a-1) / (a-2) * std::max(1 - a * phiScale / (a-1), 0.0);
	    } else {
		phiScale = 1;
	    }
            // Use the tempMat subview again to set the values (doesn't work with const)
            par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) = tempMat * phiScale;
        }
    }


    const arma::mat& gammaWeight = weightMat.submat(0, p+q, p-1, p+q+r-1);

    // Rcpp::Rcout << "GammaWeight:\n" << gammaWeight;
    
    arma::mat gammaScale = gammaWeight * tlam(3);
    arma::mat absGamma = arma::abs(par.gamma); 
    
    gammaScale /= absGamma;
    // gammaScale += 1;

    // Rcpp::Rcout << "GammaScale:\n" << gammaScale;
    
    gammaScale.transform( [&a](double val) {
			      if (val > 0.5) {
				  return std::max(1 - val, 0.0);
			      } else if (val > (1/a)) {
				  return (a-1) / (a-2) * std::max(1 - a * val / (a-1), 0.0);
			      } else {
				  return 1.0;
			      }
			  } );

    // Rcpp::Rcout << "Thresholded GammaScale:\n" << gammaScale;

    par.gamma = par.gamma % gammaScale;


    for (arma::uword i = 0; i < r; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
            double psiScale = tlam(4)*weightMat(p+j, p+q+i)/arma::norm(tempVec, 2);
	    if (psiScale > 0.5) {
		psiScale = std::max(1-psiScale, 0.0);
	    } else if (psiScale > (1/a)) {
		psiScale = (a-1) / (a-2) * std::max(1 - a * psiScale / (a-1), 0.0);
	    } else {
		psiScale = 1;
	    }
            par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * psiScale;
        }
    }



    // calcZWeights();

    // //penbeta = t(1).*(wv(1:p)'*wv(1:p));
    // //betascale=zeros(size(beta));
    // //betascale=max(0,1-penbeta./abs(beta));
    // arma::mat weightMat = weights * weights.t();

    // const arma::mat& betaWeight = weightMat.submat(0, 0, p-1, p-1);
    // arma::mat betaScale = betaWeight * -tlam(0);
    
    // betaScale /= arma::abs(par.beta);
    // betaScale += 1;
    // betaScale.transform( [](double val) { return std::max(val, 0.0); } );

    // par.beta = par.beta % betaScale;

    // //weight beta
    // //betaw = (wv(1:p)'*wv(1:p)).*beta;
    // //betanorms=sum(abs(betaw(:)));

    // /*
    // thetanorms=0;
    // for s=1:p
    //     for j=1:q
    //         tempvec=theta(Lsums(j)+1:Lsums(j+1),s);
    //         tempvec=max(0,1-t(2)*(wv(s)*wv(p+j))/norm(tempvec))*tempvec;
    //         thetanorms=thetanorms+(wv(s)*wv(p+j))*norm(tempvec);
    //         theta(Lsums(j)+1:Lsums(j+1),s)=tempvec(1:L(j));
    //     end
    // end
    // */
    // for (arma::uword i = 0; i < p; i++) {
    //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
    //         const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
    //         double thetaScale = std::max(0.0, 1 - tlam(1)*weightMat(i,p+j)/arma::norm(tempVec, 2));
    //         par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * thetaScale;
    //     }
    // }

    // /*
    // for r=1:q
    //     for j=1:q
    //         if r<j
    //             tempmat=phi(Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1));
    //             tempmat=max(0,1-t(3)*(wv(p+r)*wv(p+j))/norm(tempmat))*tempmat; % Lj by 2*Lr
    //             phinorms=phinorms+(wv(p+r)*wv(p+j))*norm(tempmat,'fro');
    //             phi( Lsums(r)+1:Lsums(r+1),Lsums(j)+1:Lsums(j+1) )=tempmat;
    //         end
    //     end
    // end
    // */
    // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
    //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
    //         const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
    //         double phiScale = std::max(0.0, 1 - tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, "fro"));
    //         // Use the tempMat subview again to set the values (doesn't work with const)
    //         par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) = tempMat * phiScale;
    //     }
    // }


    // const arma::mat& gammaWeight = weightMat.submat(0, p+q, p-1, p+q+r-1);

    // // Rcpp::Rcout << "GammaWeight:\n" << gammaWeight;
    
    // arma::mat gammaScale = gammaWeight * -tlam(3);
    // arma::mat absGamma = arma::abs(par.gamma); 
    
    // gammaScale /= absGamma;
    // gammaScale += 1;

    // // Rcpp::Rcout << "GammaScale:\n" << gammaScale;
    
    // gammaScale.transform( [](double val) {return std::max(val, 0.0); } );

    // // Rcpp::Rcout << "Thresholded GammaScale:\n" << gammaScale;

    // par.gamma = par.gamma % gammaScale;


    // for (arma::uword i = 0; i < r; i++) {
    //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
    //         arma::subview_col<double> tempVec = par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
    //         double psiScale = std::max(0.0, 1 - tlam(4)*weightMat(p+j, p+q+i)/arma::norm(tempVec, 2));
    //         par.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) = tempVec * psiScale;
    //     }
    // }

    return par.toMatrix1D();
}


void CoxMGM::iterUpdate(arma::vec& parIn) {
    CoxMGMParams par(parIn, p, lsum, r);

    // Rcpp::Rcout << "Running iterUpdate...\n";

    arma::mat eta = xDat * par.gamma + dDat * par.psi;
    eta.each_row() -= arma::mean(eta, 0);
    
    coxgrad.zeros();
    diagHess.zeros();
    
    oldCoxloss = coxGradHess(eta, coxgrad, diagHess);

    coxWeights = -diagHess;
    diagHess.replace(0, -1e-10);
    zDat = eta - coxgrad / diagHess;
    wzDat = coxWeights % zDat;
    // wzDat = eta + coxgrad; // coxWeights % zDat; 

    // fitWeight = ((nullCoxloss-oldCoxloss) / nullCoxloss).t();

    // Rcpp::Rcout << "% deviance explained: " << 100*fitWeight << std::endl;
    
    // fitWeight = arma::rowvec(r, arma::fill::zeros);
    
    // arma::mat wZDat = coxWeights % zDat;
    // arma::rowvec weightSums = arma::sum(coxWeights,0);
    // arma::rowvec zDatBar = arma::sum(wZDat) / weightSums;
    // arma::rowvec zDat2Bar = arma::sum(wZDat % zDat) / weightSums;
    // arma::rowvec zDatWSd = arma::sqrt(weightSums / (weightSums-1) * (zDat2Bar - arma::square(zDatBar)));
    
    // weights.subvec(p+q, p+q+r-1) = weightSums / ((double) n) * zDatWSd;

    // arma::rowvec eventRate = arma::mean(censMat, 0);
    // weights.subvec(p+q, p+q+r-1) = arma::sqrt(eventRate % (1-eventRate));

    // par.beta = arma::symmatu(par.beta);
    // par.beta.diag(0).zeros();

    // for (arma::uword i = 0; i < q; i++) {
    //     par.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    // }

    // par.phi = arma::symmatu(par.phi);

    // arma::mat divBetaD = arma::diagmat(arma::vec(p, arma::fill::ones) / par.betad);

    // //Xbeta=X*beta*diag(1./betad);
    // //Dtheta=D*theta*diag(1./betad);
    // arma::mat xBeta = xDat * par.beta * divBetaD;
    // arma::mat dTheta = dDat * par.theta * divBetaD;
    // arma::mat zGamma = (coxWeights % zDat) * par.gamma.t() * divBetaD;

    // // Rcpp::Rcout << "xBeta: \n" << xBeta << std::endl;
    // // Rcpp::Rcout << "dTheta: \n" << dTheta << std::endl;
    // // Rcpp::Rcout << "zGamma: \n" << dTheta << std::endl;

    // // Squared loss
    // //sqloss=-n/2*sum(log(betad))+...
    // //.5*norm((X-e*alpha1'-Xbeta-Dtheta)*diag(sqrt(betad)),'fro')^2;
    // arma::mat tempLoss(n, xDat.n_cols);

    // //wxprod=X*(theta')+D*phi+e*alpha2';
    // arma::mat wxProd = xDat * par.theta.t() + dDat * par.phi + (zDat % coxWeights) * par.psi.t();

    // for (arma::uword i = 0; i < n; i++) {
    //     for (arma::uword j = 0; j < xDat.n_cols; j++) {
    //         tempLoss(i, j) = xDat(i,j) - par.alpha1(j) - xBeta(i,j) - dTheta(i,j) - zGamma(i,j);
    //     }
    //     for (arma::uword j = 0; j < dDat.n_cols; j++) {
    //         wxProd(i,j) = wxProd(i,j) + par.alpha2(j);
    //     }
    // }

    // resid = tempLoss;

    // for (arma::uword i = 0; i < yDat.n_cols; i++) {
    //     arma::subview<double> wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));

    //     // does this need to be done in log space??
    //     wxTemp = arma::exp(wxTemp);
    //     arma::vec invDenom = arma::vec(n, arma::fill::ones) / arma::sum(wxTemp, 1);
    //     wxTemp = arma::diagmat(invDenom) * wxTemp;

    //     for (arma::uword k = 0; k < n; k++) {
    //         //wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))=wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))-1;
    //         wxTemp(k, (arma::uword) yDat(k,i)-1) -= 1;
    //     }
    // }

    // catResid = wxProd;

    // zDat = arma::mat(zDat);
    // zDat.each_row( [&par](arma::rowvec& r) { r -= par.alpha3.t(); } );

    // for (arma::uword i = 0; i < r; i++) {
    // 	zDat.col(i) /= arma::stddev(coxWeights.col(i) % zDat.col(i));
    // }
    // weights.subvec(p+q, p+q+r-1) = arma::vec(r, arma::fill::ones);

    // Rcpp::Rcout << par << std::endl;

    // wzSd = arma::rowvec(arma::stddev(coxWeights % zDat, 0, 0));

    // weights.subvec(p+q, p+q+r-1) = arma::stddev(coxWeights % zDat, 0, 0);

    // weights.subvec(p+q, p+q+r-1) = arma::stddev(coxWeights % zDat, 0, 0) / 2;
      // (2 * arma::sqrt(arma::mean(arma::conv_to<arma::mat>::from(censMat), 0)));
      // arma::vec(r, arma::fill::ones);

    

    // arma::rowvec zWeights = weights.subvec(p+q, p+q+r-1).t();

    // Rcpp::Rcout << "Z weights:\n" << weights.subvec(p+q, p+q+r-1).t() << std::endl;

    // arma::mat weightMat = weights * weights.t();

    // arma::mat gammaWeights = arma::mat(weightMat.submat(0, p+q, p-1, p+q+r-1));

    // gammaWeights.each_row( [this] (arma::rowvec& r) { r %= fitWeight; } );
    // gammaWeights.each_row( [&zWeights] (arma::rowvec& r) { r += zWeights; } );

    // Rcpp::Rcout << "Gamma weights:\n"
    // 		<< gammaWeights << std::endl;

    // arma::mat psiWeights = arma::mat(weightMat.submat(p, p+q, p+q-1, p+q+r-1));

    // psiWeights.each_row( [this] (arma::rowvec& r) { r %= fitWeight; } );
    // psiWeights.each_row( [&zWeights] (arma::rowvec& r) { r += zWeights; } );

    // Rcpp::Rcout << "Psi weights:\n"
    // 		<< psiWeights << std::endl;

    // arma::vec curParams = params.toMatrix1D();

    // smoothGradient(curParams);

    // calcZWeights();

    // Rcpp::Rcout << "Finished\n";
}

std::string CoxMGM::printParameters(arma::vec& X) {
    std::stringstream ss;
    CoxMGMParams par(X, p, lsum, r);
    ss << par;
    return ss.str();
}

/**
 *  Learn CoxMGM traditional way with objective function tolerance. Recommended for inference applications that need
 *  accurate pseudolikelihood
 *
 * @param epsilon tolerance in change of objective function
 * @param iterLimit iteration limit
 */
void CoxMGM::learn(double epsilon, int iterLimit) {
    ProximalGradient pg = ProximalGradient();
    arma::vec curParams = params.toMatrix1D();
    arma::vec newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, epsilon, iterLimit);
    params = CoxMGMParams(newParams, p, lsum, r);
}

/**
 *  Learn CoxMGM using edge convergence using default 3 iterations of no edge changes. Recommended when we only care about
 *  edge existence.
 *
 * @param iterLimit
 */
void CoxMGM::learnEdges(int iterLimit) {
    ProximalGradient pg(0.5, 0.9, true);
    arma::vec curParams = params.toMatrix1D();
    arma::vec newParams;
    if (timeout != -1)
        newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 0.0, iterLimit, timeout);
    else
        newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 0.0, iterLimit);

    params = CoxMGMParams(newParams, p, lsum, r);

    timePerIter = pg.timePerIter;
    iterCount = pg.iterComplete;
}   

/**
 *  Learn CoxMGM using edge convergence using edgeChangeTol (see ProximalGradient for documentation). Recommended when we only care about
 *  edge existence.
 *
 * @param iterLimit
 * @param edgeChangeTol
 */
void CoxMGM::learnEdges(int iterLimit, int edgeChangeTol){
    ProximalGradient pg(0.5, 0.9, true);
    arma::vec curParams = params.toMatrix1D();
    pg.setEdgeChangeTol(edgeChangeTol);
    arma::vec newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 0.0, iterLimit);
    params = CoxMGMParams(newParams, p, lsum, r);
}

/**
 * Converts CoxMGM to matrix of doubles. uses 2-norm to combine c-d edge parameters into single value and f-norm for
 * d-d edge parameters.
 *
 * @return
 */
arma::mat CoxMGM::adjMatFromCoxMGM() {
    arma::mat outMat(p+q+r-pDummy-qDummy, p+q+r-pDummy-qDummy, arma::fill::zeros);

    if (p - pDummy > 0) 
	outMat(0, 0, arma::size(p, p)) = params.beta + params.beta.t();

    for (arma::uword i = 0; i < p-pDummy; i++) {
        for (arma::uword j = 0; j < q-qDummy; j++) {
            double val = arma::norm(params.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1), 2);
            outMat(i, p+j) = val;
            outMat(p+j, i) = val;
        }
    }

    for (arma::uword i = 0; i < p-pDummy; i++) {
        for (arma::uword j = 0; j < r; j++) {
            double val = params.gamma(i,j);
            outMat(i, p+q-qDummy+j) = val;
            outMat(p+q-qDummy+j, i) = val;
        }
    }

    for (arma::uword i = 0; i < q-qDummy; i++) {
        for (arma::uword j = i+1; j < q-qDummy; j++) {
            double val = arma::norm(params.phi(lcumsum[i], lcumsum[j], arma::size(l[i], l[j])), "fro");
            outMat(p-pDummy+i, p-pDummy+j) = val;
            outMat(p-pDummy+j, p-pDummy+i) = val;
        }
    }

    for (arma::uword i = 0; i < q-qDummy; i++) {
        for (arma::uword j = 0; j < r; j++) {
            double val = arma::norm(params.psi.col(j).subvec(lcumsum[i], lcumsum[i+1]-1), 2);
            outMat(p-pDummy+i, p-pDummy+q+j) = val;
            outMat(p-pDummy+q+j, p-pDummy+i) = val;
        }
    }

    //order the adjmat to be the same as the original DataSet variable ordering
    arma::uvec varMap(p+q+r-pDummy-qDummy);
    for(arma::uword i = 0; i < p+q+r-pDummy-qDummy; i++){
        varMap(i) = std::distance(variables.begin(), std::find(variables.begin(), variables.end(), initVariables[i]));
    }
    outMat = outMat.submat(varMap, varMap);

    return outMat;
}


/**
 * Converts CoxMGM object to Graph object with edges if edge parameters are non-zero. Loses all edge param information
 *
 * @return
 */
EdgeListGraph CoxMGM::graphFromCoxMGM() {
    EdgeListGraph g(variables);

    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = i+1; j < p; j++) {
            double v1 = params.beta(i,j);

            if (std::abs(v1) > 0) {
                if(!g.isAdjacentTo(variables[i], variables[j]) &&
		   !(variables[i] == dummyVar || variables[j] == dummyVar)) {
                    g.addUndirectedEdge(variables[i], variables[j]);
                }
            }
        }
    }

    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < q; j++) {
            double v1 = arma::accu(arma::abs(params.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1)));

            if (v1 > 0) {
                if(!g.isAdjacentTo(variables[i], variables[p+j]) &&
		   !(variables[i] == dummyVar || variables[p+j] == dummyVar)) {
                    g.addUndirectedEdge(variables[i], variables[p+j]);
                }
            }
        }
    }

    for (arma::uword i = 0; i < q; i++) {
        for (arma::uword j = i+1; j < q; j++) {
            double v1 = arma::accu(arma::abs(params.phi(lcumsum[i], lcumsum[j], arma::size(l[i], l[j]))));

            if (v1 > 0) {
                if(!g.isAdjacentTo(variables[p+i], variables[p+j]) &&
		   !(variables[p+i] == dummyVar || variables[p+j] == dummyVar)) {
                    g.addUndirectedEdge(variables[p+i], variables[p+j]);
                }
            }
        }
    }

    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < r; j++) {
            double v1 = params.gamma(i,j);

            if (std::abs(v1) > 0) {
                if(!g.isAdjacentTo(variables[i], variables[p+q+j]) &&
		   !(variables[i] == dummyVar)) {
                    g.addUndirectedEdge(variables[i], variables[p+q+j]);
                }
            }
        }
    }

    for (arma::uword i = 0; i < r; i++) {
        for (arma::uword j = 0; j < q; j++) {
            double v1 = arma::accu(arma::abs(params.psi.col(i).subvec(lcumsum[j], lcumsum[j+1]-1)));

            if (v1 > 0) {
                if(!g.isAdjacentTo(variables[p+q+i], variables[p+j]) &&
		   !(variables[p+j] == dummyVar)) {
                    g.addUndirectedEdge(variables[p+q+i], variables[p+j]);
                }
            }
        }
    }

    // Set algorithm and type
    std::ostringstream alg;
    alg << "CoxMGM";

    g.setAlgorithm(alg.str());
    g.setGraphType("undirected");

    g.setHyperParam("lambda", lambda);

    return g;
}

/**
 * Simple search command for GraphSearch implementation. Uses default edge convergence, 1000 iter limit.
 *
 * @return
 */
EdgeListGraph CoxMGM::search() {
    auto start = std::chrono::high_resolution_clock::now();
    // fitWeight = arma::rowvec(r, arma::fill::zeros);
    // learn(1e-5, 500);
    // fitWeight = arma::rowvec(r, arma::fill::ones);
    if (verbose) {
	RcppThread::Rcout << "  Learning CoxMGM for lambda = { "
			  << lambda[0] << " " << lambda[1] << " "
			  << lambda[2] << " " << lambda[3] << " "
			  << lambda[4] << " }\n";
    }
    learn(1e-5, 500);
    // elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
    double elapsedTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-start).count();
    if (verbose) {
	double factor = (elapsedTime < 10) ? std::pow(10, 2 - std::ceil(std::log10(std::abs(elapsedTime)))) : 1.0;
	elapsedTime = std::round(elapsedTime * factor) / factor;
        Rcpp::Rcout << "  CoxMGM Elapsed Time =  " << elapsedTime << " s" << std::endl;
    }
    return graphFromCoxMGM();
}


std::vector<EdgeListGraph> CoxMGM::searchPath(std::vector<double> lambdas,
					      arma::vec& loglik,
					      arma::vec& nParams) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<EdgeListGraph> pathGraphs;
    std::sort(lambdas.begin(), lambdas.end(), std::greater<double>());
    for (int i = 0; i < lambdas.size(); i++) {
	if (verbose) RcppThread::Rcout << "  Learning CoxMGM for lambda = " << lambdas[i] << "\n";
	std::vector<double> lambda = { lambdas[i], lambdas[i], lambdas[i],
				       lambdas[i], lambdas[i] };
	setLambda(lambda);
	// if (verbose) RcppThread::Rcout << "    Fitting first step...";
        // fitWeight = arma::rowvec(r, arma::fill::zeros);
	learn(1e-5, 500);
	// if (verbose) RcppThread::Rcout << "done.\n    Fitting second step...";
	// fitWeight = arma::rowvec(r, arma::fill::ones);
	// learn(1e-5, 500);
	// if (verbose) RcppThread::Rcout << "done.\n";
	pathGraphs.push_back(graphFromCoxMGM());

	// Rcpp::NumericVector lambdaVec = { lambdas[i], lambdas[i], lambdas[i],
	// 				  lambdas[i], lambdas[i] };
	
	// pathGraphs[i].setHyperParam("lambda", lambdaVec);
	// Rcpp::Rcout << params << std::endl;
	
	arma::vec par(params.toMatrix1D());
	loglik(i) = -n * smoothValue(par);
	nParams(i) = arma::accu(par!=0);
	
	// nParams(i) = 2 * arma::accu(params.beta!=0);
	// nParams(i) += 2 * arma::accu(params.theta!=0);
	// nParams(i) += 2 * arma::accu(params.phi!=0);
	// nParams(i) += 2 * arma::accu(params.gamma!=0);
	// nParams(i) += 2 * arma::accu(params.psi!=0);
	// nParams(i) += params.betad.n_elem;
	// nParams(i) += params.alpha1.n_elem;
	// nParams(i) += params.alpha2.n_elem;
	// nParams(i) += params.alpha3.n_elem;

	// Rcpp::Rcout << params << std::endl;
	
	RcppThread::checkUserInterrupt();
    }
    if (verbose) RcppThread::Rcout << std::endl;;
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
    double elapsedTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-start).count();
    if (verbose) {
	double factor = (elapsedTime < 10) ? std::pow(10, 2 - std::ceil(std::log10(std::abs(elapsedTime)))) : 1.0;
	elapsedTime = std::round(elapsedTime * factor) / factor;
        Rcpp::Rcout << "  CoxMGM Path Elapsed Time =  " << elapsedTime << " s" << std::endl;
    }
    return pathGraphs;
}


std::vector<EdgeListGraph> CoxMGM::searchPath(std::vector<double> lambdas) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<EdgeListGraph> pathGraphs;
    std::sort(lambdas.begin(), lambdas.end(), std::greater<double>());
    for (int i = 0; i < lambdas.size(); i++) {
	// RcppThread::Rcout << "Learning CoxMGM for lambda = " << lambdas[i] << std::endl;
	std::vector<double> lambda = { lambdas[i], lambdas[i], lambdas[i],
				       lambdas[i], lambdas[i] };
	setLambda(lambda);
	learn(1e-3, 500);
	pathGraphs.push_back(graphFromCoxMGM());
	// pathGraphs[i].setHyperParam("lambda", Rcpp::NumericVector(lambda.begin(), lambda.end()));
	// pathGraphs[i].setHyperParam("lambda", Rcpp::NumericVector(lambda.begin(), lambda.end()));
    }
    double elapsedTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-start).count();
    if (verbose) {
	double factor = (elapsedTime < 10) ? std::pow(10, 2 - std::ceil(std::log10(std::abs(elapsedTime)))) : 1.0;
	elapsedTime = std::round(elapsedTime * factor) / factor;
        Rcpp::Rcout << "  CoxMGM Path Elapsed Time =  " << elapsedTime << " s" << std::endl;
    }
    return pathGraphs;
}
