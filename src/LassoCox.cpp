#include "LassoCox.hpp"

// #include "ContinuousVariable.hpp"
// #include "DiscreteVariable.hpp"
#include <math.h>
// #include <RcppArmadillo.h>
#include <chrono>

LassoCox::LassoCox(arma::mat& x, arma::mat& y, std::vector<Variable*> variables,
		   std::vector<int>& l, std::vector<double>& lambda) {
    
    if (l.size() != y.n_cols)
        throw std::invalid_argument("length of l doesn't match number of variables in Y");

    if (y.n_rows != x.n_rows)
        throw std::invalid_argument("different number of samples for x and y");

    //lambda should have 3 values corresponding to sc, sd
    if (lambda.size() != 2)
        throw std::invalid_argument("Lambda should have two values for sc and sd edges respectively");
    
    this->xDat = x;
    this->yDat = y;

    // for (int i; i < l.size(); i++) l[i]--;
    
    this->l = l;
    this->p = x.n_cols;
    this->q = y.n_cols;
    this->n = x.n_rows;
    this->variables = variables;
    this->initVariables = variables;

    for (Variable* var : variables)
	if (var->isCensored())
	    this->censVars.push_back((CensoredVariable*) var);

    this->r = censVars.size();

    this->order = arma::umat(n, r, arma::fill::zeros);
    this->H = arma::umat(n, r, arma::fill::zeros);
    this->censor = arma::umat(n, r, arma::fill::zeros);

    for (int k = 0; k < r; k++) {
	this->order.col(k) = censVars.at(k)->getOrder();
	this->H.submat(0, k, censVars.at(k)->getH().n_elem-1, k) = censVars.at(k)->getH();
	this->censor.col(k) = censVars.at(k)->getCensor();
    }

    this->lambda = arma::vec(lambda);
    // fixData();
    initParameters();
    makeDummy();
    calcWeights();

}

LassoCox::LassoCox(DataSet& ds, std::vector<double>& lambda) {

    this->xDat = ds.getContinuousData();
    this->yDat = ds.getDiscreteData();
    this->l = ds.getDiscLevels();
    this->p = xDat.n_cols;
    this->q = yDat.n_cols;
    this->n = xDat.n_rows;

    // for (int i; i < l.size(); i++) l[i]--;

    //the variables are now ordered continuous first then discrete
    std::vector<Variable*> cVar = ds.getContinuousVariables();
    std::vector<Variable*> dVar = ds.getDiscreteVariables();
    std::vector<Variable*> sVar = ds.getCensoredVariables();

    this->r = sVar.size();
    
    this->variables = std::vector<Variable*>();
    this->variables.reserve(p+q+r);
    this->variables.insert(this->variables.end(), cVar.begin(), cVar.end());
    this->variables.insert(this->variables.end(), dVar.begin(), dVar.end());
    this->variables.insert(this->variables.end(), sVar.begin(), sVar.end());

    for (Variable* var : sVar)
	this->censVars.push_back((CensoredVariable*) var);

    this->r = sVar.size();
    
    this->initVariables = ds.getVariables();
    this->lambda = arma::vec(lambda);

    this->order = arma::umat(n, r, arma::fill::zeros);
    this->H = arma::umat(n, r, arma::fill::zeros);
    this->censor = arma::umat(n, r, arma::fill::zeros);

    for (int k = 0; k < r; k++) {
	this->order.col(k) = censVars.at(k)->getOrder();
	this->H.submat(0, k, censVars.at(k)->getH().n_elem-1, k) = censVars.at(k)->getH();
	this->censor.col(k) = censVars.at(k)->getCensor();
    }

    //Data is checked for 0 or 1 indexing and for missing levels and N(0,1) Standardizes continuous data
    fixData();

    //Initialize all parameters to zeros
    initParameters();

    //Creates dummy variables for each category of discrete variables (stored in dDat)
    makeDummy();
    
    //Sets continuous variable weights to standard deviation and discrete variable weights to p*(1-p) for each category
    calcWeights();

    // Rcpp::Rcout << "constructor complete" << std::endl;

    // //Initialize all parameters to zeros
    // initParameters();
}

// Variables are deleted in the DataSet destructor
LassoCox::~LassoCox() {
    // for (int i = 0; i < variables.size(); i++) {
    //     delete variables[i];
    // }  
}

// init all parameters to zeros except for betad which is set to 1s
void LassoCox::initParameters() {
    lcumsum = std::vector<int>(l.size()+1);
    lcumsum[0] = 0;
    // Rcpp::Rcout << "l: "; 
    for(int i = 0; i < l.size(); i++){
	// Rcpp::Rcout << l[i] << " ";
        lcumsum[i+1] = lcumsum[i] + l[i];
    }
    // Rcpp::Rcout << std::endl;
    lsum = lcumsum[l.size()];

    params = LassoCoxParams(
        arma::mat((int) xDat.n_cols, r,                 arma::fill::zeros),  // gamma
        arma::mat(lsum,              r,                 arma::fill::zeros)   // eta
    );
}

// avoid underflow in log(sum(exp(x))) calculation
double LassoCox::logsumexp(const arma::vec& x) {
    arma::vec myX = arma::vec(x);
    double maxX = myX.max();
    myX -= maxX;
    return std::log(arma::sum(arma::exp(myX))) + maxX;
}

//calculate parameter weights as in Lee and Hastie
void LassoCox::calcWeights() {
    // Rcpp::Rcout << "p = " << p << "\nq = " << q << "\nlsum = " << lsum << "\nr = " << r << std::endl;
    // Rcpp::Rcout << "xDat dims: " << xDat.n_rows << " x " << xDat.n_cols << std::endl;
    // Rcpp::Rcout << "dDat dims: " << dDat.n_rows << " x " << dDat.n_cols << std::endl;
    weights = arma::mat(p + q, r, arma::fill::ones);
    arma::vec tempWeights(p + q, arma::fill::zeros);
    arma::uvec nEvents = arma::sum(censor, 0);
    // Rcpp::Rcout << nEvents << std::endl;

    // tempWeights = arma::vec(p+q, arma::fill::ones);

    //Continuous variable weights are standard deviations
    for (arma::uword i = 0; i < p; i++) {
        tempWeights(i) = arma::stddev(xDat.col(i));
    }

    //Discrete variable weights for each variable-category pair are p(1-p) where p is the percentage of times that category appears
    for (arma::uword j = 0; j < q; j++) {
        double curWeight = 0;
        for (int k = 0; k < l[j]; k++) {
            arma::vec equalityVec = arma::vec(yDat.col(j)).transform( [k](double val) { return val == k+1 ? 1 : 0; } );
            double curp = arma::sum(equalityVec) / (double) n;
            curWeight += curp * (1-curp);
        }
        tempWeights(p+j) = std::sqrt(curWeight);
    }

    // for (int ii = 0; ii < r; ii++) {
    // 	double rs_sum = n, sub;
    // 	arma::vec numX = arma::sum(xDat, 0).t();
    // 	arma::vec numY = arma::sum(dDat, 0).t();
    // 	arma::vec sub_numX(p), sub_numY(lsum);

    // 	// Rcpp::Rcout << "numX dims: " << numX.n_elem << " x 1\n";
    // 	// Rcpp::Rcout << "numY dims: " << numY.n_elem << " x 1\n";

    // 	int i = 0;
    // 	for (int j = 0; i < n; j++) {
    // 	    sub = 0;
    // 	    sub_numX.fill(0);
    // 	    sub_numY.fill(0);

    // 	    for (int k = 0; k < H(j,ii); k++) {
    // 		if (censor(order(i+k,ii),ii)) {
    // 		    tempWeights.submat(0, ii, p-1, ii) += arma::square(xDat.row(order(i+k,ii)).t()
    // 							      - numX / rs_sum);
    // 		    arma::vec tempP = numY / rs_sum;
    // 		    tempWeights.submat(p, ii, p+lsum-1, ii) += tempP % (1 - tempP);
    // 		}
    // 		sub++;
    // 		sub_numX += xDat.row(order(i+k,ii)).t();
    // 		sub_numY += dDat.row(order(i+k,ii)).t();
    // 	    }

    // 	    i += H(j,ii);
    // 	    rs_sum -= sub;
    // 	    numX -= sub_numX;
    // 	    numY -= sub_numY;
    // 	}
    // }
    
    // tempWeights /= (n-1);

    // Rcpp::Rcout << "tempWeights:\n" << tempWeights.t() << std::endl;

    // for (arma::uword ii = 0; ii < r; ii++) {
    // 	weights.col(ii) *= (nEvents(ii) / (double) n);
    // 	// for (arma::uword j = 0; j < q; j++) {
    // 	//     weights(p+j, ii) *= std::sqrt(l[j]);
    // 	// }
    // }

    // Continuous variable weights are the square root of the sum of square deviations risk set means
    for (arma::uword ii = 0; ii < r; ii++) {
    	// tempWeights.submat(p, ii, p+lsum-1, ii) /= nEvents(ii);
    	for (arma::uword i = 0; i < p; i++) {
    	    // weights(i, ii) = std::sqrt(tempWeights(i, ii));
    	    // weights(i, ii) = 1;
	    weights(i, ii) = tempWeights(i) * (nEvents(ii) / (double) n);
    	}
    	for (arma::uword j = 0; j < q; j++) {
    	    // weights(p+j, ii) = 1 / std::sqrt(arma::accu(tempWeights.submat(p+lcumsum[j], ii,
	    // 								   p+lcumsum[j+1]-1, ii)));
    	    // weights(p+j, ii) = std::sqrt(l[j]);
	    weights(p+j, ii) = tempWeights(p+j) * (nEvents(ii) / (double) n);
    	}
    }

    
    // std::vector<uint> idx2drop;
    // for (arma::uword j = 0; j < q; j++) {
    // 	uint maxCatW = 0;
    // 	uint maxCatCol = 0;
    // 	for (int k = 0; k < l[j]; k++) {
    // 	    if (tempWeights(p+lcumsum[j]+k) > maxCatW) {
    // 		maxCatW = tempWeights(p+lcumsum[j]+k);
    // 		maxCatCol = k;
    // 	    }
    // 	}
    // 	idx2drop.push_back(lcumsum[j]+maxCatCol);
    // 	Rcpp::Rcout << "reference for col " << j << " is cat " << maxCatCol << std::endl;
    // 	l[j]--;
    // }
    // arma::uvec dropped(idx2drop);
    // Rcpp::Rcout << "Dropped indices: " << dropped.t() << std::endl;
    
    // dDat.shed_cols(dropped);
    // tempWeights.shed_rows(dropped+p);

    // // lsum -= q;
    // // for (arma::uword j = 0; j < q; j++) {
    // // 	l[j]--;
    // // 	lcumsum[j+1] -= j;
    // // }

    // lcumsum[0] = 0;
    // Rcpp::Rcout << "l: "; 
    // for(int i = 0; i < l.size(); i++){
    // 	Rcpp::Rcout << l[i] << " ";
    //     lcumsum[i+1] = lcumsum[i] + l[i];
    // }
    // Rcpp::Rcout << std::endl;
    // lsum = lcumsum[l.size()];
    
    
    // for (arma::uword j = 0; j < q; j++) {
    //     weights(p+j) = std::sqrt(arma::accu(tempWeights.subvec(p+lcumsum[j], p+lcumsum[j+1]-1)));
    // 	// weights(p+j) = std::sqrt(l[j]);
    // }

    // Rcpp::Rcout << "LassoCox weights:" << weights.t() << std::endl;

    // weights /= arma::mean(weights);
    
    // Rcpp::Rcout << "Rescaled LassoCox weights:" << weights.t() << std::endl;

    // //Continuous variable weights are square root of the sum of square deviations risk set means
    // for (arma::uword i = 0; i < p; i++) {
    //     weights(i) = arma::stddev(xDat.col(i));
    // }

    // //Discrete variable weights for each variable-category pair are p(1-p) where p is the percentage of times that category appears
    // for (arma::uword j = 0; j < q; j++) {
    //     double curWeight = 0;
    //     for (int k = 0; k < l[j]; k++) {
    //         arma::vec equalityVec = arma::vec(yDat.col(j)).transform( [k](double val) { return val == k+1 ? 1 : 0; } );
    //         double curp = arma::sum(equalityVec) / (double) n;
    //         curWeight += curp * (1-curp);
    //     }
    //     weights(p+j) = std::sqrt(curWeight);
    // }

    // Rcpp::Rcout << "Old weights:" << weights.t() << std::endl;
}

void LassoCox::calcWeights(arma::vec& parIn) {
    weights = arma::mat(p + q, r, arma::fill::ones);
    arma::vec tempWeights(p + lsum, r, arma::fill::zeros);

    LassoCoxParams par(parIn, p, lsum, r);
    
    arma::mat X = arma::join_rows(xDat, dDat);

    double HsumTheta, m, sub, d, phi;
    // double loss = 0.0;

    arma::vec sub_num(p + lsum); // , HsumThetaVec(p + lsum), Z(p + lsum);
    // arma::vec sub_diag(p + lsum), HsumDiag(p + lsum), temp(p + lsum);

    for (int ii = 0; ii < r; ii++) {
	arma::vec beta(arma::join_cols(par.gamma.col(ii), par.eta.col(ii)));

	// Rcpp::Rcout << "beta size: " << beta.n_elem << std::endl;

	// arma::vec grad(beta.n_elem, arma::fill::zeros);
	// arma::vec hessd(beta.n_elem, arma::fill::zeros);

	arma::vec theta = arma::exp(X * beta);
	double rs_sum = arma::accu(theta);

	
	arma::vec num = arma::sum(arma::diagmat(theta) * X, 0).t();
	// arma::vec diag_num = arma::sum(arma::diagmat(theta) * arma::square(X), 0).t();

	int i = 0;
	for (int j = 0; i < n; j++) {
	    // HsumTheta = 0;
	    // m = 0;
	    sub = 0;

	    // HsumThetaVec.fill(0);
	    sub_num.fill(0);

	    // HsumDiag.fill(0);
	    // sub_diag.fill(0);

	    for (int k = 0; k < H(j,ii); k++) {
		// temp = theta[order(i+k,ii)] * arma::square(X.row(order(i+k,ii)).t());
	    
		if (censor(order(i+k,ii),ii)) {
		    // m ++;
		
		    // HsumTheta += theta[order(i+k,ii)];
		    // grad += X.row(order(i+k,ii)).t();
		    // HsumThetaVec += theta[order(i+k,ii)] * X.row(order(i+k,ii)).t();
		    tempWeights.col(ii) += arma::square(X.row(order(i+k,ii)).t() - num / rs_sum);
		    // HsumDiag += temp;
		}

		sub += theta[order(i+k,ii)];
		sub_num += theta[order(i+k,ii)] * X.row(order(i+k,ii)).t();
		// sub_diag += temp;
	    }

	    // if (HsumTheta - sub > 1e-5)
	    // 	throw std::runtime_error("Error in Lasso Cox regression smoothGradient, HsumTheta > sub: " +
	    // 				 std::to_string(HsumTheta) + " > " + std::to_string(sub));
	    if (sub - rs_sum > 1e-5) {
		if ((H(j,ii) + i) != n) {
		    rs_sum = arma::accu(theta(order.submat(i, ii, n-1, ii)));
		    num = arma::sum(arma::diagmat(theta(order.submat(i, ii, n-1, ii)))
				    * X.rows(order.submat(i, ii, n-1, ii)), 0).t();
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Lasso Cox regression calcWeights, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }

	    // for (int l = 0; l < m; l++) {
	    // 	d = l / ((double) m);
	    // 	Z = num - d * HsumThetaVec;
	    // 	phi = rs_sum - d * HsumTheta;
	    // 	grad -= Z / phi;
	    // 	// hessd -= ((diag_num - d * HsumDiag) / phi - arma::square(Z) / (phi*phi));
	    // }

	    i += H(j,ii);
	    rs_sum -= sub;
	    num -= sub_num;
	    // diag_num -= sub_diag;
	}
    }

    tempWeights /= (n-1);

    // Continuous variable weights are the square root of the sum of square deviations risk set means
    for (arma::uword ii = 0; ii < r; ii++) {
    	// tempWeights.submat(p, ii, p+lsum-1, ii) /= nEvents(ii);
    	for (arma::uword i = 0; i < p; i++) {
    	    weights(i, ii) = std::sqrt(tempWeights(i, ii));
    	    // weights(i, ii) = 1;
    	}
    	for (arma::uword j = 0; j < q; j++) {
    	    weights(p+j, ii) = std::sqrt(arma::accu(tempWeights.submat(p+lcumsum[j], ii,
    								       p+lcumsum[j+1]-1, ii)));
    	    // weights(p+j, ii) = std::sqrt(l[j]);
    	}
    }

    // Rcpp::Rcout << "LassoCox weights:" << weights.t() << std::endl;
}

// convert discrete data (in yDat) to a matrix of dummy variables (stored in dDat)
void LassoCox::makeDummy() {
    // std::vector<uint> idx2drop;
    dDat = arma::mat(n, lsum, arma::fill::zeros);
    for(int i = 0; i < q; i++) {
	// uint maxCatCount = 0;
	// uint maxCatCol = 0;
        for(int j = 0; j < l[i]; j++) {
            arma::vec curCol = arma::vec(yDat.col(i)).transform( [j](double val) { return val == j+1 ? 1 : 0; } );
	    uint curSum = (uint) arma::sum(curCol);
            if (curSum == 0)
                throw std::invalid_argument("Discrete data is missing a level: variable " + variables[p+i]->getName() + " level " + std::to_string(j));
            dDat.col(lcumsum[i] + j) = curCol;

	    // if (curSum > maxCatCount) {
	    // 	maxCatCount = curSum;
	    // 	maxCatCol = j;
	    // }
        }
	// idx2drop.push_back(lcumsum[i] + i + maxCatCol);
	// Rcpp::Rcout << "reference for col " << i << " is cat " << maxCatCol << std::endl;
    }
    // Rcpp::Rcout << "Number of columns to remove: " << idx2drop.size() << std::endl;
    // Rcpp::Rcout << "Dummy made" << std::endl;
    // dDat.shed_cols(arma::uvec(idx2drop));
}

// checks if yDat is zero indexed and converts to 1 index. zscores x
void LassoCox::fixData() {
    double ymin = yDat.min();
    if(ymin < 0 || ymin > 1)
        throw std::invalid_argument("Discrete data must be either zero or one indexed. Found min index: " + std::to_string(ymin));
    
    if (ymin == 0) {
        yDat += 1;
    }

    // z-score of columns of x
    xDat.each_col( [](arma::vec& c) {c = (c - arma::mean(c)) / arma::stddev(c); } );
}

/**
 * Calculate value of g(X) and gradient of g(X) at the same time for efficiency reasons.
 *
 * @param X input Vector
 * @param Xout gradient of g(X)
 * @return value of g(X)
 */
double LassoCox::smooth(arma::vec& parIn, arma::vec& gradOutVec) {
    // LassoCoxParams par(parIn, p, lsum);

    // Rcpp::Rcout << "par: \n" << par << std::endl;

    // Rcpp::Rcout << "smooth called" << std::endl;

    LassoCoxParams gradOut(
        arma::mat(p,     r,        arma::fill::zeros),  // gamma
        arma::mat(lsum,  r,        arma::fill::zeros)   // eta
    );

    arma::mat X = arma::join_rows(xDat, dDat);

    LassoCoxParams par(parIn, p, lsum, r);

    double HsumTheta, m, sub, d, phi;
    double loss = 0.0;

    arma::vec sub_num(p + lsum), HsumThetaVec(p + lsum), Z(p + lsum);
    // arma::vec sub_diag(p + lsum), HsumDiag(p + lsum), temp(p + lsum);

    for (int ii = 0; ii < r; ii++) {
	arma::vec beta(arma::join_cols(par.gamma.col(ii), par.eta.col(ii)));

	// Rcpp::Rcout << "beta size: " << beta.n_elem << std::endl;

	arma::vec logtheta = X * beta;

	arma::vec grad(beta.n_elem, arma::fill::zeros);
	// arma::vec hessd(beta.n_elem, arma::fill::zeros);

	arma::vec theta = arma::exp(logtheta);
	double rs_sum = arma::accu(theta);


	arma::vec num = arma::sum(arma::diagmat(theta) * X, 0).t();
	// arma::vec diag_num = arma::sum(arma::diagmat(theta) * arma::square(X), 0).t();

	int i = 0;
	for (int j = 0; i < n; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;

	    HsumThetaVec.fill(0);
	    sub_num.fill(0);

	    // HsumDiag.fill(0);
	    // sub_diag.fill(0);

	    for (int k = 0; k < H(j,ii); k++) {
		// temp = theta[order(i+k,ii)] * arma::square(X.row(order(i+k,ii)).t());
	    
		if (censor(order(i+k,ii),ii)) {
		    m ++;

		    loss += logtheta[order(i+k,ii)];
		
		    HsumTheta += theta[order(i+k,ii)];
		    grad += X.row(order(i+k,ii)).t();
		    HsumThetaVec += theta[order(i+k,ii)] * X.row(order(i+k,ii)).t();
		    // HsumDiag += temp;
		}

		sub += theta[order(i+k,ii)];
		sub_num += theta[order(i+k,ii)] * X.row(order(i+k,ii)).t();
		// sub_diag += temp;
	    }

	    if (HsumTheta - sub > 1e-5)
		throw std::runtime_error("Error in Lasso Cox regression smooth, HsumTheta > sub: " +
					 std::to_string(HsumTheta) + " > " + std::to_string(sub));
	    if (sub - rs_sum > 1e-5) {
		if ((H(j,ii) + i) != n) {
		    rs_sum = arma::accu(theta(order.submat(i, ii, n-1, ii)));
		    num = arma::sum(arma::diagmat(theta(order.submat(i, ii, n-1, ii)))
				    * X.rows(order.submat(i, ii, n-1, ii)), 0).t();
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Lasso Cox regression smooth, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }
	    
	    for (int l = 0; l < m; l++) {
		d = l / ((double) m);
		Z = num - d * HsumThetaVec;
		phi = rs_sum - d * HsumTheta;
		grad -= Z / phi;
		// hessd -= ((diag_num - d * HsumDiag) / phi - arma::square(Z) / (phi*phi));

		loss -= std::log(phi);
	    }

	    i += H(j,ii);
	    rs_sum -= sub;
	    num -= sub_num;
	    // diag_num -= sub_diag;
	}
	gradOut.gamma.col(ii) = -grad.subvec(0, p-1) / (double) n;
	gradOut.eta.col(ii) = -grad.subvec(p, p+lsum-1) / (double) n;

	// gradOut.gamma.col(ii) = grad.subvec(0, p-1) / hessd.subvec(0, p-1);
	// gradOut.eta.col(ii) = grad.subvec(p, p+lsum-1) / hessd.subvec(p, p+lsum-1);
    }

    // gradOutVec = grad / hessd;
    gradOutVec = gradOut.toMatrix1D(); // -grad / (double) n;
    
    return -loss / ((double) n);
}

double LassoCox::smoothValue(arma::vec& parIn) {

    // Rcpp::Rcout << "smoothValue called" << std::endl;
    
    LassoCoxParams par(parIn, p, lsum, r);

    arma::mat X = arma::join_rows(xDat, dDat);
    double HsumTheta, m, sub;
    double loss = 0.0;

    for (int ii = 0; ii < r; ii++) {
	arma::vec beta(arma::join_cols(par.gamma.col(ii), par.eta.col(ii)));

	arma::vec logtheta = X * beta;
  
	arma::vec theta = arma::exp(logtheta);
	double rs_sum = arma::accu(theta);

	int i = 0;
	for (int j = 0; i < n; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;
	    for (int k = 0; k < H(j,ii); k++) {
		if (censor(order(i+k,ii),ii)) {
		    m++;
		    loss += logtheta[order(i+k,ii)];
		    HsumTheta += theta[order(i+k,ii)];
		}
		sub += theta[order(i+k,ii)];
	    }

	    if (HsumTheta - sub > 1e-5)
		throw std::runtime_error("Error in Lasso Cox regression smoothValue, HsumTheta > sub: " +
					 std::to_string(HsumTheta) + " > " + std::to_string(sub));
	    
	    if (sub - rs_sum > 1e-5) {
		if ((H(j,ii) + i) != n) {
		    rs_sum = arma::accu(theta(order.submat(i, ii, n-1, ii)));
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Lasso Cox regression smoothValue, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }

	    for (int l = 0; l < m; l++) {
		loss -= std::log(rs_sum - ((double) l) / (m * HsumTheta));
	    }

	    i += H(j,ii);
	    rs_sum -= sub;
	}
    }
    
    return -loss / (double) n;
}

/**
 * Calculate value of h(X) and proxOperator of h(X) at the same time for efficiency reasons.
 *
 * @param t positive parameter for prox operator
 * @param X input vector
 * @param Xout vector solution to prox_t(X)
 * @return value of h(X)
 */
double LassoCox::nonSmooth(double t, arma::vec& X, arma::vec& pX) {

    // Rcpp::Rcout << "nonSmooth called" << std::endl;
    
    double nonSmooth = 0;

    if (t <= 0)
        throw std::invalid_argument("t must be positive: " + std::to_string(t));
    
    arma::vec tlam = lambda * t;

    //Constructor copies and checks dimension
    //par is a copy so we can update it
    LassoCoxParams par(X, p, lsum, r);

    arma::mat gammaWeight = weights.submat(0, 0, arma::size(p, r));
    arma::mat gammaScale = gammaWeight * -tlam(0);
    
    gammaScale /= arma::abs(par.gamma);
    gammaScale += 1;
    gammaScale.transform( [](double val) {return std::max(val, 0.0); } );

    par.gamma = par.gamma % gammaScale;

    double gammaNorms = arma::accu(gammaWeight % arma::abs(par.gamma));

    double etaNorms = 0;
    for (arma::uword ii = 0; ii < r; ii++) {
	for (arma::uword j = 0; j < q; j++) {
	    arma::vec tempVec = par.eta.submat(lcumsum[j], ii, lcumsum[j+1]-1, ii);
	    double etaScale = std::max(0.0, 1 - tlam(1) * weights(p+j, ii) / arma::norm(tempVec, 2));
	    par.eta.submat(lcumsum[j], ii, lcumsum[j+1]-1, ii) = tempVec * etaScale;
	    etaNorms += weights(p+j, ii) * arma::norm(tempVec * etaScale, 2);
	}
    }

    pX = par.toMatrix1D();

    // calcWeights(pX);
    
    return lambda(0)*gammaNorms + lambda(1)*etaNorms;
}

/**
 * Calculate value of h(X)
 *
 * @param X input vector
 * @return value of h(X)
 */
double LassoCox::nonSmoothValue(arma::vec& parIn) {
    // Rcpp::Rcout << "nonSmoothValue called" << std::endl;
    //Dimension checked in constructor
    //par is a copy so we can update it
    LassoCoxParams par(parIn, p, lsum, r);

    double gammaNorms = arma::accu(weights.submat(0, 0, arma::size(p, r)) % arma::abs(par.gamma));

    double etaNorms = 0;
    for (arma::uword ii = 0; ii < r; ii++) {
	for (arma::uword j = 0; j < q; j++) {
	    arma::vec tempVec = par.eta.submat(lcumsum[j], ii, lcumsum[j+1]-1, ii);
	    etaNorms += weights(p+j, ii) * arma::norm(tempVec, 2);
	}
    }

    return lambda(0)*gammaNorms + lambda(1)*etaNorms;
}

/**
 * Gradient of the pseudolikelihood
 *
 * @param parIn
 * @return
 */
arma::vec LassoCox::smoothGradient(arma::vec& parIn) {
    // int n = xDat.n_rows;
    // LassoCoxParams grad;

    // Rcpp::Rcout << "smoothGradient called" << std::endl;

    LassoCoxParams gradOut(
        arma::mat(p,    r,  arma::fill::zeros),  // gamma
        arma::mat(lsum, r,  arma::fill::zeros)   // eta
    );

    arma::mat X = arma::join_rows(xDat, dDat);

    LassoCoxParams par(parIn, p, lsum, r);

    double HsumTheta, m, sub, d, phi;
    // double loss = 0.0;

    arma::vec sub_num(p + lsum), HsumThetaVec(p + lsum), Z(p + lsum);
    // arma::vec sub_diag(p + lsum), HsumDiag(p + lsum), temp(p + lsum);

    for (int ii = 0; ii < r; ii++) {
	arma::vec beta(arma::join_cols(par.gamma.col(ii), par.eta.col(ii)));

	// Rcpp::Rcout << "beta size: " << beta.n_elem << std::endl;

	arma::vec grad(beta.n_elem, arma::fill::zeros);
	// arma::vec hessd(beta.n_elem, arma::fill::zeros);

	arma::vec theta = arma::exp(X * beta);
	double rs_sum = arma::accu(theta);

	
	arma::vec num = arma::sum(arma::diagmat(theta) * X, 0).t();
	// arma::vec diag_num = arma::sum(arma::diagmat(theta) * arma::square(X), 0).t();

	int i = 0;
	for (int j = 0; i < n; j++) {
	    HsumTheta = 0;
	    m = 0;
	    sub = 0;

	    HsumThetaVec.fill(0);
	    sub_num.fill(0);

	    // HsumDiag.fill(0);
	    // sub_diag.fill(0);

	    for (int k = 0; k < H(j,ii); k++) {
		// temp = theta[order(i+k,ii)] * arma::square(X.row(order(i+k,ii)).t());
	    
		if (censor(order(i+k,ii),ii)) {
		    m ++;
		
		    HsumTheta += theta[order(i+k,ii)];
		    grad += X.row(order(i+k,ii)).t();
		    HsumThetaVec += theta[order(i+k,ii)] * X.row(order(i+k,ii)).t();
		    // HsumDiag += temp;
		}

		sub += theta[order(i+k,ii)];
		sub_num += theta[order(i+k,ii)] * X.row(order(i+k,ii)).t();
		// sub_diag += temp;
	    }

	    if (HsumTheta - sub > 1e-5)
		throw std::runtime_error("Error in Lasso Cox regression smoothGradient, HsumTheta > sub: " +
					 std::to_string(HsumTheta) + " > " + std::to_string(sub));
	    
	    if (sub - rs_sum > 1e-5) {
		if ((H(j,ii) + i) != n) {
		    rs_sum = arma::accu(theta(order.submat(i, ii, n-1, ii)));
		    num = arma::sum(arma::diagmat(theta(order.submat(i, ii, n-1, ii)))
				    * X.rows(order.submat(i, ii, n-1, ii)), 0).t();
		    if (sub - rs_sum > 1e-5) {
			throw std::runtime_error("Error in Lasso Cox regression smoothGradient, sub > rs_sum: " + std::to_string(sub) + " > " + std::to_string(rs_sum));
		    }
		} else {
		    sub = rs_sum;
		}
	    }

	    for (int l = 0; l < m; l++) {
		d = l / ((double) m);
		Z = num - d * HsumThetaVec;
		phi = rs_sum - d * HsumTheta;
		grad -= Z / phi;
		// hessd -= ((diag_num - d * HsumDiag) / phi - arma::square(Z) / (phi*phi));
	    }

	    i += H(j,ii);
	    rs_sum -= sub;
	    num -= sub_num;
	    // diag_num -= sub_diag;
	}
	gradOut.gamma.col(ii) = -grad.subvec(0, p-1) / (double) n;
	gradOut.eta.col(ii) = -grad.subvec(p, p+lsum-1) / (double) n;

	// gradOut.gamma.col(ii) = grad.subvec(0, p-1) / hessd.subvec(0, p-1);
	// gradOut.eta.col(ii) = grad.subvec(p, p+lsum-1) / hessd.subvec(p, p+lsum-1);
    }

    // gradOutVec = grad / hessd;
    return gradOut.toMatrix1D(); // -grad / (double) n;
}

/**
 * A proximal operator is the solution to this optimization problem:
 *     prox_t(x) = argmin_z \frac{1}{2t} \|x-z\|^2_2 + h(x)
 *
 * @param t positive parameter for prox operator
 * @param X input vector
 * @return vector solution to prox_t(X)
 */
arma::vec LassoCox::proximalOperator(double t, arma::vec& X) {
    if (t <= 0)
        throw std::invalid_argument("t must be positive: " + std::to_string(t));
    
    arma::vec tlam = lambda * t;

    //Constructor copies and checks dimension
    //par is a copy so we can update it
    LassoCoxParams par(X, p, lsum, r);

    arma::mat gammaWeight = weights.submat(0, 0, arma::size(p, r));
    arma::mat gammaScale = gammaWeight * -tlam(0);
    
    gammaScale /= arma::abs(par.gamma);
    gammaScale += 1;
    gammaScale.transform( [](double val) {return std::max(val, 0.0); } );

    par.gamma = par.gamma % gammaScale;

    for (arma::uword ii = 0; ii < r; ii++) {
	for (arma::uword j = 0; j < q; j++) {
	    arma::vec tempVec = par.eta.submat(lcumsum[j], ii, lcumsum[j+1]-1, ii);
	    double etaScale = std::max(0.0, 1 - tlam(1) * weights(p+j, ii) / arma::norm(tempVec, 2));
	    par.eta.submat(lcumsum[j], ii, lcumsum[j+1]-1, ii) = tempVec * etaScale;
	}
    }

    return par.toMatrix1D();
}

/**
 *  Learn LassoCox traditional way with objective function tolerance. Recommended for inference applications that need
 *  accurate pseudolikelihood
 *
 * @param epsilon tolerance in change of objective function
 * @param iterLimit iteration limit
 */
void LassoCox::learn(double epsilon, int iterLimit) {
    ProximalGradient pg = ProximalGradient();
    // ProximalGradient pg(0.5, 0.9, true);
    arma::vec curParams = params.toMatrix1D();
    arma::vec newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, epsilon, iterLimit);
    params = LassoCoxParams(newParams, p, lsum, r);
    // Rcpp::Rcout << "Iter count: " << pg.iterComplete << std::endl;
    // Rcpp::Rcout << "Time per iter: " << pg.timePerIter << " ms\n\n";
}

/**
 *  Learn LassoCox using edge convergence using default 3 iterations of no edge changes. Recommended when we only care about
 *  edge existence.
 *
 * @param iterLimit
 */
void LassoCox::learnEdges(int iterLimit) {
    ProximalGradient pg(0.5, 0.8, true);
    arma::vec curParams = params.toMatrix1D();
    arma::vec newParams;
    if (timeout != -1)
        newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 1e-10, iterLimit, timeout);
    else
        newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 1e-10, iterLimit);

    params = LassoCoxParams(newParams, p, lsum, r);

    timePerIter = pg.timePerIter;
    iterCount = pg.iterComplete;
    // Rcpp::Rcout << "Iter count: " << pg.iterComplete << std::endl;
    // Rcpp::Rcout << "Time per iter: " << pg.timePerIter << " ms\n\n";
}   

/**
 *  Learn LassoCox using edge convergence using edgeChangeTol (see ProximalGradient for documentation). Recommended when we only care about
 *  edge existence.
 *
 * @param iterLimit
 * @param edgeChangeTol
 */
void LassoCox::learnEdges(int iterLimit, int edgeChangeTol){
    // Rcpp::Rcout << "learnEdges called\n";
    ProximalGradient pg(0.5, 0.8, true);
    arma::vec curParams = params.toMatrix1D();
    // Rcpp::Rcout << "toMatrix1D succeeded\n";
    pg.setEdgeChangeTol(edgeChangeTol);
    arma::vec newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 1e-10, iterLimit);
    params = LassoCoxParams(newParams, p, lsum, r);
    // Rcpp::Rcout << "Iter count: " << pg.iterComplete << std::endl;
    // Rcpp::Rcout << "Time per iter: " << pg.timePerIter << " ms\n\n";
}

// /**
//  * Converts LassoCox to matrix of doubles. uses 2-norm to combine c-d edge parameters into single value and f-norm for
//  * d-d edge parameters.
//  *
//  * @return
//  */
// arma::mat LassoCox::adjMatFromLassoCox() {
//     arma::mat outMat(p+q, p+q, arma::fill::zeros);

//     outMat(0, 0, arma::size(p, p)) = params.beta + params.beta.t();

//     for (arma::uword i = 0; i < p; i++) {
//         for (arma::uword j = 0; j < q; j++) {
//             double val = arma::norm(params.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1), 2);
//             outMat(i, p+j) = val;
//             outMat(p+j, i) = val;
//         }
//     }

//     for (arma::uword i = 0; i < q; i++) {
//         for (arma::uword j = i+1; j < q; j++) {
//             double val = arma::norm(params.phi(lcumsum[i], lcumsum[j], arma::size(l[i], l[j])), "fro");
//             outMat(p+i, p+j) = val;
//             outMat(p+j, p+i) = val;
//         }
//     }

//     //order the adjmat to be the same as the original DataSet variable ordering
//     arma::uvec varMap(p+q);
//     for(arma::uword i = 0; i < p+q; i++){
//         varMap(i) = std::distance(variables.begin(), std::find(variables.begin(), variables.end(), initVariables[i]));
//     }
//     outMat = outMat.submat(varMap, varMap);

//     return outMat;
// }

/**
 * Adds edges from LassoCox object to EdgeListGraph object with edges if edge parameters 
 * are non-zero. Loses all edge param information
 *
 * @return
 */
void LassoCox::addEdgesFromLassoCox(EdgeListGraph& graph) {
    for (arma::uword k = 0; k < r; k++) {
	for (arma::uword i = 0; i < p; i++) {
	    if (std::abs(params.gamma(i, k)) > 0) {
		if(!graph.isAdjacentTo(censVars[k], variables[i])) {
		    graph.addUndirectedEdge(censVars[k], variables[i]);
		    // Rcpp::Rcout << "Adding " << variables[i]->getName()
		    // 		<< " --- " << censVars[k]->getName() << std::endl;
		}
	    }
	}

	for (arma::uword j = 0; j < q; j++) {
	    if (arma::norm(params.eta.submat(lcumsum[j], k, lcumsum[j+1]-1, k), 2) > 0) {
		if(!graph.isAdjacentTo(censVars[k], variables[p+j])) {
		    graph.addUndirectedEdge(censVars[k], variables[p+j]);
		    // Rcpp::Rcout << "Adding " << variables[p+j]->getName()
		    // 		<< " --- " << censVars[k]->getName() << std::endl;
		}
	    }
	}
    }
}

// /**
//  * Converts LassoCox object to Graph object with edges if edge parameters are non-zero. Loses all edge param information
//  *
//  * @return
//  */
// EdgeListGraph LassoCox::graphFromLassoCox() {
//     EdgeListGraph g(variables);

//     for (arma::uword i = 0; i < p; i++) {
//         for (arma::uword j = i+1; j < p; j++) {
//             double v1 = params.beta(i,j);

//             if (std::abs(v1) > 0) {
//                 if(!g.isAdjacentTo(variables[i], variables[j])) {
//                     g.addUndirectedEdge(variables[i], variables[j]);
//                 }
//             }
//         }
//     }

//     for (arma::uword i = 0; i < p; i++) {
//         for (arma::uword j = 0; j < q; j++) {
//             double v1 = arma::accu(arma::abs(params.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1)));

//             if (v1 > 0) {
//                 if(!g.isAdjacentTo(variables[i], variables[p+j])) {
//                     g.addUndirectedEdge(variables[i], variables[p+j]);
//                 }
//             }
//         }
//     }

//     for (arma::uword i = 0; i < q; i++) {
//         for (arma::uword j = i+1; j < q; j++) {
//             double v1 = arma::accu(arma::abs(params.phi(lcumsum[i], lcumsum[j], arma::size(l[i], l[j]))));

//             if (v1 > 0) {
//                 if(!g.isAdjacentTo(variables[p+i], variables[p+j])) {
//                     g.addUndirectedEdge(variables[p+i], variables[p+j]);
//                 }
//             }
//         }
//     }

//     // Set algorithm and type
//     std::ostringstream alg;
//     alg << "LassoCox: lambda = [" 
//         << lambda(0) << ", " << lambda(1) << ", " << lambda(2) << "]";

//     g.setAlgorithm(alg.str());
//     g.setGraphType("undirected");

//     return g;
// }

// /**
//  * Simple search command for GraphSearch implementation. Uses default edge convergence, 1000 iter limit.
//  *
//  * @return
//  */
// EdgeListGraph LassoCox::search() {
//     auto start = std::chrono::high_resolution_clock::now();
//     learnEdges(500);
//     elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
//     EdgeListGraph g(variables);
//     return g;
// }

// [[Rcpp::export]]
void LassoCoxTest(const Rcpp::DataFrame &df, const int maxDiscrete, std::vector<double> lambda) {
    DataSet ds = DataSet(df, maxDiscrete);

    // Variable* target = ds.getVariable(targetName);

    // if (!target->isCensored())
    // 	throw std::runtime_error("Target variable for Lasso Cox regression must be censored");
    
    LassoCox lassoCox(ds, lambda);
    lassoCox.learnEdges(500, 5);
    Rcpp::Rcout << lassoCox.getParams() << std::endl;

    // lassoCox = LassoCox(ds, lambda);
    
    // lassoCox.learn(1e-5, 500);
    // Rcpp::Rcout << lassoCox.getParams() << std::endl;

    EdgeListGraph g(ds.getVariables());
    lassoCox.addEdgesFromLassoCox(g);
}
