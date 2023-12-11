#include "MGM.hpp"


MGM::MGM(arma::mat& x, arma::mat& y, std::vector<Node>& variables, std::vector<int>& l, std::vector<double>& lambda) : xDat(x), yDat(y) {
    
    if (l.size() != y.n_cols)
        throw std::invalid_argument("length of l doesn't match number of variables in Y");

    if (y.n_rows != x.n_rows)
        throw std::invalid_argument("different number of samples for x and y");

    //lambda should have 3 values corresponding to cc, cd, and dd
    if (lambda.size() != 3)
        throw std::invalid_argument("Lambda should have three values for cc, cd, and dd edges respectively");

    bool mixed = true;

    Node dummyVar = Node(new DiscreteVariable("dummy.gLpkx1Hs6x", 2));

    if (std::find(variables.begin(), variables.end(), dummyVar) != variables.end()) {
	this->dummyVar = dummyVar;
	mixed = false;
	qDummy = 1;
    }

    dummyVar = Node(new ContinuousVariable("dummy.qCm6jaC1VK"));

    if (std::find(variables.begin(), variables.end(), dummyVar) != variables.end()) {
	this->dummyVar = dummyVar;
	mixed = false;
	pDummy = 1;
    }
    
    // this->xDat = x;
    // this->yDat = y;
    this->l = l;
    this->p = x.n_cols;
    this->q = y.n_cols;
    this->n = x.n_rows;
    this->variables = variables;
    this->initVariables = variables;

    this->lambda = arma::vec(lambda);
    // fixData();
    initParameters();
    calcWeights();
    makeDummy();

}

MGM::MGM(arma::mat&& x, arma::mat&& y, std::vector<Node> variables, std::vector<int> l, std::vector<double> lambda) : xDat(x), yDat(y) {
    
    if (l.size() != y.n_cols)
        throw std::invalid_argument("length of l doesn't match number of variables in Y");

    if (y.n_rows != x.n_rows)
        throw std::invalid_argument("different number of samples for x and y");

    //lambda should have 3 values corresponding to cc, cd, and dd
    if (lambda.size() != 3)
        throw std::invalid_argument("Lambda should have three values for cc, cd, and dd edges respectively");

    bool mixed = true;

    Node dummyVar = Node(new DiscreteVariable("dummy.gLpkx1Hs6x", 2));

    if (std::find(variables.begin(), variables.end(), dummyVar) != variables.end()) {
	this->dummyVar = dummyVar;
	mixed = false;
	qDummy = 1;
    }

    dummyVar = Node(new ContinuousVariable("dummy.qCm6jaC1VK"));

    if (std::find(variables.begin(), variables.end(), dummyVar) != variables.end()) {
	this->dummyVar = dummyVar;
	mixed = false;
	pDummy = 1;
    }
    
    // std::swap(this->xDat, x);
    // std::swap(this->yDat, y);
    this->l = l;
    this->p = x.n_cols;
    this->q = y.n_cols;
    this->n = x.n_rows;
    this->variables = variables;
    this->initVariables = variables;

    this->lambda = arma::vec(lambda);
    // fixData();
    initParameters();
    calcWeights();
    makeDummy();

}


MGM::MGM(DataSet ds) {

    bool mixed = true;

    if (ds.isContinuous()) {
	this->dummyVar = Node(new DiscreteVariable("dummy.gLpkx1Hs6x", 2));
	ds.addVariable(dummyVar);
	arma::uword j = ds.getColumn(dummyVar);
	for (arma::uword i = 0; i < ds.getNumRows(); i++) {
	    ds.set(i, j, std::floor(R::runif(0,2)));
	}
	mixed = false;
	qDummy = 1;
    }

    if (ds.isDiscrete()) {
	this->dummyVar = Node(new ContinuousVariable("dummy.qCm6jaC1VK"));
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
    this->l = ds.getDiscLevels();
    this->p = xDat.n_cols;
    this->q = yDat.n_cols;
    this->n = xDat.n_rows;

    //the variables are now ordered continuous first then discrete
    std::vector<Node> cVar = ds.getContinuousVariables();
    std::vector<Node> dVar = ds.getDiscreteVariables();
    this->variables = std::vector<Node>();
    this->variables.reserve(p+q);
    this->variables.insert(this->variables.end(), cVar.begin(), cVar.end());
    this->variables.insert(this->variables.end(), dVar.begin(), dVar.end());
    
    this->initVariables = ds.getVariables();
    this->lambda = 5 * std::sqrt(std::log10(p+q)/((double)n)) * arma::vec(3, arma::fill::ones);

    // if (!mixed) {
    // 	variables.erase(std::find(variables.begin(), variables.end(), dummyVar));
    // 	initVariables.erase(std::find(initVariables.begin(), initVariables.end(), dummyVar));
    // }
    //Data is checked for 0 or 1 indexing and for missing levels and N(0,1) Standardizes continuous data
    fixData();

    //Initialize all parameters to zeros
    initParameters();

    //Sets continuous variable weights to standard deviation and discrete variable weights to p*(1-p) for each category
    calcWeights();

    //Creates dummy variables for each category of discrete variables (stored in dDat)
    makeDummy();
}


MGM::MGM(DataSet ds, std::vector<double>& lambda) {

    bool mixed = true;

    if (ds.isContinuous()) {
	this->dummyVar = Node(new DiscreteVariable("dummy.gLpkx1Hs6x", 2));
	ds.addVariable(dummyVar);
	arma::uword j = ds.getColumn(dummyVar);
	for (arma::uword i = 0; i < ds.getNumRows(); i++) {
	    ds.set(i, j, std::floor(R::runif(0,2)));
	}
	mixed = false;
	qDummy = 1;
    }

    if (ds.isDiscrete()) {
	this->dummyVar = Node(new ContinuousVariable("dummy.qCm6jaC1VK"));
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
    this->l = ds.getDiscLevels();
    this->p = xDat.n_cols;
    this->q = yDat.n_cols;
    this->n = xDat.n_rows;

    //the variables are now ordered continuous first then discrete
    std::vector<Node> cVar = ds.getContinuousVariables();
    std::vector<Node> dVar = ds.getDiscreteVariables();
    this->variables = std::vector<Node>();
    this->variables.reserve(p+q);
    this->variables.insert(this->variables.end(), cVar.begin(), cVar.end());
    this->variables.insert(this->variables.end(), dVar.begin(), dVar.end());
    
    this->initVariables = ds.getVariables();
    this->lambda = arma::vec(lambda);

    // if (!mixed) {
    // 	variables.erase(std::find(variables.begin(), variables.end(), dummyVar));
    // 	initVariables.erase(std::find(initVariables.begin(), initVariables.end(), dummyVar));
    // }
    //Data is checked for 0 or 1 indexing and for missing levels and N(0,1) Standardizes continuous data
    fixData();

    //Initialize all parameters to zeros
    initParameters();

    //Sets continuous variable weights to standard deviation and discrete variable weights to p*(1-p) for each category
    calcWeights();

    //Creates dummy variables for each category of discrete variables (stored in dDat)
    makeDummy();
}


// init all parameters to zeros except for betad which is set to 1s
void MGM::initParameters() {
    lcumsum = std::vector<int>(l.size()+1);
    lcumsum[0] = 0;
    for(int i = 0; i < l.size(); i++){
        lcumsum[i+1] = lcumsum[i] + l[i];
    }
    lsum = lcumsum[l.size()];

    // arma::mat beta((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros);

    params = MGMParams(
        arma::mat((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros),  // beta
        arma::vec((int) xDat.n_cols,                    arma::fill::ones),   // betad
        arma::mat(lsum,              (int) xDat.n_cols, arma::fill::zeros),  // theta
        arma::mat(lsum,              lsum,              arma::fill::zeros),  // phi
        arma::vec((int) xDat.n_cols,                    arma::fill::zeros),  // alpha1
        arma::vec(lsum,                                 arma::fill::zeros)   // alpha2
    );
}

// avoid underflow in log(sum(exp(x))) calculation
double MGM::logsumexp(const arma::vec& x) {
    arma::vec myX = arma::vec(x);
    double maxX = myX.max();
    myX -= maxX;
    return std::log(arma::sum(arma::exp(myX))) + maxX;
}

//calculate parameter weights as in Lee and Hastie
void MGM::calcWeights() {
    weights = arma::vec(p+q);

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
}

// convert discrete data (in yDat) to a matrix of dummy variables (stored in dDat)
void MGM::makeDummy() {
    dDat = arma::mat(n, lsum, arma::fill::zeros);
    for(int i = 0; i < q; i++) {
        for(int j = 0; j < l[i]; j++) {
            arma::vec curCol = arma::vec(yDat.col(i)).transform( [j](double val) { return val == j+1 ? 1 : 0; } );
            // if (arma::sum(curCol) == 0)
            //     throw std::invalid_argument("Discrete data is missing a level: variable " + variables[p+i].getName() + " level " + std::to_string(j));
            dDat.col(lcumsum[i]+j) = curCol;
        }
    }
}

// checks if yDat is zero indexed and converts to 1 index. zscores x
void MGM::fixData() {
    double ymin = yDat.min();
    if(ymin < 0 || ymin > 1)
        throw std::invalid_argument("Discrete data must be either zero or one indexed. Found min index: " + std::to_string(ymin));
    
    if (ymin == 0) {
        yDat += 1;
    }

    // z-score of columns of x
    xDat.each_col( [](arma::vec& c) {c = (c - arma::mean(c)) / arma::stddev(c); } );
}

double MGM::calcLambdaMax() {
    lcumsum = std::vector<int>(l.size()+1);
    lcumsum[0] = 0;
    for(int i = 0; i < l.size(); i++){
        lcumsum[i+1] = lcumsum[i] + l[i];
    }
    lsum = lcumsum[l.size()];

    arma::mat beta((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros);

    MGMParams nullParams(
        arma::mat((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros),  // beta
        arma::vec((int) xDat.n_cols,                    arma::fill::ones),   // betad
        arma::mat(lsum,              (int) xDat.n_cols, arma::fill::zeros),  // theta
        arma::mat(lsum,              lsum,              arma::fill::zeros),  // phi
        arma::vec((int) xDat.n_cols,                    arma::fill::zeros),  // alpha1
        arma::vec(lsum,                                 arma::fill::zeros)   // alpha2
    );

    arma::vec nullParams1D = nullParams.toMatrix1D();
    arma::vec nullGrad1D = smoothGradient(nullParams1D);

    MGMParams nullGrad(nullGrad1D, p, lsum);

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

    // Rcpp::Rcout << "lambdaMax after phi = " << lambdaMax << std::endl;

    return lambdaMax;

}

/**
 * Calculate value of g(X) and gradient of g(X) at the same time for efficiency reasons.
 *
 * @param X input Vector
 * @param Xout gradient of g(X)
 * @return value of g(X)
 */
double MGM::smooth(arma::vec& parIn, arma::vec& gradOutVec) {
    MGMParams par(parIn, p, lsum);

    // Rcpp::Rcout << "par: \n" << par << std::endl;

    MGMParams gradOut;

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

    // Squared loss
    //tempLoss =  (X-e*alpha1'-Xbeta-Dtheta) = -res (in gradient code)
    arma::mat tempLoss(n, xDat.n_cols);

    //wxprod=X*(theta')+D*phi+e*alpha2';
    arma::mat wxProd = xDat * par.theta.t() + dDat * par.phi;

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < xDat.n_cols; j++) {
	    tempLoss(i, j) = xDat(i,j) - par.alpha1(j) / par.betad(j) - xBeta(i,j) - dTheta(i,j);
        }
        for (arma::uword j = 0; j < dDat.n_cols; j++) {
            wxProd(i,j) = wxProd(i,j) + par.alpha2(j);
        }
    }

    //sqloss=-n/2*sum(log(betad))+...
    //.5*norm((X-e*alpha1'-Xbeta-Dtheta)*diag(sqrt(betad)),'fro')^2;
    double sqloss = -n/2.0*arma::sum(arma::log(arma::vec(par.betad))) +
                    0.5 * std::pow(arma::norm(tempLoss * arma::diagmat(arma::sqrt(par.betad)), "fro"), 2);

    //ok now tempLoss = res
    tempLoss *= -1;

    //gradbeta=X'*(res);
    gradOut.beta = xDat.t() * tempLoss;

    //gradbeta=gradbeta-diag(diag(gradbeta)); % zero out diag
    //gradbeta=tril(gradbeta)'+triu(gradbeta);
    gradOut.beta = arma::trimatl(gradOut.beta, -1).t() + arma::trimatu(gradOut.beta, 1);

    //gradalpha1=diag(betad)*sum(res,0)';
    gradOut.alpha1 = arma::sum(tempLoss, 0).t();

    //gradtheta=D'*(res);
    gradOut.theta = dDat.t() * tempLoss;

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
        wxTemp = arma::exp(wxTemp);
        arma::vec invDenom = arma::vec(n, arma::fill::ones) / arma::sum(wxTemp, 1);
        wxTemp = arma::diagmat(invDenom) * wxTemp;

        for (arma::uword k = 0; k < n; k++) {
            const arma::vec& curRow0 = wxTemp0.row(k);

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

    //gradphi=tril(gradphi)'+triu(gradphi);
    gradOut.phi = arma::trimatl(gradOut.phi, 0).t() + arma::trimatu(gradOut.phi, 0);

    //zero out gradphi diagonal and ensure each group in theta and phi sum to 0 
    //for r=1:q
    //gradphi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    for (arma::uword i = 0; i < q; i++) {
        gradOut.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
        // for (arma::uword j = 0; j < p; i++) {
        //     arma::subview<double> tempVec = gradOut.theta.col(j).subvec(lcumsum[i], lcumsum[i+1]-1);
        //     gradOut.theta.col(j).subvec(lcumsum[i], lcumsum[i+1]-1) -= arma::mean(tempVec);
        // }
	// for (arma::uword j = i+1; j < q; j++) {
        //     arma::subview<double> tempMat = gradOut.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
        //     gradOut.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) -= arma::mean(arma::mean(tempMat));
        // }
    }

    // for (arma::uword i = 0; i < p; i++) {
    //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
    //         const arma::vec& tempVec = gradOut.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
    //         gradOut.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) -= arma::mean(tempVec);
    //     }
    // }

    // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
    //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
    //         const arma::mat& tempMat = gradOut.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
    //         // Use the tempMat subview again to set the values (doesn't work with const)
    //         gradOut.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) -= arma::mean(arma::mean(tempMat));
    //     }
    // }


    /*
    for s=1:p
        gradbetad(s)=-n/(2*betad(s))+1/2*norm(res(:,s))^2-res(:,s)'*(Xbeta(:,s)+Dtheta(:,s));
    end
    */
    gradOut.betad = arma::vec(xDat.n_cols);
    for(arma::uword i = 0; i < p; i++){
        gradOut.betad(i) = -n / (2.0 * par.betad(i)) + std::pow(arma::norm(tempLoss.col(i), 2), 2) / 2.0 -
                           arma::as_scalar(tempLoss.col(i).t() * (xBeta.col(i) + dTheta.col(i)));
    }

    gradOut.alpha1 /= (double) n;
    gradOut.alpha2 /= (double) n;
    gradOut.betad /= (double) n;
    gradOut.beta /= (double) n;
    gradOut.theta /= (double) n;
    gradOut.phi /= (double) n;

    // Rcpp::Rcout << "gradOut: \n" << gradOut << std::endl;

    if (pDummy != 0) {
	gradOut.alpha1.fill(0);
	gradOut.betad.fill(0);
	gradOut.beta.fill(0);
	gradOut.theta.fill(0);
    }

    if (qDummy != 0) {
	gradOut.alpha2.fill(0);
	gradOut.theta.fill(0);
	gradOut.phi.fill(0);
    }

    gradOutVec = gradOut.toMatrix1D();

    return (sqloss + catloss)/((double) n);
}

double MGM::smoothValue(arma::vec& parIn) {
    MGMParams par(parIn, p, lsum);

    // Rcpp::Rcout << "par: \n" << par << std::endl;

    for(arma::uword i = 0; i < par.betad.size(); i++) {
        if(par.betad(i) <= 0)
            return std::numeric_limits<double>::infinity();
    }

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

    // Rcpp::Rcout << "xBeta: \n" << xBeta << std::endl;
    // Rcpp::Rcout << "dTheta: \n" << dTheta << std::endl;

    // Squared loss
    //sqloss=-n/2*sum(log(betad))+...
    //.5*norm((X-e*alpha1'-Xbeta-Dtheta)*diag(sqrt(betad)),'fro')^2;
    arma::mat tempLoss(n, xDat.n_cols);

    //wxprod=X*(theta')+D*phi+e*alpha2';
    arma::mat wxProd = xDat * par.theta.t() + dDat * par.phi;

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < xDat.n_cols; j++) {
            tempLoss(i, j) = xDat(i,j) - par.alpha1(j)/par.betad(j) - xBeta(i,j) - dTheta(i,j);
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
            arma::vec curRow = wxTemp.row(k);

            catloss -= curRow((arma::uword) yDat(k,i) - 1);
            catloss += logsumexp(curRow);
        }
    }

    // Rcpp::Rcout << "catloss: " << catloss << std::endl;

    return (sqloss + catloss) / (double) n;
}

/**
 * Calculate value of h(X) and proxOperator of h(X) at the same time for efficiency reasons.
 *
 * @param t positive parameter for prox operator
 * @param X input vector
 * @param Xout vector solution to prox_t(X)
 * @return value of h(X)
 */
double MGM::nonSmooth(double t, arma::vec& X, arma::vec& pX) {
    double nonSmooth = 0;

    arma::vec tlam = lambda * t;

    //Constructor copies and checks dimension
    //par is a copy so we can update it
    MGMParams par(X, p, lsum);

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

    pX = par.toMatrix1D();
    return lambda(0)*betaNorms + lambda(1)*thetaNorms + lambda(2)*phiNorms;
}

/**
 * Calculate value of h(X)
 *
 * @param X input vector
 * @return value of h(X)
 */
double MGM::nonSmoothValue(arma::vec& parIn) {
    //Dimension checked in constructor
    //par is a copy so we can update it
    MGMParams par(parIn, p, lsum);

    //penbeta = t(1).*(wv(1:p)'*wv(1:p));
    //betascale=zeros(size(beta));
    //betascale=max(0,1-penbeta./abs(beta));
    arma::mat weightMat = weights * weights.t();

    //weight beta
    //betaw = (wv(1:p)'*wv(1:p)).*abs(beta);
    //betanorms=sum(betaw(:));
    double betaNorms = arma::accu(arma::mat(weightMat.submat(0, 0, p-1, p-1)) % arma::abs(par.beta));

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
            thetaNorms += weightMat(i, p+j) * arma::norm(tempVec, 2);
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
            phiNorms += weightMat(p+i, p+j) * arma::norm(tempMat, "fro");
        }
    }

    // Rcpp::Rcout << "NSV phiNorms = " << phiNorms << std::endl;

    return lambda(0)*betaNorms + lambda(1)*thetaNorms + lambda(2)*phiNorms;
}

/**
 * Gradient of the pseudolikelihood
 *
 * @param parIn
 * @return
 */
arma::vec MGM::smoothGradient(arma::vec& parIn) {
    int n = xDat.n_rows;
    MGMParams grad;

    MGMParams par(parIn, p, lsum);

    // Rcpp::Rcout << "par: \n" << par << std::endl;

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
    // Rcpp::Rcout << "xBeta: \n" << xBeta << std::endl;
    // Rcpp::Rcout << "dTheta: \n" << dTheta << std::endl;

    //res=Xbeta-X+e*alpha1'+Dtheta;
    //wxprod=X*(theta')+D*phi+e*alpha2';
    arma::mat negLoss(n, xDat.n_cols);

    arma::mat wxProd = xDat*par.theta.t() + dDat*par.phi;
    // Rcpp::Rcout << "wxProd1: \n" << wxProd << std::endl;

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < p; j++) {
            negLoss(i,j) = xBeta(i,j) - xDat(i,j) + par.alpha1(j)/par.betad(j) + dTheta(i,j);
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
    grad.alpha1 = arma::sum(negLoss, 0).t();

    // Rcpp::Rcout << "grad.alpha1: \n" << grad.alpha1 << std::endl;

    //gradtheta=D'*(res);
    grad.theta = dDat.t() * negLoss;
    // Rcpp::Rcout << "grad.theta: \n" << grad.theta << std::endl;

    for (arma::uword i = 0; i < yDat.n_cols; i++) {
        arma::subview<double> wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));

        // does this need to be done in log space??
        wxTemp = arma::exp(wxTemp);
        arma::vec invDenom = arma::vec(n, arma::fill::ones) / arma::sum(wxTemp, 1);
        wxTemp = arma::diagmat(invDenom) * wxTemp;

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

    //zero out gradphi diagonal
    //for r=1:q
    //gradphi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    // for (arma::uword i = 0; i < q; i++) {
    //     grad.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    // }
    
    //gradphi=tril(gradphi)'+triu(gradphi);
    grad.phi = arma::trimatl(grad.phi, 0).t() + arma::trimatu(grad.phi, 0);


    //zero out gradphi diagonal and ensure each group in gradtheta and gradphi sum to 0 
    //for r=1:q
    //gradphi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    for (arma::uword i = 0; i < q; i++) {
        grad.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
        // for (arma::uword j = 0; j < p; i++) {
        //     arma::subview<double> tempVec = grad.theta.col(j).subvec(lcumsum[i], lcumsum[i+1]-1);
        //     tempVec -= arma::mean(tempVec);
        // }
	// for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
        //     arma::subview<double> tempMat = grad.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
        //     tempMat -= arma::mean(arma::mean(tempMat));
        // }
    }

    // for (arma::uword i = 0; i < p; i++) {
    //     for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
    //         const arma::vec& tempVec = grad.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1);
    //         grad.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1) -= arma::mean(tempVec);
    //     }
    // }

    // for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
    //     for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
    //         const arma::mat& tempMat = grad.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1);
    //         // Use the tempMat subview again to set the values (doesn't work with const)
    //         grad.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) -= arma::mean(arma::mean(tempMat));
    //     }
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
                           arma::as_scalar(negLoss.col(i).t() * (xBeta.col(i) + dTheta.col(i)));
    }

    grad.alpha1 /= (double) n;
    grad.alpha2 /= (double) n;
    grad.betad /= (double) n;
    grad.beta /= (double) n;
    grad.theta /= (double) n;
    grad.phi /= (double) n;

    if (pDummy != 0) {
	grad.alpha1.fill(0);
	grad.betad.fill(0);
	grad.beta.fill(0);
	grad.theta.fill(0);
    }

    if (qDummy != 0) {
	grad.alpha2.fill(0);
	grad.theta.fill(0);
	grad.phi.fill(0);
    }

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
arma::vec MGM::proximalOperator(double t, arma::vec& X) {
    if (t <= 0)
        throw std::invalid_argument("t must be positive: " + std::to_string(t));
    
    arma::vec tlam = lambda * t;

    //Constructor copies and checks dimension
    //par is a copy so we can update it
    MGMParams par(X, p, lsum);

    //penbeta = t(1).*(wv(1:p)'*wv(1:p));
    //betascale=zeros(size(beta));
    //betascale=max(0,1-penbeta./abs(beta));
    arma::mat weightMat = weights * weights.t();

    const arma::mat& betaWeight = weightMat.submat(0, 0, p-1, p-1);
    arma::mat betaScale = betaWeight * -tlam(0);
    
    betaScale /= arma::abs(par.beta);
    betaScale += 1;
    betaScale.transform( [](double val) {return std::max(val, 0.0); } );

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
            double thetaScale = std::max(0.0, 1 - tlam(1)*weightMat(i,p+j)/arma::norm(tempVec, 2));
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
            double phiScale = std::max(0.0, 1 - tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, "fro"));
            // Use the tempMat subview again to set the values (doesn't work with const)
            par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1]-1, lcumsum[j+1]-1) = tempMat * phiScale;
        }
    }

    return par.toMatrix1D();
}

/**
 *  Learn MGM traditional way with objective function tolerance. Recommended for inference applications that need
 *  accurate pseudolikelihood
 *
 * @param epsilon tolerance in change of objective function
 * @param iterLimit iteration limit
 */
void MGM::learn(double epsilon, int iterLimit) {
    ProximalGradient pg = ProximalGradient();
    arma::vec curParams = params.toMatrix1D();
    arma::vec newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, epsilon, iterLimit);
    params = MGMParams(newParams, p, lsum);
}

/**
 *  Learn MGM using edge convergence using default 3 iterations of no edge changes. Recommended when we only care about
 *  edge existence.
 *
 * @param iterLimit
 */
void MGM::learnEdges(int iterLimit) {
    ProximalGradient pg(0.5, 0.9, true);
    arma::vec curParams = params.toMatrix1D();
    arma::vec newParams;
    if (timeout != -1)
        newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 0.0, iterLimit, timeout);
    else
        newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 0.0, iterLimit);

    params = MGMParams(newParams, p, lsum);

    timePerIter = pg.timePerIter;
    iterCount = pg.iterComplete;
}   

/**
 *  Learn MGM using edge convergence using edgeChangeTol (see ProximalGradient for documentation). Recommended when we only care about
 *  edge existence.
 *
 * @param iterLimit
 * @param edgeChangeTol
 */
void MGM::learnEdges(int iterLimit, int edgeChangeTol){
    ProximalGradient pg(0.5, 0.9, true);
    arma::vec curParams = params.toMatrix1D();
    pg.setEdgeChangeTol(edgeChangeTol);
    arma::vec newParams = pg.learnBackTrack((ConvexProximal *) this, curParams, 0.0, iterLimit);
    params = MGMParams(newParams, p, lsum);
}

/**
 * Converts MGM to matrix of doubles. uses 2-norm to combine c-d edge parameters into single value and f-norm for
 * d-d edge parameters.
 *
 * @return
 */
arma::mat MGM::adjMatFromMGM() {
    arma::mat outMat(p+q-pDummy-qDummy, p+q-pDummy-qDummy, arma::fill::zeros);

    if (p - pDummy > 0) 
	outMat(0, 0, arma::size(p, p)) = params.beta + params.beta.t();

    for (arma::uword i = 0; i < p-pDummy; i++) {
        for (arma::uword j = 0; j < q-qDummy; j++) {
            double val = arma::norm(params.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]-1), 2);
            outMat(i, p+j) = val;
            outMat(p+j, i) = val;
        }
    }

    for (arma::uword i = 0; i < q-qDummy; i++) {
        for (arma::uword j = i+1; j < q-qDummy; j++) {
            double val = arma::norm(params.phi(lcumsum[i], lcumsum[j], arma::size(l[i], l[j])), "fro");
            outMat(p+i, p+j) = val;
            outMat(p+j, p+i) = val;
        }
    }

    //order the adjmat to be the same as the original DataSet variable ordering
    arma::uvec varMap(p+q-pDummy-qDummy);
    for(arma::uword i = 0; i < p+q-pDummy-qDummy; i++){
        varMap(i) = std::distance(variables.begin(), std::find(variables.begin(), variables.end(), initVariables[i]));
    }
    outMat = outMat.submat(varMap, varMap);

    return outMat;
}


/**
 * Converts MGM object to Graph object with edges if edge parameters are non-zero. Loses all edge param information
 *
 * @return
 */
EdgeListGraph MGM::graphFromMGM() {
    std::vector<Node> initVars(initVariables);
    if (!dummyVar.isNull()) {
    	initVars.erase(std::find(initVars.begin(), initVars.end(), dummyVar));
    }
    
    EdgeListGraph g(initVars);

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

    // Set algorithm and type
    std::ostringstream alg;
    alg << "MGM";// : lambda = [" 
        // << lambda(0) << ", " << lambda(1) << ", " << lambda(2) << "]";

    // RcppThread::Rcout << "lambda = [" << lambda(0) << ", " << lambda(1) << ", " << lambda(2) << "]\n";
    g.setAlgorithm(alg.str());
    g.setGraphType("undirected");
    
    // Rcpp::NumericVector rLambda(lambda.begin(), lambda.end());
    // for (double l : lambda) rLambda.push_back(l);
    
    // g.setHyperParam("lambda", Rcpp::NumericVector(lambda.begin(), lambda.end()));

    g.setHyperParam("lambda", lambda);

    std::vector<std::string> names;

    for (int i = 0; i < p; i++) {
	// Rcpp::Rcout << variables[i].getName() << std::endl;
	names.push_back(variables[i].getName());
    }

    for (int i = 0; i < q; i++) {
	// Rcpp::Rcout << variables[p+i].getName() << std::endl << "  ";
	for (std::string cat : variables[p+i].getCategories()) {
	    // Rcpp::Rcout << variables[p+i].getName() + "." + cat << " ";
	    names.push_back(variables[p+i].getName() + ":" + cat);
	}
	// Rcpp::Rcout << std::endl;
    }

    params.setNames(names);

    g.setParams(params);

    return g;
}

/**
 * Simple search command for GraphSearch implementation. Uses default edge convergence, 1000 iter limit.
 *
 * @return
 */
EdgeListGraph MGM::search() {
    auto start = std::chrono::high_resolution_clock::now();
    if (verbose) {
	RcppThread::Rcout << "  Learning MGM for lambda = { "
			  << lambda[0] << " " << lambda[1] << " "
			  << lambda[2] << " }\n";
    }
    learn(1e-5, 500);
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
    return graphFromMGM();
}


std::vector<EdgeListGraph> MGM::searchPath(std::vector<double> lambdas,
					  arma::vec& loglik,
					  arma::vec& nParams) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<EdgeListGraph> pathGraphs;
    std::sort(lambdas.begin(), lambdas.end(), std::greater<double>());
    for (int i = 0; i < lambdas.size(); i++) {
	if (verbose) RcppThread::Rcout << "  Learning MGM for lambda = " << lambdas[i] << "\n";
	std::vector<double> lambda = { lambdas[i], lambdas[i], lambdas[i] };
	setLambda(lambda);
	learn(1e-5, 500);
	pathGraphs.push_back(graphFromMGM());

	arma::vec par(params.toMatrix1D());
	loglik(i) = -n * smoothValue(par);
	nParams(i) = arma::accu(par!=0);

	RcppThread::checkUserInterrupt();
		
    }
    if (verbose) RcppThread::Rcout << std::endl;;
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
    return pathGraphs;
}


std::vector<EdgeListGraph> MGM::searchPathCV(std::vector<double> lambdas,
					     int nfolds,
					     arma::uvec& foldid,
					     arma::mat& loglik) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<EdgeListGraph> cvGraphs;
    // std::vector<MGM> trainMGMs;
    // std::vector<MGM> testMGMs;
    std::sort(lambdas.begin(), lambdas.end(), std::greater<double>());
    std::vector<double> lambda = { lambdas[0], lambdas[0], lambdas[0] };

    for (int k = 1; k <= nfolds; k++) {
	arma::uvec trainIdxs = arma::find(foldid != k);
	arma::uvec testIdxs = arma::find(foldid == k);

	std::vector<Node> variables(this->variables);
	std::vector<int> l(this->l);

	lambda = { lambdas[0], lambdas[0], lambdas[0] };
	
	MGM trainMGM(arma::mat(xDat.rows(trainIdxs)),
		     arma::mat(yDat.rows(trainIdxs)),
		     variables, l, lambda);

	trainMGM.setVerbose(false);
	
        MGM testMGM(arma::mat(xDat.rows(testIdxs)),
		    arma::mat(yDat.rows(testIdxs)),
		    variables, l, lambda);
	
	for (int i = 0; i < lambdas.size(); i++) {
	    if (verbose) {
		RcppThread::Rcout << "  Fold " << k << ": Learning MGM for lambda = "
				  << lambdas[i] << "\r";
	    }
	    
	    lambda = { lambdas[i], lambdas[i], lambdas[i] };

	    // Rcpp::Rcout << "Setting Lambda...\n";
	    
	    trainMGM.setLambda(lambda);
	    
	    // Rcpp::Rcout << "Training...\n";
	    
	    trainMGM.learn(1e-5, 500);
	    // pathGraphs.push_back(graphFromMGM());

	    arma::vec par(trainMGM.getParams().toMatrix1D());

	    // Rcpp::Rcout << "Evaluating on test data...\n";
	    loglik(i,k-1) = testMGM.smoothValue(par);
	    // nParams(i) = arma::accu(par!=0);

	    RcppThread::checkUserInterrupt();
		
	}

	if (verbose) {
	  RcppThread::Rcout << std::endl;
	}
    }

    if (verbose) RcppThread::Rcout << std::endl;

    // Rcpp::Rcout << "Test LogLiks:\n" << loglik << std::endl;

    arma::vec meanLoglik = arma::mean(loglik, 1);
    arma::vec seLoglik = arma::stddev(loglik, 0, 1);

    // Rcpp::Rcout << "Mean Test LogLiks:\n" << meanLoglik.t() << std::endl;
    // Rcpp::Rcout << "SE Test LogLiks:\n" << seLoglik.t() << std::endl;

    arma::uword minIdx = arma::index_min(meanLoglik);

    arma::uword seIdx = minIdx;
    for (arma::uword i = minIdx; i >= 0; i--) {
	if (meanLoglik(i) == meanLoglik(minIdx)) minIdx = i;
	else if (meanLoglik(i) > meanLoglik(minIdx) + seLoglik(minIdx)) break;
	seIdx = i;
    }

    lambda = { lambdas[seIdx], lambdas[seIdx], lambdas[seIdx] };

    if (verbose) RcppThread::Rcout << "lambda (1 SE) = " << lambdas[seIdx] << "\n";

    setLambda(lambda);
    learn(1e-5, 500);

    cvGraphs.push_back(graphFromMGM());

    lambda = { lambdas[minIdx], lambdas[minIdx], lambdas[minIdx] };

    if (verbose) RcppThread::Rcout << "lambda (min) = " << lambdas[minIdx] << "\n";

    setLambda(lambda);
    learn(1e-5, 500);

    cvGraphs.push_back(graphFromMGM());
    
    if (verbose) RcppThread::Rcout << std::endl;
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count();
    return cvGraphs;
}
