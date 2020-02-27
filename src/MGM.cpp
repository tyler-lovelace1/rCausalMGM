#include "MGM.hpp"

MGM::MGM(arma::mat& x, arma::mat& y, std::vector<Variable*>& variables, std::vector<int>& l, std::vector<double>& lambda) {
    
    if (l.size() != y.n_cols)
        throw std::invalid_argument("length of l doesn't match number of variables in Y");

    if (y.n_rows != x.n_rows)
        throw std::invalid_argument("different number of samples for x and y");

    //lambda should have 3 values corresponding to cc, cd, and dd
    if (lambda.size() != 3)
        throw std::invalid_argument("Lambda should have three values for cc, cd, and dd edges respectively");
    
    this->xDat = x;
    this->yDat = y;
    this->l = l;
    this->p = x.n_cols;
    this->q = y.n_cols;
    this->n = x.n_rows;
    this->variables = variables;

    this->lambda = arma::vec(lambda);
    fixData();
    initParameters();
    calcWeights();
    makeDummy();

}

MGM::MGM(DataSet& ds, std::vector<double>& lambda) {

    this->xDat = ds.getContinuousData();
    this->yDat = ds.getDiscreteData();

    //TODO 
    // this->l = ???

    this->p = xDat.n_cols;
    this->q = yDat.n_cols;
    this->n = xDat.n_rows;

    std::vector<Variable*> cVar = ds.copyContinuousVariables();
    std::vector<Variable*> dVar = ds.copyDiscreteVariables();

    this->variables = std::vector<Variable*>();
    this->variables.reserve(p+q);
    this->variables.insert(this->variables.end(), cVar.begin(), cVar.end());
    this->variables.insert(this->variables.end(), dVar.begin(), dVar.end());
    
    this->initVariables = ds.copyVariables();

    this->lambda = arma::vec(lambda);

    //Data is checked for 0 or 1 indexing and for missing levels and N(0,1) Standardizes continuous data
    fixData();

    //Initialize all parameters to zeros
    initParameters();

    //Sets continuous variable weights to standard deviation and discrete variable weights to p*(1-p) for each category
    calcWeights();

    //Creates dummy variables for each category of discrete variables (stored in dDat)
    makeDummy();
}

MGM::~MGM() {
    for (int i = 0; i < variables.size(); i++) {
        delete variables[i];
        delete initVariables[i];
    }  
}

void MGM::initParameters() {
    lcumsum = std::vector<int>(l.size()+1);
    lcumsum[0] = 0;
    for(int i = 0; i < l.size(); i++){
        lcumsum[i+1] = lcumsum[i] + l[i];
    }
    lsum = lcumsum[l.size()];

    arma::mat beta((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros);

    params = MGMParams(
        arma::mat((int) xDat.n_cols, (int) xDat.n_cols, arma::fill::zeros),  // beta
        arma::vec(arma::size(xDat),                     arma::fill::ones),   // betad
        arma::mat(lsum, (int) xDat.n_cols,              arma::fill::zeros),  // theta
        arma::mat(lsum, lsum,                           arma::fill::zeros),  // phi
        arma::vec(arma::size(xDat),                     arma::fill::zeros),  // alpha1
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

/**
 * Convert discrete data (in yDat) to a matrix of dummy variables (stored in dDat)
 */
void MGM::makeDummy() {
    dDat = arma::mat(n, lsum, arma::fill::zeros);
    for(int i = 0; i < q; i++) {
        for(int j = 0; j < l[i]; j++) {
            arma::vec curCol = arma::vec(yDat.col(i)).transform( [j](double val) { return val == j+1 ? 1 : 0; } );
            if (arma::sum(curCol) == 0)
                throw std::invalid_argument("Discrete data is missing a level: variable " + variables[p+i]->getName() + " level " + std::to_string(j));
            dDat.col(lcumsum[i]+j) = curCol;
        }
    }
}

void MGM::fixData() {
    double ymin = yDat.min();
    if(ymin < 0 || ymin > 1)
        throw std::invalid_argument("Discrete data must be either zero or one indexed. Found min index: " + std::to_string(ymin));
    
    if (ymin == 0) {
        yDat += 1;
    }

    //TODO??
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
double MGM::smooth(arma::vec& parIn, arma::vec& gradOutVec) {
    MGMParams par(parIn, p, lsum);

    MGMParams gradOut;

    for(arma::uword i = 0; i < par.betad.size(); i++) {
        if(par.betad(i) <= 0)
            return std::numeric_limits<double>::infinity();
    }

    //beta=beta-diag(diag(beta));
    //for r=1:q
    //  phi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    //beta=triu(beta); phi=triu(phi);
    //beta=beta+beta';
    //phi=phi+phi';
    par.beta = arma::symmatu(par.beta);
    par.beta.diag(0).zeros();  // TODO: needed?

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
            tempLoss(i, j) = xDat(i,j) - par.alpha1(j) - xBeta(i,j) - dTheta(i,j);
        }
        for (arma::uword j = 0; j < dDat.n_cols; j++) {
            wxProd(i,j) = wxProd(i,j) + par.alpha2(j);
        }
    }

    //sqloss=-n/2*sum(log(betad))+...
    //.5*norm((X-e*alpha1'-Xbeta-Dtheta)*diag(sqrt(betad)),'fro')^2;
    double sqloss = -n*2.0*arma::sum(arma::log(arma::vec(par.betad))) +
                    0.5 * std::pow(arma::norm(tempLoss * arma::diagmat(arma::sqrt(arma::vec(par.betad))), "fro"), 2);

    //ok now tempLoss = res
    tempLoss *= -1;

    //gradbeta=X'*(res);
    gradOut.beta = xDat.t() * tempLoss;

    //gradbeta=gradbeta-diag(diag(gradbeta)); % zero out diag
    //gradbeta=tril(gradbeta)'+triu(gradbeta);
    // TODO - is this equivalent?
    gradOut.beta = arma::trimatu(gradOut.beta, 1) * 2;

    //gradalpha1=diag(betad)*sum(res,0)';
    gradOut.alpha1 = arma::diagmat(par.betad) * arma::sum(tempLoss, 0);

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
        arma::mat wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));

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
    gradOut.alpha2 = arma::sum(wxProd, 0);

    //gradw=X'*wxprod;
    arma::mat gradW = xDat.t() * wxProd;

    //gradtheta=gradtheta+gradw';
    gradOut.theta += gradW.t();

    //gradphi=D'*wxprod;
    gradOut.phi = dDat.t() * wxProd;

    //zero out gradphi diagonal
    //for r=1:q
    //gradphi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    for (arma::uword i = 0; i < q; i++) {
        gradOut.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    //gradphi=tril(gradphi)'+triu(gradphi);
    gradOut.phi = arma::trimatu(gradOut.phi, 0) * 2;
    // TODO - is this good?

    /*
    for s=1:p
        gradbetad(s)=-n/(2*betad(s))+1/2*norm(res(:,s))^2-res(:,s)'*(Xbeta(:,s)+Dtheta(:,s));
    end
        */
    gradOut.betad = arma::vec(xDat.n_cols);
    for(arma::uword i = 0; i < p; i++){
        gradOut.betad(i) = -n / (2.0 * par.betad(i)) + arma::norm(tempLoss.col(i), 2) / 2.0 -
                           arma::as_scalar(tempLoss.col(i) * (xBeta.col(i) + dTheta.col(i)));
    }

    gradOut.alpha1 /= (double) n;
    gradOut.alpha2 /= (double) n;
    gradOut.betad /= (double) n;
    gradOut.beta /= (double) n;
    gradOut.theta /= (double) n;
    gradOut.phi /= (double) n;

    gradOutVec = gradOut.toMatrix1D();
    return (sqloss + catloss)/((double) n);
}

double MGM::smoothValue(arma::vec& parIn) {
    MGMParams par(parIn, p, lsum);

    for(arma::uword i = 0; i < par.betad.size(); i++) {
        if(par.betad(i) <= 0)
            return std::numeric_limits<double>::infinity();
    }

    //double nll = 0;
    //int n = xDat.rows();
    //beta=beta+beta';
    //phi=phi+phi';
    par.beta = arma::symmatu(par.beta);
    par.beta.diag(0).zeros();  // TODO: needed?

    for (arma::uword i = 0; i < q; i++) {
        par.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    par.phi = arma::symmatu(par.phi);

    arma::mat divBetaD = arma::diagmat(arma::vec(p, arma::fill::ones) / par.betad);

    //Xbeta=X*beta*diag(1./betad);
    //Dtheta=D*theta*diag(1./betad);
    arma::mat xBeta = xDat * par.beta * divBetaD;
    arma::mat dTheta = dDat * par.theta * divBetaD;

    // Squared loss
    //sqloss=-n/2*sum(log(betad))+...
    //.5*norm((X-e*alpha1'-Xbeta-Dtheta)*diag(sqrt(betad)),'fro')^2;
    arma::mat tempLoss(n, xDat.n_cols);

    //wxprod=X*(theta')+D*phi+e*alpha2';
    arma::mat wxProd = xDat * par.theta.t() + dDat * par.phi;

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < xDat.n_cols; j++) {
            tempLoss(i, j) = xDat(i,j) - par.alpha1(j) - xBeta(i,j) - dTheta(i,j);
        }
        for (arma::uword j = 0; j < dDat.n_cols; j++) {
            wxProd(i,j) = wxProd(i,j) + par.alpha2(j);
        }
    }

    double sqloss = -n*2.0*arma::sum(arma::log(arma::vec(par.betad))) +
                    0.5 * std::pow(arma::norm(tempLoss * arma::diagmat(arma::sqrt(arma::vec(par.betad))), "fro"), 2);
    
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
        arma::mat wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));
        for (arma::uword k = 0; k < n; k++) {
            arma::vec curRow = wxTemp.row(k);

            catloss -= curRow((arma::uword) yDat(k,i) - 1);
            catloss += logsumexp(curRow);
        }
    }

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

    //penbeta = t(1).*(wv(1:p)'*wv(1:p));
    //betascale=zeros(size(beta));
    //betascale=max(0,1-penbeta./abs(beta));
    arma::mat weightMat = weights * weights.t();

    const arma::mat& betaWeight = weightMat.submat(0, 0, p, p);
    arma::mat betaScale = betaWeight * -tlam(0);
    arma::mat absBeta = arma::abs(par.beta); 
    
    betaScale /= absBeta;
    betaScale += 1;
    betaScale.transform( [](double val) {return std::max(val, 0.0); } );

    double betaNorms = 0;

    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < p; j++) {
            double curVal = par.beta(i,j);
            if (curVal != 0) {
                curVal *= betaScale(i,j);
                par.beta(i,j) = curVal;
                betaNorms += std::abs(betaWeight(i,j)*curVal);
            }
        }
    }

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
            const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]);
            double thetaScale = std::max(0.0, 1 - tlam(0)*weightMat(i,p+j)/arma::norm(tempVec, 2));
            par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]) = tempVec * thetaScale;
            thetaNorms += weightMat(i, p+j) * std::sqrt(arma::norm(tempVec, 2));
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
    double phiNorms = 0;
    for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
        for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
            const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1], lcumsum[j+1]);
            double phiScale = std::max(0.0, 1 - tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, 2));
            // Use the tempMat subview again to set the values (doesn't work with const)
            par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1], lcumsum[j+1]) = tempMat * phiScale;
            phiNorms += weightMat(p+i, p+j) * arma::norm(tempMat, "fro");
        }
    }

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
    double betaNorms = arma::accu(arma::mat(weightMat.submat(0, 0, p, p)) % arma::abs(par.beta));

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
           const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]);
           thetaNorms += weightMat(i, p+j) * std::sqrt(arma::norm(tempVec, 2));
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
   double phiNorms = 0;
   for (arma::uword i = 0; i < lcumsum.size()-1; i++) {
       for (arma::uword j = i+1; j < lcumsum.size()-1; j++) {
           const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1], lcumsum[j+1]);
           phiNorms += weightMat(p+i, p+j) * arma::norm(tempMat, "fro");
       }
   }

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

    par.beta = arma::symmatu(par.beta);
    par.beta.diag(0).zeros();  // TODO: needed?

    for (arma::uword i = 0; i < q; i++) {
        par.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    par.phi = arma::symmatu(par.phi);

    arma::mat divBetaD = arma::diagmat(arma::vec(p, arma::fill::ones) / par.betad);

    //Xbeta=X*beta*diag(1./betad);
    //Dtheta=D*theta*diag(1./betad);
    arma::mat xBeta = xDat * par.beta * divBetaD;
    arma::mat dTheta = dDat * par.theta * divBetaD;

    arma::mat negLoss(n, xDat.n_cols);

    arma::mat wxProd = xDat*par.theta.t() + dDat*par.phi;

    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < p; j++) {
            negLoss(i,j) = xBeta(i,j) - xDat(i,j) + par.alpha1(j) + dTheta(i,j);
        }
        for (arma::uword j = 0; j < dDat.n_cols; j++) {
            wxProd(i, j) += par.alpha2(j);
        }
    }

    //gradbeta=X'*(res);
    grad.beta = xDat.t() * negLoss;

    //gradbeta=gradbeta-diag(diag(gradbeta)); % zero out diag
    //gradbeta=tril(gradbeta)'+triu(gradbeta);
    // TODO - is this equivalent?
    grad.beta = arma::trimatu(grad.beta, 1) * 2;

    //gradalpha1=diag(betad)*sum(res,0)';
    grad.alpha1 = arma::diagmat(par.betad) * arma::sum(negLoss, 0);

    //gradtheta=D'*(res);
    grad.theta = dDat.t() * negLoss;

    for (arma::uword i = 0; i < yDat.n_cols; i++) {
        arma::mat wxTemp = wxProd(0, lcumsum[i], arma::size(n, l[i]));

        // does this need to be done in log space??
        wxTemp = arma::exp(wxTemp);
        arma::vec invDenom = arma::vec(n, arma::fill::ones) / arma::sum(wxTemp, 1);
        wxTemp = arma::diagmat(invDenom) * wxTemp;

        for (arma::uword k = 0; k < n; k++) {
            //wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))=wxtemp(sub2ind(size(wxtemp),(1:n)',Y(:,r)))-1;
            wxTemp(k, (arma::uword) yDat(k,i)-1) -= 1;
        }
    }

    //gradalpha2=sum(wxprod,0)';
    grad.alpha2 = arma::sum(wxProd, 0);

    //gradw=X'*wxprod;
    arma::mat gradW = xDat.t() * wxProd;

    //gradtheta=gradtheta+gradw';
    grad.theta += gradW.t();

    //gradphi=D'*wxprod;
    grad.phi = dDat.t() * wxProd;

    //zero out gradphi diagonal
    //for r=1:q
    //gradphi(Lsum(r)+1:Lsum(r+1),Lsum(r)+1:Lsum(r+1))=0;
    //end
    for (arma::uword i = 0; i < q; i++) {
        grad.phi(lcumsum[i], lcumsum[i], arma::size(l[i], l[i])).zeros();
    }

    //gradphi=tril(gradphi)'+triu(gradphi);
    grad.phi = arma::trimatu(grad.phi, 0) * 2;
    // TODO - is this good?

    /*
    for s=1:p
        gradbetad(s)=-n/(2*betad(s))+1/2*norm(res(:,s))^2-res(:,s)'*(Xbeta(:,s)+Dtheta(:,s));
    end
        */
    grad.betad = arma::vec(xDat.n_cols);
    for(arma::uword i = 0; i < p; i++){
        grad.betad(i) = -n / (2.0 * par.betad(i)) + arma::norm(negLoss.col(i), 2) / 2.0 -
                           arma::as_scalar(negLoss.col(i) * (xBeta.col(i) + dTheta.col(i)));
    }

    grad.alpha1 /= (double) n;
    grad.alpha2 /= (double) n;
    grad.betad /= (double) n;
    grad.beta /= (double) n;
    grad.theta /= (double) n;
    grad.phi /= (double) n;

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

    const arma::mat& betaWeight = weightMat.submat(0, 0, p, p);
    arma::mat betaScale = betaWeight * -tlam(0);
    
    betaScale /= arma::abs(par.beta);
    betaScale += 1;
    betaScale.transform( [](double val) {return std::max(val, 0.0); } );

    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < p; j++) {
            double curVal = par.beta(i,j);
            if (curVal != 0) {
                par.beta(i,j) = curVal*betaScale(i,j);
            }
        }
    }

    // TODO - same as par.beta = par.beta % betaScale ??

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
    for (arma::uword i = 0; i < p; i++) {
        for (arma::uword j = 0; j < lcumsum.size()-1; j++) {
            const arma::vec& tempVec = par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]);
            double thetaScale = std::max(0.0, 1 - tlam(0)*weightMat(i,p+j)/arma::norm(tempVec, 2));
            par.theta.col(i).subvec(lcumsum[j], lcumsum[j+1]) = tempVec * thetaScale;
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
            const arma::mat& tempMat = par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1], lcumsum[j+1]);
            double phiScale = std::max(0.0, 1 - tlam(2)*weightMat(p+i,p+j)/arma::norm(tempMat, 2));
            // Use the tempMat subview again to set the values (doesn't work with const)
            par.phi.submat(lcumsum[i], lcumsum[j], lcumsum[i+1], lcumsum[j+1]) = tempMat * phiScale;
        }
    }

    return par.toMatrix1D();
}
