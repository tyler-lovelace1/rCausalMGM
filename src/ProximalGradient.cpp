#include "ProximalGradient.hpp"

ProximalGradient::ProximalGradient(double beta, double alpha, bool edgeConverge) {
    if(beta <= 0 || beta >=1)
        throw std::invalid_argument("beta must be (0,1): " + std::to_string(beta));

    if(alpha <= 0 || alpha >=1)
        throw std::invalid_argument("alpha must be (0,1): " + std::to_string(alpha));

    this->beta = beta;
    this->alpha = alpha;
    this->edgeConverge = edgeConverge;
}

arma::vec ProximalGradient::learnBackTrack(ConvexProximal *cp, arma::vec& Xin, double epsilon, int iterLimit, long time) {
    auto start = std::chrono::high_resolution_clock::now();

    arma::vec X(Xin); // = cp->proximalOperator(1.0, Xin);
    arma::vec Y(X);
    arma::vec Z(X);
    arma::vec Xold(X);

    cp->iterUpdate(X);

    arma::vec GrY(cp->smoothGradient(Y));
    arma::vec GrX(cp->smoothGradient(X));

    int iterCount = 0;
    int noEdgeChangeCount = 0;

    double theta = std::numeric_limits<double>::infinity();
    double thetaOld = theta;
    double L = 1.0;
    double Lold = L;

    bool backtrackSwitch = true;
    double dx;
    double Fx = std::numeric_limits<double>::infinity();
    double Gx = std::numeric_limits<double>::infinity();
    double Fy;
    double obj;

    while(true) {
        auto lastStart = std::chrono::high_resolution_clock::now();
        Lold = L;
        L = L*alpha;
        thetaOld = theta;
	cp->iterUpdate(Xold);
        Xold = arma::vec(X);
	obj = Fx + Gx;

        while(true) {
            theta = 2.0 / (1.0 + std::sqrt(1.0+(4.0*L)/(Lold*std::pow(thetaOld,2))));

            if (theta < 1) {
                Y = Xold * (1 - theta);
                Y += Z * theta;
            }

            Fy = cp->smooth(Y, GrY);
            arma::vec temp = Y - (GrY * (1.0/L));
            Gx = cp->nonSmooth(1.0 / L, temp, X);
            if (backtrackSwitch) {
                Fx = cp->smoothValue(X);
            } else {
                Fx = cp->smooth(X, GrX);
            }

            arma::vec XmY = X - Y;

            double normXY = std::pow(arma::norm(XmY, 2), 2);
            if (normXY == 0)
                break;

            double Qx;
            double LocalL;

            if (backtrackSwitch) {
                Qx = Fy + arma::as_scalar(XmY.t() * GrY) + (L / 2.0) * normXY;
                LocalL = L + 2*std::max(Fx - Qx, 0.0)/normXY;
                backtrackSwitch = std::abs(Fy - Fx) >= backtrackTol * std::max(std::abs(Fx), std::abs(Fy));
            } else {
                LocalL = 2 * arma::as_scalar(XmY.t() * (GrX - GrY)) / normXY;
            }

            if (LocalL <= L) {
                break;
            } else if (LocalL != std::numeric_limits<double>::infinity()) {
                L = LocalL;
            } else {
                LocalL = L;
            }

            L = std::max(LocalL, L/beta);
        }

        int diffEdges = 0;
        for (int i = 0; i < X.n_elem; i++) {
            double a = X(i);
            double b = Xold(i);
            if ((a == 0 || b == 0) && a != b) {
                diffEdges++;
            }
        }

	dx = arma::norm((X - Xold), 2) / std::max(1.0, arma::norm(X, 2));

        //sometimes there are more edge changes after initial 0, so may want to do two zeros in a row...
        if (diffEdges == 0 && edgeConverge) {
            noEdgeChangeCount++;
            if (noEdgeChangeCount >= noEdgeChangeTol) {
                break;
            } 
        } else if (noEdgeChangeTol < 0 && diffEdges <= std::abs(noEdgeChangeTol)) {
            break;
        } else {
            noEdgeChangeCount = 0;
        }

        //edge converge should happen before params converge, unless epsilon is big
        if (dx < epsilon && !edgeConverge) {
            break;
        }

        //restart acceleration if objective got worse
        if (Fx + Gx > obj) {
            theta = std::numeric_limits<double>::infinity();
            Y = X;
            Z = X;
        } else if (theta == 1) {
            Z = X;
        } else {
            Z = X * (1/theta);
            Z += Xold * (1 - (1.0 / theta));
        }

        printIter = 10;
        if (iterCount % printIter == 0) {
	    if (RcppThread::isInterrupted()) {
		break;
	    }
            // Rcpp::Rcout << "Iter: " << iterCount << 
            //     " |dx|/|x|: " << dx << 
            //     " normX: " << arma::norm(X, 2) << 
            //     " Fx: " << Fx << 
            //     " Fy: " << Fy <<
            //     " Gx: " << Gx << 
            //     " DiffEdges: " << diffEdges << 
            //     " L: " << L << 
            //     " theta: " << theta << std::endl;
            // Rcpp::Rcout << "X = \n" << X.t() << std::endl;
        }

        iterCount++;
        if (iterCount >= iterLimit) {
            break;
        }

        auto finish = std::chrono::high_resolution_clock::now();
        timePerIter += std::chrono::duration_cast<std::chrono::seconds>(finish-start).count();
        iterComplete++;
        if (std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() >= time) {
            return X;
        }
    }
    return X;
 }

arma::vec ProximalGradient::learnBackTrack(ConvexProximal *cp, arma::vec& Xin, double epsilon, int iterLimit) {
    arma::vec X(Xin); // = cp->proximalOperator(1.0, Xin);

    // Rcpp::Rcout << "Original X = \n" << X.t() << std::endl;

    arma::vec Y(X);
    arma::vec Z(X);
    arma::vec Xold(X);

    cp->iterUpdate(X);

    arma::vec GrY(cp->smoothGradient(Y));
    arma::vec GrX(cp->smoothGradient(X));

    int iterCount = 0;
    int noEdgeChangeCount = 0;
    int increaseObjCount = 0;

    double theta = std::numeric_limits<double>::infinity();
    double thetaOld = theta;
    double L = 1.0;
    double Lcur = L;
    double Lold = L;
    double LocalL;

    bool backtrackSwitch = true;
    double dx;
    double Fx = std::numeric_limits<double>::infinity();
    double Gx = std::numeric_limits<double>::infinity();
    double Gxold;
    double Fxold;
    double Fy;
    double obj;
    // double minObj = std::numeric_limits<double>::infinity();
    // arma::vec minX = arma::vec(X);

    while(true) {
        Lold = L;
        L = L*alpha;
        thetaOld = theta;
	Gxold = Gx;
	Fxold = Fx;
	cp->iterUpdate(Xold);
        Xold = arma::vec(X);
	obj = Fx + Gx;
	
        while(true) {
            theta = 2.0 / (1.0 + std::sqrt(1.0+(4.0*L)/(Lold*std::pow(thetaOld,2))));

            if (theta < 1) {
                Y = Xold * (1 - theta);
                Y += Z * theta;
            }

            Fy = cp->smooth(Y, GrY);
            arma::vec temp = Y - (GrY * (1.0/L));

            // Rcpp::Rcout << "Y: " << Y.t() << std::endl;
            // Rcpp::Rcout << "GrY: " << GrY.t() << std::endl;
            // MGMParams tempParams(temp, 5, 20);
            // Rcpp::Rcout << "tempParams: \n" << tempParams << std::endl;

	    Lcur = L;

            Gx = cp->nonSmooth(1.0 / L, temp, X);
            if (backtrackSwitch) {
                Fx = cp->smoothValue(X);
            } else {
                Fx = cp->smooth(X, GrX);
            }

            arma::vec XmY = X - Y;

            double normXY = std::pow(arma::norm(XmY, 2), 2);
            if (normXY == 0)
                break;

            double Qx = 0.0;
            // double LocalL;
	    
            if (backtrackSwitch) {
                Qx = Fy + arma::as_scalar(XmY.t() * GrY) + (L / 2.0) * normXY;
                LocalL = L + 2*std::max(Fx - Qx, 0.0)/normXY;
                backtrackSwitch = std::abs(Fy - Fx) >= backtrackTol * std::max(std::abs(Fx), std::abs(Fy));
            } else {
                LocalL = 2 * arma::as_scalar(XmY.t() * (GrX - GrY)) / normXY;
            }

            // Rcpp::Rcout << "   LocalL: " << LocalL << " L: " << L << std::endl;
	    // Rcpp::Rcout << "Qx: " << Qx << " Fx: " << Fx << std::endl;
            if (LocalL <= L) {
		// L = (LocalL/alpha < 1) ? LocalL/alpha : L;
		break;
            } else if (LocalL != std::numeric_limits<double>::infinity()) {
                L = LocalL;
            } else {
                LocalL = L;
            }

	    // if (Qx - Fx <= backtrackTol * Fx) {
	    // 	break;
	    // }

            L = std::max(LocalL, L/beta);
        }

	// L = std::min(10*LocalL, L);

        int diffEdges = 0;
        for (int i = 0; i < X.n_elem; i++) {
            double a = X(i);
            double b = Xold(i);
            if ((a == 0 || b == 0) && a != b) {
                diffEdges++;
            }
        }

	dx = arma::norm((X - Xold), 2) / std::max(1.0, arma::norm(X, 2));

	// dx = (Fx+Gx - obj) / std::min(Fx+Gx, obj);

	// if (dx > 0) {
	//     Rcpp::Rcout << "Objective function increase at Iter: " << iterCount << 
	// 	    " |dx|/|x|: " << dx << 
	// 	    " Fx: " << Fx << 
	// 	    " Gx: " << Gx <<
	// 	    " loss: " << Fx + Gx <<
	// 	    " DiffEdges: " << diffEdges << 
	// 	    " L: " << L << 
	// 	    " theta: " << theta << std::endl;
	// }

        // sometimes there are more edge changes after initial 0, so may want to do two zeros in a row...
	if (edgeConverge) {
	    if (diffEdges == 0) {
		noEdgeChangeCount++;
		if (noEdgeChangeCount >= noEdgeChangeTol) {
		    break;
		} 
	    } else if (noEdgeChangeTol < 0 && diffEdges <= std::abs(noEdgeChangeTol)) {
		break;
	    } else {
		noEdgeChangeCount = 0;
	    }
	} else {
	    if (diffEdges == 0) {
		noEdgeChangeCount++;
	    } else {
		noEdgeChangeCount = 0;
	    }
	}

        //edge converge should happen before params converge, unless epsilon is big
        if (dx < epsilon && L < 1/epsilon && !edgeConverge && noEdgeChangeCount >= noEdgeChangeTol) {
            break;
        }

	// if (increaseObjCount > 10) {
	//     break;
	// }

	// Rcpp::Rcout << "restart if "
	// 	    << arma::as_scalar((X-Xold).t() * GrY) + Gx - Gxold
	// 	    << " > 0\n";
        //restart acceleration if objective got worse
	// if (arma::as_scalar((X-Xold).t() * GrY) + Gx - Gxold > 0) {
	if (Fx + Gx > obj) {
	    // Rcpp::Rcout << "      restarting acceleration: "
	    // 		<< Fx + Gx << " > " << obj << "\n";
            theta = std::numeric_limits<double>::infinity();
            Y = X;
            Z = X;
	    // Y = Xold;
	    // Z = X * (1.0/theta);
            // Z += Xold * (1 - (1.0 / theta));
        } else if (theta == 1) {
            Z = X;
        } else {
            Z = X * (1.0/theta);
            Z += Xold * (1 - (1.0 / theta));
        }

	// if (arma::as_scalar((X-Xold).t() * GrY) + Gx - Gxold > 0) {
	//     theta = std::numeric_limits<double>::infinity();
	// }

        printIter = 10;
        if (iterCount % printIter == 0) {
	    if (RcppThread::isInterrupted()) {
		break;
	    }
	    
	    if (cp->isVerbose()) {
		RcppThread::Rcout << "\r    Iter: " << iterCount << 
		    " ||dx||/||x||: " << dx << 
		    " loss: " << Fx + Gx;
		// RcppThread::Rcout << "Params = \n" << cp->printParameters(X) << std::endl;
	    }
        }

	L = (L > 1/epsilon || L/LocalL > 100) ? std::max(LocalL/alpha, std::sqrt(L)) : L;

	// if (obj < minObj) {
	//     increaseObjCount = 0;
	//     minObj = obj;
	//     minX = X;
	// } else {
	//     increaseObjCount++;
	// }

        iterCount++;
        if (iterCount >= iterLimit) {
            break;
        }
    }

    if (cp->isVerbose()) {
	RcppThread::Rcout << "\r    Iter: " << iterCount << 
	    " ||dx||/||x||: " << dx << 
	    " loss: " << Fx + Gx << "\n";
	// Rcpp::Rcout << "X = \n" << X.t() << std::endl;
    }

    // MGMParams XP(X, 5, 20);
    // Rcpp::Rcout << "XP: \n" << XP << std::endl;

    return X;
 }
