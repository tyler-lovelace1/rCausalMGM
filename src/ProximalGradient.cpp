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

    // RcppThread::Rcout << "Original X = \n" << Xin.t() << std::endl;
    // RcppThread::Rcout << "Proximal Operator X = \n" << X.t() << std::endl;
    
    arma::vec Y = arma::vec(X);
    arma::vec Z = arma::vec(X);

    arma::vec GrY = cp->smoothGradient(Y);
    arma::vec GrX = cp->smoothGradient(X);

    int iterCount = 0;
    int noEdgeChangeCount = 0;

    double theta = std::numeric_limits<double>::infinity();
    double thetaOld = theta;
    double L = 1.0;
    double Lold = L;
    double LocalL;
    double stepSize;

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
        arma::vec Xold = arma::vec(X);
	cp->iterUpdate(Xold);
        obj = Fx + Gx;

        while(true) {
            theta = 2.0 / (1.0 + std::sqrt(1.0+(4.0*L)/(Lold*std::pow(thetaOld,2))));

            if (theta < 1) {
                Y = Xold * (1 - theta);
                Y += Z * theta;
            }

            Fy = cp->smooth(Y, GrY);

	    stepSize = (1.0/L);
	    
            arma::vec temp = Y - (GrY * stepSize);
            Gx = cp->nonSmooth(stepSize, temp, X);
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
            // double LocalL;

            if (backtrackSwitch) {
                Qx = Fy + arma::as_scalar(XmY.t() * GrY) + (L / 2.0) * normXY;
                LocalL = L + 2*std::max(Fx - Qx, 0.0)/normXY;
                backtrackSwitch = std::abs(Fy - Fx) >= backtrackTol * std::max(std::abs(Fx), std::abs(Fy));
            } else {
                LocalL = 2 * arma::as_scalar(XmY.t() * (GrX - GrY)) / normXY;
            }

            if (LocalL <= L) {
		if (L/LocalL > 5) L = std::min(LocalL/beta, L);
                break;
            } else if (LocalL == std::numeric_limits<double>::infinity()) {
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

        dx = arma::norm((X - Xold) / std::max(1.0, arma::norm(X, 2)), 2);

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
        if (dx < epsilon && !edgeConverge && (1/stepSize)/LocalL < 10) {
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

        printIter = 5;
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

    // RcppThread::Rcout << "Original X = \n" << Xin.t() << std::endl;
    // RcppThread::Rcout << "Proximal Operator X = \n" << X.t() << std::endl;

    arma::vec Y = arma::vec(X);
    arma::vec Z = arma::vec(X);

    cp->iterUpdate(X);

    arma::vec GrY = cp->smoothGradient(Y);
    arma::vec GrX = cp->smoothGradient(X);

    int iterCount = 0;
    int noEdgeChangeCount = 0;

    double theta = std::numeric_limits<double>::infinity();
    double thetaOld = theta;
    double L = 1.0;
    double Lold = L;
    double LocalL;
    double stepSize;

    bool backtrackSwitch = true;
    double dx;
    double Fx = std::numeric_limits<double>::infinity();
    double Gx = std::numeric_limits<double>::infinity();
    double Fy;
    double obj;

    while(true) {
        Lold = L;
        L = L*alpha;
        thetaOld = theta;
        arma::vec Xold = arma::vec(X);
	cp->iterUpdate(Xold);
        obj = Fx + Gx;

        while(true) {
            theta = 2.0 / (1.0 + std::sqrt(1.0+(4.0*L)/(Lold*std::pow(thetaOld,2))));

            if (theta < 1) {
                Y = Xold * (1 - theta);
                Y += Z * theta;
            }

            Fy = cp->smooth(Y, GrY);
	    
	    stepSize = (1.0/L);
	    
            arma::vec temp = Y - (GrY * stepSize);

            // Rcpp::Rcout << "Y: " << Y.t() << std::endl;
            // Rcpp::Rcout << "GrY: " << GrY.t() << std::endl;
            // MGMParams tempParams(temp, 5, 20);
            // Rcpp::Rcout << "tempParams: \n" << tempParams << std::endl;

            Gx = cp->nonSmooth(stepSize, temp, X);
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
            // double LocalL;

            if (backtrackSwitch) {
                Qx = Fy + arma::as_scalar(XmY.t() * GrY) + (L / 2.0) * normXY;
                LocalL = L + 2*std::max(Fx - Qx, 0.0)/normXY;
                backtrackSwitch = std::abs(Fy - Fx) >= backtrackTol * std::max(std::abs(Fx), std::abs(Fy));
            } else {
                LocalL = 2 * arma::as_scalar(XmY.t() * (GrX - GrY)) / normXY;
            }

            // RcppThread::Rcout << "   LocalL: " << LocalL << " L: " << L << std::endl;
            // if (LocalL <= L) {
	    // 	// RcppThread::Rcout << "   L/LocalL ratio: " << L/LocalL << std::endl;
	    // 	if (L/LocalL > 5) L = std::min(LocalL/beta, L);
            //     break;
            // } else if (LocalL == std::numeric_limits<double>::infinity()) {
            //     LocalL = L;
            // }

	    if (LocalL <= L) {
	        // L = std::min(LocalL/beta, L);
		L = (LocalL/alpha < 1) ? LocalL/alpha : L;
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

        dx = arma::norm((X - Xold) / std::max(1.0, arma::norm(X, 2)), 2);

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

	// RcppThread::Rcout << "   L/LocalL ratio: " << (1/stepSize)/LocalL << std::endl;
        //edge converge should happen before params converge, unless epsilon is big
        if (dx < epsilon && !edgeConverge && (1/stepSize)/LocalL < 10) {
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

	    RcppThread::Rcout << "    Iter: " << iterCount << 
		" |dx|/|x|: " << dx << 
		" loss: " << Fx + Gx << "\n";
		

	    // RcppThread::Rcout << "  Iter: " << iterCount << 
            //     "\n    obj: " << obj <<
            //     "\n    Fx: " << Fx << 
            //     "\n    Fy: " << Fy <<
            //     "\n    Gx: " << Gx << 
            //     "\n    |dx|/|x|: " << dx <<  
            //     "\n    DiffEdges: " << diffEdges << 
            //     "\n    L: " << L << 
            //     "\n    theta: " << theta << std::endl;
            // Rcpp::Rcout << "X = \n" << X.t() << std::endl;
        }

        iterCount++;
        if (iterCount >= iterLimit) {
            break;
        }

    }

    RcppThread::Rcout << "    Iter: " << iterCount << 
	" |dx|/|x|: " << dx << 
	" loss: " << Fx + Gx << "\n";
	

    // MGMParams XP(X, 5, 20);
    // Rcpp::Rcout << "XP: \n" << XP << std::endl;

    return X;
 }
