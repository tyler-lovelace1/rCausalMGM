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

    arma::vec X = arma::vec(Xin); // cp->proximalOperator(1.0, Xin);
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

    bool backtrackSwitch = false;
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
        obj = Fx + Gx;

        while(true) {
	    if (!std::isinf(thetaOld))
		theta = 2.0 / (1.0 + std::sqrt(1.0+(4.0*L) / (Lold*std::pow(thetaOld,2))));
	    else
		theta = 1.0;

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

            L = std::max(1.0, L/beta);
        }

        int diffEdges = 0;
        for (int i = 0; i < X.n_elem; i++) {
            double a = X(i);
            double b = Xold(i);
            if ((a == 0 || b == 0) && a != b) {
                diffEdges++;
            }
        }

        dx = arma::norm(X - Xold, 2) / std::max(1e-10, arma::norm(X, 2));

        printIter = 1;
        if (iterCount % printIter == 0) {
            Rcpp::Rcout << "Iter: " << iterCount <<
		" obj: " << obj << 
                " |dx|/|x|: " << dx << 
                " normX: " << arma::norm(X, 2) << 
                " Fx: " << Fx << 
                " Fy: " << Fy <<
                " Gx: " << Gx << 
                " DiffEdges: " << diffEdges << 
                " L: " << L << 
                " theta: " << theta <<
		"\n GrY: " << GrY.t() << std::endl;
            // Rcpp::Rcout << "X = \n" << X.t() << std::endl;
        }

        auto finish = std::chrono::high_resolution_clock::now();
        timePerIter += std::chrono::duration_cast<std::chrono::milliseconds>(finish-lastStart).count();
        iterComplete++;

	iterCount++;
        if (iterCount >= iterLimit) {
            break;
        }
	
        if (std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() >= time) {
	    break;
        }

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
	
    }
    timePerIter /= iterComplete;
    return X;
 }

arma::vec ProximalGradient::learnBackTrack(ConvexProximal *cp, arma::vec& Xin, double epsilon, int iterLimit) {
    arma::vec X = arma::vec(Xin); // cp->proximalOperator(1.0, Xin);

    // Rcpp::Rcout << "Original X = \n" << X.t() << std::endl;

    arma::vec Y = arma::vec(X);
    arma::vec Z = arma::vec(X);

    arma::vec XmY;

    arma::vec GrY = cp->smoothGradient(Y);
    arma::vec GrX = cp->smoothGradient(X);

    int iterCount = 0;
    int noEdgeChangeCount = 0;

    double theta = std::numeric_limits<double>::infinity();
    double thetaOld = theta;
    double L = 1.0 / alpha;
    double Lold = L;

    bool backtrackSwitch = false;
    double dx;
    double Fx = std::numeric_limits<double>::infinity();
    double Gx = std::numeric_limits<double>::infinity();
    double Fy;
    double obj;
    double normXY;

    while(true) {
	auto lastStart = std::chrono::high_resolution_clock::now();
        Lold = L;
        L = L*alpha;
        thetaOld = theta;
        arma::vec Xold = arma::vec(X);
        obj = Fx + Gx;

        while(true) {
	    if (!std::isinf(thetaOld))
		theta = 2.0 / (1.0 + std::sqrt(1.0+(4.0*L) / (Lold*std::pow(thetaOld,2))));
	    else
		theta = 1.0;
            // theta = 2.0 / (1.0 + std::sqrt(1.0+(4.0*L)/(Lold*std::pow(thetaOld,2))));

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

            Gx = cp->nonSmooth(1.0 / L, temp, X);
            if (backtrackSwitch) {
                Fx = cp->smoothValue(X);
            } else {
		// Rcpp::Rcout << " (X - Y):\n" << (X.subvec(150, 199) - Y.subvec(150, 199)).t() << std::endl;
                Fx = cp->smooth(X, GrX);
            }

            XmY = X - Y;

            normXY = std::pow(arma::norm(XmY, 2), 2);
            if (normXY == 0)
                break;

            double Qx;
            double LocalL;

            if (backtrackSwitch) {
                Qx = Fy + arma::as_scalar(XmY.t() * GrY) + (L / 2.0) * normXY;
                LocalL = L + 2 * std::max(Fx - Qx, 0.0) / normXY;
                backtrackSwitch = std::abs(Fy - Fx) >= backtrackTol * std::max(std::abs(Fx), std::abs(Fy));
            } else {
                LocalL = 2 * arma::as_scalar(XmY.t() * (GrX - GrY)) / normXY;
            }

            // Rcpp::Rcout << "   LocalL: " << LocalL << " L: " << L << std::endl;
	    // Rcpp::Rcout << "   backtrackSwitch: " << backtrackSwitch << std::endl;
	    
            if (LocalL <= L) {
                break;
            } else if (LocalL != std::numeric_limits<double>::infinity()) {
                L = LocalL;
            } else {
                LocalL = L;
            }

	    // Rcpp::Rcout << "   LocalL: " << LocalL << " L: " << L << std::endl;

            L = std::max(1.0, L/beta);
        }

        int diffEdges = 0;
        for (int i = 0; i < X.n_elem; i++) {
            double a = X(i);
            double b = Xold(i);
            if ((a == 0 || b == 0) && a != b) {
                diffEdges++;
            }
        }

        dx = arma::norm(X - Xold, 2) / std::max(1e-10, arma::norm(X, 2));

        //sometimes there are more edge changes after initial 0, so may want to do two zeros in a row..

        printIter = 1;
        if (iterCount % printIter == 0) {
            // Rcpp::Rcout << "Iter: " << iterCount <<
	    // 	" obj: " << obj << 
            //     " |dx|/|x|: " << dx << 
            //     " normX: " << arma::norm(X, 2) << 
            //     " Fx: " << Fx << 
            //     " Fy: " << Fy <<
            //     " Gx: " << Gx << 
            //     " DiffEdges: " << diffEdges << 
            //     " L: " << L << 
            //     " theta: " << theta <<
	    // 	" normXY: " << normXY << std::endl;
		// "\n GrXmGrY: " << (GrX.subvec(150, 199) - GrY.subvec(150, 199)).t() <<
		// "\n XmY: " << XmY.subvec(150, 199).t() << std::endl;
            // Rcpp::Rcout << "X = \n" << X.t() << std::endl;
        }

        iterCount++;
	auto finish = std::chrono::high_resolution_clock::now();
        timePerIter += std::chrono::duration_cast<std::chrono::milliseconds>(finish-lastStart).count();
        iterComplete++;

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
	
        if (iterCount >= iterLimit) {
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
    }

    // MGMParams XP(X, 5, 20);
    // Rcpp::Rcout << "XP: \n" << XP << std::endl;

    timePerIter /= iterComplete;

    return X;
 }
