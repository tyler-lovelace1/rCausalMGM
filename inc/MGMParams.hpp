#ifndef MGMPARAMS_HPP_
#define MGMPARAMS_HPP_

#include <RcppArmadillo.h>

class MGMParams {

private:
    arma::mat beta;
    arma::vec betad;
    arma::mat theta;
    arma::mat phi;
    arma::vec alpha1;
    arma::vec alpha2;

public:
    MGMParams() {}

    MGMParams(arma::mat& beta, arma::mat& betad, arma::mat& theta, arma::mat& phi, arma::vec& alpha1, arma::vec& alpha2);

    MGMParams(MGMParams& parIn);

    MGMParams(arma::vec& vec, int p, int ltot);

    arma::vec& toMatrix1D();

    arma::mat& getBeta() {
        return beta;
    }

    arma::vec& getBetad() {
        return betad;
    }

    arma::mat& getTheta() {
        return theta;
    }

    arma::mat& getPhi() {
        return phi;
    }

    arma::vec& getAlpha1() {
        return alpha1;
    }

    arma::vec& getAlpha2() {
        return alpha2;
    }

    void setBeta(arma::mat& beta) {
        this->beta = beta;
    }

    void setBetad(arma::vec& betad) {
        this->betad = betad;
    }

    void setTheta(arma::mat& theta) {
        this->theta = theta;
    }

    void setPhi(arma::mat& phi) {
        this->phi = phi;
    }

    void setAlpha1(arma::vec& alpha1) {
        this->alpha1 = alpha1;
    }

    void setAlpha2(arma::vec& alpha2) {
        this->alpha2 = alpha2;
    }
};

#endif /* MGMPARAMS_HPP_ */
