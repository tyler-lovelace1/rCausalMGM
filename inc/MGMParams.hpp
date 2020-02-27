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

    friend class MGM;

public:
    MGMParams() {}
    MGMParams(const arma::mat& beta, const arma::vec& betad, const arma::mat& theta, const arma::mat& phi, const arma::vec& alpha1, const arma::vec& alpha2);
    MGMParams(MGMParams& parIn);
    MGMParams(arma::vec& vec, int p, int ltot);

    /**
     * Copy all params into a single vector
     */
    arma::vec toMatrix1D();

    arma::mat& getBeta() { return beta; }
    arma::vec& getBetad() { return betad; }
    arma::mat& getTheta() { return theta; }
    arma::mat& getPhi() { return phi; }
    arma::vec& getAlpha1() { return alpha1; }
    arma::vec& getAlpha2() { return alpha2; }
    void setBeta(arma::mat& beta) { this->beta = beta; }
    void setBetad(arma::vec& betad) { this->betad = betad; }
    void setTheta(arma::mat& theta) { this->theta = theta; }
    void setPhi(arma::mat& phi) { this->phi = phi; }
    void setAlpha1(arma::vec& alpha1) { this->alpha1 = alpha1; }
    void setAlpha2(arma::vec& alpha2) { this->alpha2 = alpha2; }

    friend std::ostream& operator<<(std::ostream& os, MGMParams& params);
};

#endif /* MGMPARAMS_HPP_ */
