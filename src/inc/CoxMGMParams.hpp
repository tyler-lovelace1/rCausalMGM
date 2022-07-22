#ifndef COXMGMPARAMS_HPP_
#define COXMGMPARAMS_HPP_

#include "armaLapack.hpp"

class CoxMGMParams {

private:
    arma::mat beta;   //continuous-continuous   p x p
    arma::vec betad;  //cont squared node pot   p x 1
    arma::mat theta;  //continuous-discrete     ltot x p
    arma::mat phi;    //discrete-discrete       ltot x ltot
    arma::mat gamma;  //cox-continuous          p x r
    arma::mat psi;    //cox-discrete            ltot x r
    arma::vec alpha1; //cont linear node pot    p x 1
    arma::vec alpha2; //disc node pot           ltot x 1
    arma::vec alpha3; //cox node pot            r x 1

    friend class CoxMGM;

public:
    CoxMGMParams() {}
    CoxMGMParams(const arma::mat& beta, const arma::vec& betad, const arma::mat& theta,
		 const arma::mat& phi, const arma::mat& gamma, const arma::mat& psi,
		 const arma::vec& alpha1, const arma::vec& alpha2, const arma::vec& alpha3);
    // CoxMGMParams(CoxMGMParams& parIn);
    CoxMGMParams(arma::vec& vec, int p, int ltot, int r);  //copy params from flattened vector

    /**
     * Copy all params into a single vector
     */
    arma::vec toMatrix1D();

    arma::mat& getBeta() { return beta; }
    arma::vec& getBetad() { return betad; }
    arma::mat& getTheta() { return theta; }
    arma::mat& getPhi() { return phi; }
    arma::mat& getGamma() { return gamma; }
    arma::mat& getPsi() { return psi; }
    arma::vec& getAlpha1() { return alpha1; }
    arma::vec& getAlpha2() { return alpha2; }
    arma::vec& getAlpha3() { return alpha3; }
    void setBeta(arma::mat& beta) { this->beta = beta; }
    void setBetad(arma::vec& betad) { this->betad = betad; }
    void setTheta(arma::mat& theta) { this->theta = theta; }
    void setPhi(arma::mat& phi) { this->phi = phi; }
    void setGamma(arma::mat& gamma) { this->gamma = gamma; }
    void setPsi(arma::mat& psi) { this->psi = psi; }
    void setAlpha1(arma::vec& alpha1) { this->alpha1 = alpha1; }
    void setAlpha2(arma::vec& alpha2) { this->alpha2 = alpha2; }
    void setAlpha3(arma::vec& alpha3) { this->alpha3 = alpha3; }

    friend std::ostream& operator<<(std::ostream& os, CoxMGMParams& params);
};

#endif /* COXMGMPARAMS_HPP_ */
