#ifndef MGMPARAMS_HPP_
#define MGMPARAMS_HPP_

#include "armaLapack.hpp"

class MGMParams {

private:
    std::vector<std::string> names;
    arma::mat beta;   //continuous-continuous
    arma::vec betad;  //cont squared node pot
    arma::mat theta;  //continuous-discrete
    arma::mat phi;    //discrete-discrete
    arma::vec alpha1; //cont linear node pot
    arma::vec alpha2; //disc node pot

    friend class MGM;

public:
    MGMParams() {}
    MGMParams(int p, int ltot);
    MGMParams(const arma::mat& beta, const arma::vec& betad, const arma::mat& theta,
	      const arma::mat& phi, const arma::vec& alpha1, const arma::vec& alpha2);
    // MGMParams(const std::vector<std::string>& names, const arma::mat& beta,
    // 	      const arma::vec& betad, const arma::mat& theta, const arma::mat& phi,
    // 	      const arma::vec& alpha1, const arma::vec& alpha2);
    // MGMParams(MGMParams& parIn);
    MGMParams(arma::vec& vec, int p, int ltot);  //copy params from flattened vector

    // MGMParams(const MGMParams& other) = default;
    // MGMParams& operator=(const MGMParams& other) = default;
    // MGMParams(const MGMParams&& other) = default;
    // MGMParams& operator=(const MGMParams&& other) = default;
    // ~MGMParams() = default;


    /**
     * Copy all params into a single vector
     */
    arma::vec toMatrix1D();

    arma::mat getBeta() { return beta; }
    arma::vec getBetad() { return betad; }
    arma::mat getTheta() { return theta; }
    arma::mat getPhi() { return phi; }
    arma::vec getAlpha1() { return alpha1; }
    arma::vec getAlpha2() { return alpha2; }
    void setBeta(arma::mat& beta) { this->beta = beta; }
    void setBetad(arma::vec& betad) { this->betad = betad; }
    void setTheta(arma::mat& theta) { this->theta = theta; }
    void setPhi(arma::mat& phi) { this->phi = phi; }
    void setAlpha1(arma::vec& alpha1) { this->alpha1 = alpha1; }
    void setAlpha2(arma::vec& alpha2) { this->alpha2 = alpha2; }

    void setNames(std::vector<std::string>& names) { this->names = names; }
    std::vector<std::string> getNames() { return names; }

    Rcpp::List toList() const;

    friend std::ostream& operator<<(std::ostream& os, MGMParams params);
};

#endif /* MGMPARAMS_HPP_ */
