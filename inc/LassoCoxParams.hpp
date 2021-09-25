#ifndef LASSOCOXPARAMS_HPP_
#define LASSOCOXPARAMS_HPP_

#include <RcppArmadillo.h>

class LassoCoxParams {

private:
    arma::mat gamma;  //continuous-censored
    arma::mat eta;    //discrete-censored
    
    friend class LassoCox;

public:
    LassoCoxParams() {}
    LassoCoxParams(const arma::mat& gamma, const arma::mat& eta);
    LassoCoxParams(LassoCoxParams& parIn);
    LassoCoxParams(arma::vec& vec, int p, int ltot, int r);  //copy params from flattened vector

    /**
     * Copy all params into a single vector
     */
    arma::vec toMatrix1D();

    arma::mat& getGamma() { return gamma; }
    arma::mat& getEta() { return eta; }
    
    void setGamma(arma::mat& gamma) { this->gamma = gamma; }
    void setEta(arma::mat& eta) { this->eta = eta; }

    friend std::ostream& operator<<(std::ostream& os, LassoCoxParams& params);
    friend std::ostream& operator<<(std::ostream& os, LassoCoxParams&& params);
};

#endif /* LASSOCOXPARAMS_HPP_ */
