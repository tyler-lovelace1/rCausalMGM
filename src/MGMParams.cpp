#include "MGMParams.hpp"

MGMParams::MGMParams(const arma::mat& beta, const arma::vec& betad, const arma::mat& theta, const arma::mat& phi, const arma::vec& alpha1, const arma::vec& alpha2) {
    this->beta = beta;
    this->betad = betad;
    this->theta = theta;
    this->phi = phi;
    this->alpha1 = alpha1;
    this->alpha2 = alpha2;
}

// Copy constructor
MGMParams::MGMParams(MGMParams& parIn) {
    this->beta = arma::mat(parIn.beta);
    this->betad = arma::vec(parIn.betad);
    this->theta = arma::mat(parIn.theta);
    this->phi = arma::mat(parIn.phi);
    this->alpha1 = arma::vec(parIn.alpha1);
    this->alpha2 = arma::vec(parIn.alpha2);
}

MGMParams::MGMParams(arma::vec& vec, int p, int ltot) {
    std::vector<int> lens{p*p, p, p*ltot, ltot*ltot, p, ltot};
    std::vector<arma::uword> lenSums(lens.size());
    lenSums[0] = lens[0];

    for (int i = 1; i < lenSums.size(); i++) {
        lenSums[i] = lens[i] + lenSums[i-1];
    }

    if (vec.n_elem != lenSums[5])
        throw std::invalid_argument("Param vector dimension doesn't match: Found " + std::to_string(vec.n_elem) + " need " + std::to_string(lenSums[5]));

    beta = arma::reshape(vec.subvec(0, lenSums[0]-1), p, p);
    betad = arma::vec(vec.subvec(lenSums[0], lenSums[1]-1));
    theta = arma::reshape(vec.subvec(lenSums[1], lenSums[2]-1), ltot, p);
    phi = arma::reshape(vec.subvec(lenSums[2], lenSums[3]-1), ltot, ltot);
    alpha1 = arma::vec(vec.subvec(lenSums[3], lenSums[4]-1));
    alpha2 = arma::vec(vec.subvec(lenSums[4], lenSums[5]-1));
}

/**
 * Copy all params into a single vector
 */
arma::vec MGMParams::toMatrix1D() {
    int p = alpha1.n_elem;
    int ltot = alpha2.n_elem;
    std::vector<int> lens{p*p, p, p*ltot, ltot*ltot, p, ltot};
    std::vector<arma::uword> lenSums(lens.size());
    lenSums[0] = lens[0];
    for (int i = 1; i < lenSums.size(); i++) {
        lenSums[i] = lens[i] + lenSums[i-1];
    }

    arma::vec outVec(p*p + p + p*ltot + ltot*ltot + p + ltot);
    
    outVec.subvec(0,          arma::size(beta.as_col())) = beta.as_col();
    outVec.subvec(lenSums[0], arma::size(betad)) = betad;
    outVec.subvec(lenSums[1], arma::size(theta.as_col())) = theta.as_col();
    outVec.subvec(lenSums[2], arma::size(phi.as_col())) = phi.as_col();
    outVec.subvec(lenSums[3], arma::size(alpha1)) = alpha1;
    outVec.subvec(lenSums[4], arma::size(alpha2)) = alpha2;
    
    return outVec;
}

std::ostream& operator<<(std::ostream& os, MGMParams& params) {
    os << "alpha1:\n" << params.alpha1.t() << "\n";
    os << "alpha2:\n" << params.alpha2.t() << "\n";
    os << "beta:\n"   << params.beta       << "\n";
    os << "betad:\n"  << params.betad.t()  << "\n";
    os << "theta:\n"  << params.theta      << "\n";
    os << "phi:\n"    << params.phi;
    return os;
}

