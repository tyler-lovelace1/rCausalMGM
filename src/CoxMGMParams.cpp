#include "CoxMGMParams.hpp"

CoxMGMParams::CoxMGMParams(const arma::mat& beta, const arma::vec& betad,
			   const arma::mat& theta, const arma::mat& phi,
			   const arma::mat& gamma, const arma::mat& psi,
			   const arma::vec& alpha1, const arma::vec& alpha2,
			   const arma::vec& alpha3) {
    this->beta = beta;
    this->betad = betad;
    this->theta = theta;
    this->phi = phi;
    this->gamma = gamma;
    this->psi = psi;
    this->alpha1 = alpha1;
    this->alpha2 = alpha2;
    this->alpha3 = alpha3;
}

// // Copy constructor
// CoxMGMParams::CoxMGMParams(CoxMGMParams& parIn) {
//     this->beta = arma::mat(parIn.beta);
//     this->betad = arma::vec(parIn.betad);
//     this->theta = arma::mat(parIn.theta);
//     this->phi = arma::mat(parIn.phi);
//     this->gamma = arma::mat(parIn.gamma);
//     this->psi = arma::mat(parIn.psi);
//     this->alpha1 = arma::vec(parIn.alpha1);
//     this->alpha2 = arma::vec(parIn.alpha2);
//     this->alpha3 = arma::vec(parIn.alpha3);
// }

CoxMGMParams::CoxMGMParams(arma::vec& vec, int p, int ltot, int r) {
    std::vector<int> lens{p*p, p, p*ltot, ltot*ltot, p*r, ltot*r, p, ltot, r};
    std::vector<arma::uword> lenSums(lens.size());
    lenSums[0] = lens[0];

    for (int i = 1; i < lenSums.size(); i++) {
        lenSums[i] = lens[i] + lenSums[i-1];
    }

    if (vec.n_elem != lenSums[8])
        throw std::invalid_argument("Param vector dimension doesn't match: Found " + std::to_string(vec.n_elem) + " need " + std::to_string(lenSums[8]));

    beta = arma::reshape(vec.subvec(0, lenSums[0]-1), p, p);
    betad = arma::vec(vec.subvec(lenSums[0], lenSums[1]-1));
    theta = arma::reshape(vec.subvec(lenSums[1], lenSums[2]-1), ltot, p);
    phi = arma::reshape(vec.subvec(lenSums[2], lenSums[3]-1), ltot, ltot);
    gamma = arma::reshape(vec.subvec(lenSums[3], lenSums[4]-1), p, r);
    psi = arma::reshape(vec.subvec(lenSums[4], lenSums[5]-1), ltot, r);
    alpha1 = arma::vec(vec.subvec(lenSums[5], lenSums[6]-1));
    alpha2 = arma::vec(vec.subvec(lenSums[6], lenSums[7]-1));
    alpha3 = arma::vec(vec.subvec(lenSums[7], lenSums[8]-1));
}

/**
 * Copy all params into a single vector
 */
arma::vec CoxMGMParams::toMatrix1D() {
    int p = alpha1.n_elem;
    int ltot = alpha2.n_elem;
    int r = alpha3.n_elem;
    std::vector<int> lens{p*p, p, p*ltot, ltot*ltot, p*r, ltot*r, p, ltot, r};
    std::vector<arma::uword> lenSums(lens.size());
    lenSums[0] = lens[0];
    for (int i = 1; i < lenSums.size(); i++) {
        lenSums[i] = lens[i] + lenSums[i-1];
    }

    arma::vec outVec(p*p + p + p*ltot + ltot*ltot + p*r + ltot*r + p + ltot + r);
    
    outVec.subvec(0,          arma::size(beta.as_col())) = beta.as_col();
    outVec.subvec(lenSums[0], arma::size(betad)) = betad;
    outVec.subvec(lenSums[1], arma::size(theta.as_col())) = theta.as_col();
    outVec.subvec(lenSums[2], arma::size(phi.as_col())) = phi.as_col();
    outVec.subvec(lenSums[3], arma::size(gamma.as_col())) = gamma.as_col();
    outVec.subvec(lenSums[4], arma::size(psi.as_col())) = psi.as_col();
    outVec.subvec(lenSums[5], arma::size(alpha1)) = alpha1;
    outVec.subvec(lenSums[6], arma::size(alpha2)) = alpha2;
    outVec.subvec(lenSums[7], arma::size(alpha3)) = alpha3;
    
    return outVec;
}

std::ostream& operator<<(std::ostream& os, CoxMGMParams& params) {
    os << "alpha1:\n" << params.alpha1.t() << "\n";
    os << "alpha2:\n" << params.alpha2.t() << "\n";
    os << "alpha3:\n" << params.alpha3.t() << "\n";
    os << "beta:\n"   << params.beta       << "\n";
    os << "betad:\n"  << params.betad.t()  << "\n";
    os << "theta:\n"  << params.theta      << "\n";
    os << "phi:\n"    << params.phi        << "\n";
    os << "gamma:\n"  << params.gamma      << "\n";
    os << "psi:\n"    << params.psi;
    return os;
}
