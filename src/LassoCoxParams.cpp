#include "LassoCoxParams.hpp"

LassoCoxParams::LassoCoxParams(const arma::mat& gamma, const arma::mat& eta) {
    this->gamma = gamma;
    this->eta = eta;
}

// Copy constructor
LassoCoxParams::LassoCoxParams(LassoCoxParams& parIn) {
    this->gamma = arma::mat(parIn.gamma);
    this->eta = arma::mat(parIn.eta);
}

LassoCoxParams::LassoCoxParams(arma::vec& vec, int p, int ltot, int r) {
    std::vector<int> lens{p*r, ltot*r};
    std::vector<arma::uword> lenSums(lens.size());
    
    lenSums[0] = lens[0];
    for (int i = 1; i < lenSums.size(); i++) {
        lenSums[i] = lens[i] + lenSums[i-1];
    }

    if (vec.n_elem != lenSums[1])
        throw std::invalid_argument("Param vector dimension doesn't match: Found " + std::to_string(vec.n_elem) + " need " + std::to_string(lenSums[1]));

    gamma = arma::reshape(vec.subvec(0, lenSums[0]-1), p, r);
    eta = arma::reshape(vec.subvec(lenSums[0], lenSums[1]-1), ltot, r);
}

/**
 * Copy all params into a single vector
 */
arma::vec LassoCoxParams::toMatrix1D() {
    int p = gamma.n_rows;
    int ltot = eta.n_rows;
    int r = gamma.n_cols;
    std::vector<int> lens{p*r, ltot*r};
    std::vector<arma::uword> lenSums(lens.size());
    
    lenSums[0] = lens[0];
    for (int i = 1; i < lenSums.size(); i++) {
        lenSums[i] = lens[i] + lenSums[i-1];
    }

    arma::vec outVec(r * (p + ltot));
    
    outVec.subvec(0,          arma::size(gamma.as_col())) = gamma.as_col();
    outVec.subvec(lenSums[0], arma::size(eta.as_col())) = eta.as_col();

    return outVec;
}

std::ostream& operator<<(std::ostream& os, LassoCoxParams& params) {
    os << "gamma:\n"  << params.gamma.t() << "\n";
    os << "eta:\n"    << params.eta.t();
    return os;
}

std::ostream& operator<<(std::ostream& os, LassoCoxParams&& params) {
    os << "gamma:\n"  << params.gamma.t() << "\n";
    os << "eta:\n"    << params.eta.t();
    return os;
}
