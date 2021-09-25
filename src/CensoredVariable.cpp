#include "CensoredVariable.hpp"

bool CensoredVariable::checkValue(const std::string& val) {
  try {
    double value = stod(val);
  } catch (std::exception& e) {
    return false;
  }
  return true;
}

bool CensoredVariable::setCensor(arma::vec& values, arma::uvec& censor) {
    this->censor = censor;
    n = values.n_elem;
    if (censor.size() != n)
	throw std::runtime_error("Number of values and number of censoring indicators must be equal.");
    
    order = arma::sort_index(values);
    
    std::vector<uint> H( {1} );
    
    std::vector<uint> CC;
    if (censor[0])
	CC.push_back(0);
    
    for (arma::uword i = 1; i < n; i++) {
	if (censor[i])
	    CC.push_back(i);
	if (values[order[i]] == values[order[i-1]])
	    H[H.size()-1]++;
	else
	    H.push_back(1);
    }
    unique = H.size();
    this->H = arma::uvec(H);
    this->CC = arma::uvec(CC);
    nEvents = arma::sum(censor);
    return true;
}

bool CensoredVariable::setCensor(arma::vec&& values, arma::uvec&& censor) {
    this->censor = censor;
    n = values.n_elem;
    if (censor.size() != n)
	throw std::runtime_error("Number of values and number of censoring indicators must be equal.");
    
    order = arma::sort_index(values);
    
    std::vector<uint> H( {1} );
    
    std::vector<uint> CC;
    if (censor[0])
	CC.push_back(0);
    
    for (arma::uword i = 1; i < n; i++) {
	if (censor[i])
	    CC.push_back(i);
	if (values[order[i]] == values[order[i-1]])
	    H[H.size()-1]++;
	else
	    H.push_back(1);
    }
    unique = H.size();
    this->H = arma::uvec(H);
    this->CC = arma::uvec(CC);
    nEvents = arma::sum(censor);
    return true;
}


// [[Rcpp::export]]
void test_censored(const std::string& name, arma::vec& values, arma::uvec& censor) {
  CensoredVariable cv = CensoredVariable(name);
  Rcpp::Rcout << std::boolalpha;
  Rcpp::Rcout << "Censored Variable Test" << std::endl;
  Rcpp::Rcout << "Name: " << cv.getName() << std::endl;
  Rcpp::Rcout << "Type: " << cv.getType() << std::endl;
  Rcpp::Rcout << "Missing Value: " << cv.getMissingValueMarker() << std::endl;
  Rcpp::Rcout << "checkValue Test: " << cv.checkValue(std::to_string(values[0])) << std::endl;
  Rcpp::Rcout << "setCensor Test: " << cv.setCensor(values, censor) << std::endl;
  Rcpp::Rcout << "Censored Variable Details:\n" << cv << std::endl;
  Rcpp::Rcout << "Values:\n" << values.t() << std::endl;
  Rcpp::Rcout << "Ordered Values:\n" << values(cv.order).t() << std::endl;
}

std::ostream& operator<<(std::ostream &os, CensoredVariable& cv) {
    os << "Name: " << cv.getName() << std::endl;
    os << "   N:\t" << cv.n << std::endl;
    os << "   Number unique:\t" << cv.unique << std::endl;
    os << "   Order:\n" << cv.order.t() << std::endl;
    // std::vector<uint> censor(cv.n);
    // for (int i = 0; i < cv.n; i++) {
    // 	if (cv.censor[i])
    // 	    censor.push_back(1);
    // 	else
    // 	    censor.push_back(0);
    // }
    os << "   Censor:\n" << cv.censor.t() << std::endl;
    os << "   H:\n" << cv.H.t() << std::endl;
    os << "   H sum:\t" << arma::sum(cv.H) << std::endl;
    os << "   CC:\n" << cv.CC.t() << std::endl;

    return os;
}    
