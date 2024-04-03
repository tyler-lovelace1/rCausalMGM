#include "CensoredVariable.hpp"

bool CensoredVariable::checkValue(const std::string& val) {
  try {
    double value = stod(val);
  } catch (std::exception& e) {
    return false;
  }
  return true;
}

bool CensoredVariable::setCensor(arma::vec& values, arma::uvec& censor, arma::uvec& strata) {
    this->censorVec = censor;
    this->strata = strata;
    n = values.n_elem;
    if (censor.size() != n || strata.size() != n)
	throw std::runtime_error("Number of values, number of censoring indicators, and number of strata labels must be equal.");

    arma::uvec uniqStrata = arma::unique(strata);
    numStrata = uniqStrata.n_elem;

    this->order.clear();
    this->censor.clear();
    this->index.clear();
    this->H.clear();

    // Rcpp::Rcout << name << std::endl;
    // Rcpp::Rcout << "Unique strata: " << uniqStrata.t(); // << std::endl;
    // Rcpp::Rcout << "Num strata: " << numStrata << std::endl;

    for (uint strat = 0; strat < numStrata; strat++) {

	// Rcpp::Rcout << "  Strata " << strat << ":" << std::endl;

	arma::uvec index = arma::find(strata == strat);
	this->index.push_back(index);

	// Rcpp::Rcout << "  Index: " << index.t();
	
	this->censor.push_back(this->censorVec(index));

	// Rcpp::Rcout << "  Censor: " << this->censor[strat].t();
    
	arma::uvec order = arma::sort_index(values(index));
	this->order.push_back(order);

	// Rcpp::Rcout << "  Order: " << this->order[strat].t();
    
	std::vector<uint> H( {1} );
    
	// std::vector<uint> CC;
	// if (censor[0])
	// 	CC.push_back(0);
    
	for (arma::uword i = 1; i < index.n_elem; i++) {
	    // if (censor[i])
	    //     CC.push_back(i);
	    if (values[index[order[i]]] == values[index[order[i-1]]])
		H[H.size()-1]++;
	    else
		H.push_back(1);
	}

	this->H.push_back(arma::uvec(H));

	// Rcpp::Rcout << "  H: " << this->H[strat].t();
    }

    arma::vec uniqVals = arma::unique(values);
    unique = uniqVals.n_elem;
    // this->CC = arma::uvec(CC);
    nEvents = arma::sum(censorVec);
    return true;
}

bool CensoredVariable::setCensor(arma::vec&& values, arma::uvec&& censor, arma::uvec&& strata) {
    this->censorVec = censor;
    this->strata = strata;
    n = values.n_elem;
    if (censor.size() != n || strata.size() != n)
	throw std::runtime_error("Number of values, number of censoring indicators, and number of strata labels must be equal.");

    arma::uvec uniqStrata = arma::unique(strata);
    numStrata = uniqStrata.n_elem;

    this->order.clear();
    this->censor.clear();
    this->index.clear();
    this->H.clear();

    // Rcpp::Rcout << name << std::endl;
    // Rcpp::Rcout << "Unique strata: " << uniqStrata.t(); // << std::endl;
    // Rcpp::Rcout << "Num strata: " << numStrata << std::endl;

    for (uint strat = 0; strat < numStrata; strat++) {

	// Rcpp::Rcout << "  Strata " << strat << ":" << std::endl;

	arma::uvec index = arma::find(strata == strat);
	this->index.push_back(index);

	// Rcpp::Rcout << "  Index: " << index.t();
	
	this->censor.push_back(this->censorVec(index));

	// Rcpp::Rcout << "  Censor: " << this->censor[strat].t();
    
	arma::uvec order = arma::sort_index(values(index));
	this->order.push_back(order);

	// Rcpp::Rcout << "  Order: " << this->order[strat].t();
    
	std::vector<uint> H( {1} );
    
	// std::vector<uint> CC;
	// if (censor[0])
	// 	CC.push_back(0);
    
	for (arma::uword i = 1; i < index.n_elem; i++) {
	    // if (censor[i])
	    //     CC.push_back(i);
	    if (values[index[order[i]]] == values[index[order[i-1]]])
		H[H.size()-1]++;
	    else
		H.push_back(1);
	}

	this->H.push_back(arma::uvec(H));

	// Rcpp::Rcout << "  H: " << this->H[strat].t();
    }

    arma::vec uniqVals = arma::unique(values);
    unique = uniqVals.n_elem;
    // this->CC = arma::uvec(CC);
    nEvents = arma::sum(censorVec);
    return true;
}


// bool CensoredVariable::setCensor(arma::vec& values, arma::uvec& censor) {
//     this->censor = censor;
//     n = values.n_elem;
//     if (censor.size() != n)
// 	throw std::runtime_error("Number of values and number of censoring indicators must be equal.");
    
//     order = arma::sort_index(values);
    
//     std::vector<uint> H( {1} );
    
//     // std::vector<uint> CC;
//     // if (censor[0])
//     // 	CC.push_back(0);
    
//     for (arma::uword i = 1; i < n; i++) {
// 	// if (censor[i])
// 	//     CC.push_back(i);
// 	if (values[order[i]] == values[order[i-1]])
// 	    H[H.size()-1]++;
// 	else
// 	    H.push_back(1);
//     }
//     unique = H.size();
//     this->H = arma::uvec(H);
//     // this->CC = arma::uvec(CC);
//     nEvents = arma::sum(censor);
//     return true;
// }

// bool CensoredVariable::setCensor(arma::vec&& values, arma::uvec&& censor) {
//     this->censor = censor;
//     n = values.n_elem;
//     if (censor.size() != n)
// 	throw std::runtime_error("Number of values and number of censoring indicators must be equal.");
    
//     order = arma::sort_index(values);
    
//     std::vector<uint> H( {1} );
    
//     // std::vector<uint> CC;
//     // if (censor[0])
//     // 	CC.push_back(0);
    
//     for (arma::uword i = 1; i < n; i++) {
// 	// if (censor[i])
// 	//     CC.push_back(i);
// 	if (values[order[i]] == values[order[i-1]])
// 	    H[H.size()-1]++;
// 	else
// 	    H.push_back(1);
//     }
//     unique = H.size();
//     this->H = arma::uvec(H);
//     // this->CC = arma::uvec(CC);
//     nEvents = arma::sum(censor);
//     return true;
// }


// no export // [[Rcpp::export]]
void test_censored(const std::string& name, arma::vec& values, arma::uvec& censor, arma::uvec& strata) {
  CensoredVariable cv = CensoredVariable(name);
  Rcpp::Rcout << std::boolalpha;
  Rcpp::Rcout << "Censored Variable Test" << std::endl;
  Rcpp::Rcout << "Name: " << cv.getName() << std::endl;
  Rcpp::Rcout << "Type: " << cv.getType() << std::endl;
  Rcpp::Rcout << "Missing Value: " << cv.getMissingValueMarker() << std::endl;
  Rcpp::Rcout << "checkValue Test: " << cv.checkValue(std::to_string(values[0])) << std::endl;
  Rcpp::Rcout << "setCensor Test: " << cv.setCensor(values, censor, strata) << std::endl;
  Rcpp::Rcout << "Censored Variable Details:\n" << cv << std::endl;
  Rcpp::Rcout << "Values:\n" << values.t() << std::endl;
  Rcpp::Rcout << "Ordered Values:\n";
  for (uint strat = 0; strat < cv.getNumStrata(); strat++) {
      Rcpp::Rcout << "   Strata " << strat << ":\n";
      Rcpp::Rcout << "      Ordered Values:\n" << values(cv.getIndex(strat)(cv.getOrder(strat))).t() << std::endl;
  }
}

std::ostream& operator<<(std::ostream &os, CensoredVariable& cv) {
    os << "Name: " << cv.getName() << std::endl;
    os << "   N:                 " << cv.n << std::endl;
    os << "   Number unique:     " << cv.unique << std::endl;
    os << "   Number of strata:  " << cv.numStrata << std::endl;
    os << "   Strata:            " << cv.strata.t() << std::endl;
    for (uint strat = 0; strat < cv.numStrata; strat++) {
	os << "   Strata " << strat << ":\n";
	os << "      Order:          " << cv.order[strat].t() << std::endl;
	os << "      Censor:         " << cv.censor[strat].t() << std::endl;
        os << "      H:              " << cv.H[strat].t() << std::endl;
        os << "      H sum:          " << arma::sum(cv.H[strat]) << std::endl;
    }

    return os;
}    
