#include "ContinuousVariable.hpp"

bool ContinuousVariable::checkValue(const std::string& val) {
    try {
	double value = stod(val);
    } catch (std::exception& e) {
	return false;
    }
    return true;
}


// // no export // [[Rcpp::export]]
// void test_continuous(const std::string& name, const std::string& val) {
//   ContinuousVariable cv = ContinuousVariable(name);
//   Rcpp::Rcout << std::boolalpha;
//   Rcpp::Rcout << "Continuous Variable Test" << std::endl;
//   Rcpp::Rcout << "Name: " << cv.getName() << std::endl;
//   Rcpp::Rcout << "Type: " << cv.getType() << std::endl;
//   Rcpp::Rcout << "Missing Value: " << cv.getMissingValueMarker() << std::endl;
//   Rcpp::Rcout << "checkValue Test: " << cv.checkValue(val) << std::endl;
// }
