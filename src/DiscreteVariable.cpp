#include "DiscreteVariable.hpp"

void DiscreteVariable::setCategories(const int& numCategories) {
  for (int i = 0; i < numCategories; i++)
    this->categories.push_back("category" + std::to_string(i));
}

int DiscreteVariable::getIndex(const std::string& category) {
  for (int i = 0; i < categories.size(); i++) {
    if (categories.at(i)==category) {
      return i;
    }
  }
  throw std::runtime_error("category is not in categories for " + getName());
}

bool DiscreteVariable::checkValue(const std::string& val) {
  for (int i = 0; i < categories.size(); i++) {
    if (categories.at(i)==val)
      return true;
  }
  return false;
}


// // [[Rcpp::export]]
// void test_discrete(const std::string& name, const int& numCats, const std::string& val) {
//   DiscreteVariable dv = DiscreteVariable(name, numCats);
//   Rcpp::Rcout << std::boolalpha;
//   Rcpp::Rcout << "Discrete Variable Test" << std::endl;
//   Rcpp::Rcout << "Name: " << dv.getName() << std::endl;
//   Rcpp::Rcout << "Type: " << dv.getType() << std::endl;
//   Rcpp::Rcout << "Missing Value: " << dv.getMissingValueMarker() << std::endl;
//   Rcpp::Rcout << "checkValue Test: " << dv.checkValue(val) << std::endl;
// }
