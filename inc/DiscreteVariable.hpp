/*
 * Represents a real-valued variable. The values are doubles, and the default
 * missing value marker for is numeric_limits<double>::quiet_NaN().
 *
 * @author Willie Wheeler 7/99
 * @author Joseph Ramsey modifications 12/00
`* @author Tyler Lovelace conversion to C++ 1/20
 */

#ifndef DISCRETEVARIABLE_HPP_
#define DISCRETEVARIABLE_HPP_

#include "Variable.hpp"
#include <vector>
#include <exception>
#include <iostream>

class DiscreteVariable : public Variable {
private:
  const int MISSING_VALUE = -99;
  const std::string MISSING_VALUE_STRING = "*";
  std::vector<std::string> categories;

  void setCategories(const int& numCategories);

public:
  DiscreteVariable(const std::string& name) {
    this->name = name;
    this->type = DISCRETE;
  }

  DiscreteVariable(const std::string& name, const int& numCategories) {
    this->name = name;
    this->type = DISCRETE;
    setCategories(numCategories);
  }

  DiscreteVariable(const std::string& name, const std::vector<std::string>& categories) {
    this->name = name;
    this->type = DISCRETE;
    for (int i = 0; i < categories.size(); i++)
      this->categories.push_back(categories.at(i));
  }

  DiscreteVariable(const DiscreteVariable& var) {
    this->name = var.name;
    this->type = var.type;
    for (int i = 0; i < var.categories.size(); i++)
      this->categories.push_back(var.categories.at(i));
  }

  int getMissingValueMarker() { return MISSING_VALUE; }

  int getIndex(const std::string& category);

  int getNumCategories() { return categories.size(); }

  std::vector<std::string> getCategories() { return categories; }

  bool isMissingValue(const int& val) { return val==MISSING_VALUE; }

  bool isMissingValue(const std::string& val) { return val==MISSING_VALUE_STRING || val=="-99"; }

  // bool checkValue(const int& category) { return (category >= 0) && (category < getNumCategories()); }

  bool checkValue(const std::string& val);

  DiscreteVariable* like(const std::string& name) { return new DiscreteVariable(name); }

  // friend void test_discrete(const std::string& name, const int& numCats, const std::string& val);

};

#endif /* DISCRETEVARIABLE_HPP_ */
