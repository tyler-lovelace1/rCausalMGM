/*
 * Represents a real-valued variable. The values are doubles, and the default
 * missing value marker for is numeric_limits<double>::quiet_NaN().
 *
 * @author Willie Wheeler 7/99
 * @author Joseph Ramsey modifications 12/00
 * @author Tyler Lovelace conversion to C++ 1/20
*/

#ifndef CONTINUOUSVARIABLE_HPP_
#define CONTINUOUSVARIABLE_HPP_

#include "Variable.hpp"
#include <cmath>
#include <limits>
#include <exception>
#include <iostream>

class ContinuousVariable : public Variable {
private:
    const double MISSING_VALUE = std::numeric_limits<double>::quiet_NaN();
    const std::string MISSING_VALUE_STRING = "NaN";

public:
    ContinuousVariable(const std::string& name) {
	this->name = name;
	this->type = CONTINUOUS;
    }

    ContinuousVariable(Variable* var) {
	if (var->getType() == CONTINUOUS) {
	    this->name = var->getName();
	    this->type = CONTINUOUS;
	} else {
	    throw std::runtime_error("Trying to construct continuous variable from a variable of a different type");
	}
    }

    ContinuousVariable(const ContinuousVariable& var) = default;

    ContinuousVariable& operator=(const ContinuousVariable& var) = default;
    
    ContinuousVariable(ContinuousVariable&& var) = default;

    ContinuousVariable& operator=(ContinuousVariable&& var) = default;     
  
    ~ContinuousVariable() = default;

    double getMissingValueMarker() { return this->MISSING_VALUE; }

    bool isMissingValue(const double& val) { return std::isnan(val); }

    bool isMissingValue(const std::string& val) { return val==MISSING_VALUE_STRING; }

    // bool checkValue(const double& val) { return true; }

    bool checkValue(const std::string& val);
  
    // ContinuousVariable* like(const std::string& name) { return new ContinuousVariable(name); }

    // friend void test_continuous(const std::string& name, const std::string& val);
};

#endif /* CONTINUOUSVARIABLE_HPP_ */
