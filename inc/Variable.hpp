/*
 * Base class for variable specifications for DataSet. These objects govern the
 * types of values which may be recorded in a Column of data and provide
 * information about the interpretation of these values. Variables of every
 * type must provide a marker which is recorded in a column of data for that
 * variable when the value is missing; this missing data marker should not be
 * used for other purposes.
 *
 * @author Willie Wheeler 7/99
 * @author Joseph Ramsey modifications 12/00
 * @author Tyler Lovelace conversion to C++ 1/20
 */

#ifndef VARIABLE_HPP_
#define VARIABLE_HPP_

#include "armaLapack.hpp"
#include <string>

enum DataType { CONTINUOUS, DISCRETE };

class Variable {
    
protected:
  std::string name = "??";
  DataType type;

public:
  Variable() {}
    
  virtual ~Variable() {}
  
  void setName(const std::string& name) { this->name = name; }
  
  std::string getName() { return name; }
  
  DataType getType() { return type; }

  bool isContinuous() { return type==CONTINUOUS; }

  bool isDiscrete() { return type==DISCRETE; }
  
  // template <typename T>
  // T getMissingValueMarker();
  
  virtual bool isMissingValue(const std::string& val) = 0;
  
  virtual bool checkValue(const std::string& val) = 0;
  
  virtual Variable* like(const std::string& name) = 0;
    
};

#endif /* VARIABLE_HPP_ */
