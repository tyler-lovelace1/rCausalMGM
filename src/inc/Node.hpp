/*
 * Base class for Node specifications for DataSet. These objects govern the
 * types of values which may be recorded in a Column of data and provide
 * information about the interpretation of these values. Nodes of every
 * type must provide a marker which is recorded in a column of data for that
 * Node when the value is missing; this missing data marker should not be
 * used for other purposes.
 *
 * @author Willie Wheeler 7/99
 * @author Joseph Ramsey modifications 12/00
 * @author Tyler Lovelace conversion to C++ 1/20
 */

#ifndef NODE_HPP_
#define NODE_HPP_

#include "armaLapack.hpp"
#include <string>

enum DataType { CONTINUOUS, DISCRETE };

class Node {
    
    std::string name = "??";
    DataType type;
    Variable* var;

public:
    Node() {}
    
    ~Node() { delete var; }
  
    void setName(const std::string& name) { this->name = name; }
  
    std::string getName() { return name; }
  
    DataType getType() { return type; }

    bool isContinuous() { return type==CONTINUOUS; }

    bool isDiscrete() { return type==DISCRETE; }
  
    // template <typename T>
    // T getMissingValueMarker();
  
    bool isMissingValue(const std::string& val) {
	return false;
    }
  
    bool checkValue(const std::string& val) {
	return true;
    }
  
    // virtual Node* like(const std::string& name) = 0;

    friend struct std::hash<Node>;
    
};

// template<> struct std::hash<DataType> {
// public:
//     std::size_t operator()(const DataType t) const {
// 	std::size_t res = 17;
// 	res = res * 31 + std::hash<std::string>()(k.getName());
// 	res = res * 43 + std::hash<DataType>()(n.getType());
// 	return res;
//     }
// };


template<> struct std::hash<Node> {
public:
    std::size_t operator()(const Node n) const {
	std::size_t res = 17;
	res = res * 31 + std::hash<std::string>()(n.getName());
	res = res * 43 + std::hash<DataType>()(n.getType());
	return res;
    }
};


#endif /* NODE_HPP_ */
