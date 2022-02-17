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
#include "RcppThread.h"
#include "Variable.hpp"
#include "ContinuousVariable.hpp"
#include "DiscreteVariable.hpp"
#include <string>
// #include <memory>
#include <exception>

class Node {
    
    std::string name = "??";
    DataType type = NONE;
    // std::unique_ptr<Variable> var{};
    Variable* var = nullptr;

public:
    Node() {
	// RcppThread::Rcout << "new\n";
    }

    Node(Variable* var) :
	name(var->getName()),
	type(var->getType()),
	var(var) {
	// RcppThread::Rcout << "new\n";
    }
    
    Node(ContinuousVariable* var) :
	name(var->getName()),
	type(var->getType()),
	var(var) {
	// RcppThread::Rcout << "new\n";
    }

    Node(DiscreteVariable* var) :
	name(var->getName()),
	type(var->getType()),
	var(var) {
    	// RcppThread::Rcout << "new\n";
    }
    
    ~Node() {
	// RcppThread::Rcout << "delete\n";
	delete var;
    }

    Node(const Node& other) {
	name = other.name;
	type = other.type;
	if (type==CONTINUOUS) {
	    var = (Variable*) new ContinuousVariable(other.var);
	} else if (type==DISCRETE) {
	    var = (Variable*) new DiscreteVariable(other.var);
	} else if (other.isNull()) {
	    // RcppThread::Rcout << "copying null node\n";
	    var = nullptr;
	} else {
	    throw std::runtime_error("Invalid variable type detected for Node " + name);
	}
	// RcppThread::Rcout << "new\n";
    }

    Node& operator=(const Node& other) {
	if (this != &other) {
	    delete var;
	    name = other.name;
	    type = other.type;
	    if (type==CONTINUOUS) {
		var = (Variable*) new ContinuousVariable(other.var);
	    } else if (type==DISCRETE) {
		var = (Variable*) new DiscreteVariable(other.var);
	    } else if (other.isNull()) {
		// RcppThread::Rcout << "copying null node\n";
		var = nullptr;
	    }else {
		throw std::runtime_error("Invalid variable type detected for Node " + name);
	    }
	}
	// RcppThread::Rcout << "new\n";
	return *this;
    }
    
    Node(Node&& other) : name(std::move(other.name)),
			 type(std::exchange(other.type, NONE)),
			 var(std::exchange(other.var, nullptr)) {
	// if (other.isNull()) RcppThread::Rcout << "moving null node\n";
	// else                RcppThread::Rcout << "moving " << *this << "\n";
	// RcppThread::Rcout << "move\n";
    }

    Node& operator=(Node&& other) {
	if (this != &other) {
	    // if (other.isNull()) RcppThread::Rcout << "moving null node\n";
	    // else                RcppThread::Rcout << "moving " << other << "\n";
	    delete var;
	    name = std::move(other.name);
	    type = std::exchange(other.type, NONE);
	    var = std::exchange(other.var, nullptr);
	}
	// RcppThread::Rcout << "move\n";
	return *this;
    }
    
    void setName(const std::string& name) { this->name = name; }
  
    std::string getName() const { return name; }
  
    DataType getType() const { return type; }

    bool isNull() const { return var == nullptr; }
    
    bool isContinuous() const { return type==CONTINUOUS; }

    bool isDiscrete() const { return type==DISCRETE; }
  
    bool isMissingValue(const std::string& val) const {
        switch (type) {
	    
    	case CONTINUOUS:
	    return ((ContinuousVariable*)var)->isMissingValue(val);
    	case DISCRETE:
	    return ((DiscreteVariable*)var)->isMissingValue(val);
    	default:
	    throw std::runtime_error("Invalid variable type detected for Node " + name);
	      
    	}
    	return false;
    }
  
    bool checkValue(const std::string& val) const {
    	switch (type) {
	    
    	case CONTINUOUS:
	    return ((ContinuousVariable*)var)->checkValue(val);
    	case DISCRETE:
	    return ((DiscreteVariable*)var)->checkValue(val);
    	default:
	    throw std::runtime_error("Invalid variable type detected for Node " + name);
	      
    	}
    	return false;
    }

    int getIndex(const std::string& category) const {
    	if (type!=DISCRETE)
    	    throw std::runtime_error("Node " + name + " is not discrete");
    	return ((DiscreteVariable*)var)->getIndex(category);
    }

    int getNumCategories() const {
    	if (type!=DISCRETE)
    	    throw std::runtime_error("Node " + name + " is not discrete");
    	return ((DiscreteVariable*)var)->getNumCategories();
    }

    std::vector<std::string> getCategories() const {
    	if (type!=DISCRETE)
    	    throw std::runtime_error("Node " + name + " is not discrete");
    	return ((DiscreteVariable*)var)->getCategories();
    }
    
    friend struct std::hash<Node>;
    friend std::size_t hash_value(const Node& n);
    
    friend std::ostream& operator<<(std::ostream& os, const Node& n);
    friend bool operator==(const Node& n1, const Node& n2);
    friend bool operator!=(const Node& n1, const Node& n2);
    friend bool operator> (const Node& n1, const Node& n2);
    friend bool operator<= (const Node& n1, const Node& n2);
    friend bool operator< (const Node& n1, const Node& n2);
    friend bool operator>= (const Node& n1, const Node& n2);

};

// bool operator==(const Node& n1, const Node& n2) {
//     return n1.name == n2.name && n1.type == n2.type;
// }

// bool operator!=(const Node& n1, const Node& n2) {
//     return !(n1 == n2);
// }

// bool operator< (const Node& n1, const Node& n2) {
//     return n1.name < n2.name;
// }

// bool operator<=(const Node& n1, const Node& n2) {
//     return n1 < n2 || n1 == n2;
// }

// bool operator> (const Node& n1, const Node& n2) {
//     return n1.name > n2.name;
// }

// bool operator>=(const Node& n1, const Node& n2) {
//     return n1 > n2 || n1 == n2;
// }

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
    std::size_t operator()(const Node& n) const {
	std::size_t res = 17;
	res = res * 31 + std::hash<std::string>()(n.name);
	// res = res * 43 + std::hash<DataType>()(n.type);
	return res;
    }
};

#endif /* NODE_HPP_ */
