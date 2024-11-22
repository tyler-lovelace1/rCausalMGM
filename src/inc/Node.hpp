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
// #include "CensoredVariable.hpp"
#include <string>
// #include <memory>
#include <exception>
// #include <algorithm>
// #include <iterator>

class Node {
    
    std::string name = "??";
    DataType type = NONE;
    // std::unique_ptr<Variable> var{};
    Variable* var = nullptr;
    bool observed = true;

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

    // Node(CensoredVariable* var) :
    // 	name(var->getName()),
    // 	type(var->getType()),
    // 	var(var) {
    // 	// RcppThread::Rcout << "new censored variable " + name + ", type = " << type << "\n";
    // }
    
    ~Node() {
	// RcppThread::Rcout << "delete\n";
	delete var;
    }

    Node(const Node& other) {
	// RcppThread::Rcout << "copying (constructor) node " + other.name + " of type " << other.type << "\n";
	name = other.name;
	type = other.type;
	observed = other.observed;
	if (type==CONTINUOUS) {
	    var = (Variable*) new ContinuousVariable(other.var);
	} else if (type==DISCRETE) {
	    var = (Variable*) new DiscreteVariable(other.var);
	// } else if (type==CENSORED) {
	//     var = (Variable*) new CensoredVariable(other.var);
	} else if (other.isNull()) {
	    // RcppThread::Rcout << "copying null node\n";
	    var = nullptr;
	} else {
	    // RcppThread::Rcout << "copying (constructor) node " + name + "\n";
	    throw std::runtime_error("Invalid variable type detected for Node " + name);
	}
	// RcppThread::Rcout << "new\n";
    }

    Node& operator=(const Node& other) {
	// RcppThread::Rcout << "copying (assignment) node " + other.name + " of type " << other.type << "\n";
	if (this != &other) {
	    delete var;
	    name = other.name;
	    type = other.type;
	    observed = other.observed;
	    if (type==CONTINUOUS) {
		var = (Variable*) new ContinuousVariable(other.var);
	    } else if (type==DISCRETE) {
		var = (Variable*) new DiscreteVariable(other.var);
	    // } else if (type==CENSORED) {
	    // 	var = (Variable*) new CensoredVariable(other.var);
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
    
    Node(Node&& other) : Node() {
	// if (other.isNull()) RcppThread::Rcout << "moving null node\n";
	// else                RcppThread::Rcout << "moving " << other << "\n";
	// RcppThread::Rcout << "move\n";
	std::swap(name, other.name);
	std::swap(type, other.type);
	std::swap(observed, other.observed);
	std::swap(var, other.var);
    }

    Node& operator=(Node&& other) {
	if (this != &other) {
	    // if (other.isNull()) RcppThread::Rcout << "moving null node\n";
	    // else                RcppThread::Rcout << "moving " << other << "\n";
	    // delete var;
	    std::swap(name, other.name);
	    std::swap(type, other.type);
	    std::swap(observed, other.observed);
	    std::swap(var, other.var);
	}
	// RcppThread::Rcout << "move\n";
	return *this;
    }
    
    void setName(const std::string& name) { this->name = name; }

    void setObserved(bool observed) { this->observed = observed; }
  
    std::string getName() const { return name; }
  
    DataType getType() const { return type; }

    bool isNull() const { return var == nullptr; }
    
    bool isContinuous() const { return type==CONTINUOUS; }

    bool isDiscrete() const { return type==DISCRETE; }

    // bool isCensored() const { return type==CENSORED; }

    bool isObserved() const { return observed; }
  
    bool isMissingValue(const std::string& val) const {
	// RcppThread::Rcout << "checking missing value for node " + name + "\n";

        switch (type) {
	    
    	case CONTINUOUS:
	    return ((ContinuousVariable*)var)->isMissingValue(val);
    	case DISCRETE:
	    return ((DiscreteVariable*)var)->isMissingValue(val);
	// case CENSORED:
	//     return ((CensoredVariable*)var)->isMissingValue(val);
    	default:
	    throw std::runtime_error("Invalid variable type detected for Node " + name);
	      
    	}
    	return false;
    }
  
    bool checkValue(const std::string& val) const {
	// RcppThread::Rcout << "checking value for node " + name + "\n";

    	switch (type) {
	    
    	case CONTINUOUS:
	    return ((ContinuousVariable*)var)->checkValue(val);
    	case DISCRETE:
	    return ((DiscreteVariable*)var)->checkValue(val);
	// case CENSORED:
	//     return ((CensoredVariable*)var)->checkValue(val);
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

    // bool setCensor(arma::vec& values, arma::uvec& censor, arma::uvec& strata) {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->setCensor(values, censor, strata);
    // }
  
    // bool setCensor(arma::vec&& values, arma::uvec&& censor, arma::uvec&& strata) {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->setCensor(values, censor, strata);
    // }

    // std::vector<arma::uvec> getOrder() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getOrder();
    // }

    // arma::uvec getOrder(int strata) const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getOrder(strata);
    // }

    // arma::uvec getStrata() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getStrata();
    // }

    // arma::uvec getCensorVec() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getCensorVec();
    // }

    // std::vector<arma::uvec> getCensor() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getCensor();
    // }

    // arma::uvec getCensor(int strata) const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getCensor(strata);
    // }

    // std::vector<arma::uvec> getH() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getH();
    // }

    // arma::uvec getH(int strata) const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getH(strata);
    // }

    // std::vector<arma::uvec> getIndex() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getIndex();
    // }

    // arma::uvec getIndex(int strata) const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getIndex(strata);
    // }

    // arma::uvec getCC() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getCC();
    // }

    // arma::uword getCC(int i) const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getCC(i);
    // }

    // arma::uword getNumStrata() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getNumStrata();
    // }

    // arma::uword getNEvents() const {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getNEvents();
    // }

    // void setNeighbors(std::vector<std::string>& neighbors) {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    //     ((CensoredVariable*)var)->setNeighbors(neighbors);
    // }
    
    // std::vector<std::string> getNeighbors() {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getNeighbors();
    // }
    
    // void setWZ(arma::vec& WZ) {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	((CensoredVariable*)var)->setWZ(WZ);
    // }
    
    // arma::vec getWZ() {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	return ((CensoredVariable*)var)->getWZ();
    // }

    // void resetCoxModel(CoxIRLSRegression& coxRegression) {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	((CensoredVariable*)var)->resetCoxModel(coxRegression);
    // 	std::vector<Node> regressors(coxRegression.getVariables());
    // 	((CensoredVariable*)var)->coxResult = regress(*this, regressors);
    // }
    
    // arma::vec getWZ(std::vector<Node> xyList) {
    // 	if (type!=CENSORED)
    // 	    throw std::runtime_error("Node " + name + " is not censored");
    // 	std::vector<Node> regressors(((CensoredVariable*)var)->coxRegression.getVariables());
    // 	std::vector<Node> diff;
    // 	std::sort(xyList.begin(), xyList.end());
    // 	std::sort(regressors.begin(), regressors.end());
    // 	std::set_difference(regressors.begin(), regressors.end(),
    // 			    xyList.begin(), xyList.end(),
    // 			    std::inserter(diff, diff.begin()));
    // 	if (diff != regressors) {
    // 	    CoxRegressionResult
    // 	}
    // 	return ((CensoredVariable*)var)->getWZ();
    // }
    
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
