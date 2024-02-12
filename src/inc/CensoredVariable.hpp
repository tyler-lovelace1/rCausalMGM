/*
 * Represents a right censored variable with a real value and a censoring indicator. 
 * The values are doubles, and the default missing value marker for is 
 * numeric_limits<double>::quiet_NaN(). The censoring indicator is boolean value, 
 * false when the value is censored and true when uncensored.
 *
 * @author Tyler Lovelace conversion to C++ 7/21
 */

#ifndef CENSOREDVARIABLE_HPP_
#define CENSOREDVARIABLE_HPP_

// class Node;

#include "Variable.hpp"
#include <cmath>
#include <limits>
#include <exception>
#include <iostream>

class CensoredVariable : public Variable {
private:
    const double MISSING_VALUE = std::numeric_limits<double>::quiet_NaN();
    const std::string MISSING_VALUE_STRING = "NaN";
    uint n;
    uint unique;
    uint nEvents;
    uint numStrata;
    arma::uvec strata;
    arma::uvec censorVec;
    std::vector<arma::uvec> order;
    std::vector<arma::uvec> censor;
    std::vector<arma::uvec> H;
    std::vector<arma::uvec> index;
    // arma::uvec CC;
    arma::vec WZ;
    std::vector<std::string> neighbors;

public:
    CensoredVariable(const std::string& name) {
	this->name = name;
	this->type = CENSORED;
    }

    CensoredVariable(Variable* var) {
	if (var->getType() == CENSORED) {
	    name = ((CensoredVariable*) var)->name;
	    type = ((CensoredVariable*) var)->type;
	    n = ((CensoredVariable*) var)->n;
	    unique = ((CensoredVariable*) var)->unique;
	    nEvents = ((CensoredVariable*) var)->nEvents;
	    numStrata = ((CensoredVariable*) var)->numStrata;
	    strata = ((CensoredVariable*) var)->strata;
	    censorVec = ((CensoredVariable*) var)->censorVec;
	    order = ((CensoredVariable*) var)->order;
	    censor = ((CensoredVariable*) var)->censor;
	    H = ((CensoredVariable*) var)->H;
	    index = ((CensoredVariable*) var)->index;
	    // CC = arma::uvec(((CensoredVariable*) var)->CC);
	    WZ = ((CensoredVariable*) var)->WZ;
	    neighbors = ((CensoredVariable*) var)->neighbors;
	} else {
	    throw std::runtime_error("Trying to construct censored variable from a variable of a different type");
	}
    }

    CensoredVariable(CensoredVariable* var) {
	name = var->name;
	type = var->type;
	n = var->n;
	unique = var->unique;
	nEvents = var->nEvents;
	numStrata = var->numStrata;
	strata = var->strata;
	censorVec = var->censorVec;
	order = var->order;
	censor = var->censor;
	H = var->H;
	index = var->index;
	// CC = arma::uvec(var->CC);
	WZ = var->WZ;
	neighbors = var->neighbors;
    }

    // CensoredVariable(const CensoredVariable& var) {
    // 	name = var.name;
    // 	type = var.type;
    // 	n = var.n;
    // 	unique = var.unique;
    // 	nEvents = var.nEvents;
    // 	numStrata = var.numStrata;
    // 	strata = var.strata;
    // 	censorVec = var.censorVec;
    // 	order = var.order;
    // 	censor = var.censor;
    // 	H = var.H;
    // 	// CC = var.CC;
    // 	WZ = var.WZ;
    // 	neighbors = var.neighbors;
    // }

    CensoredVariable(const CensoredVariable& var) = default;

    CensoredVariable& operator=(const CensoredVariable& var) = default;
    
    CensoredVariable(CensoredVariable&& var) = default;

    CensoredVariable& operator=(CensoredVariable&& var) = default;     
  
    ~CensoredVariable() = default;
  
    double getMissingValueMarker() { return this->MISSING_VALUE; }

    bool setCensor(arma::vec& values, arma::uvec& censor, arma::uvec& strata);

    bool setCensor(arma::vec&& values, arma::uvec&& censor, arma::uvec&& strata);
    
    // bool setCensor(arma::vec& values, arma::uvec& censor);
    
    // bool setCensor(arma::vec&& values, arma::uvec&& censor);

    std::vector<arma::uvec> getOrder() { return this->order; }

    arma::uvec getOrder(int strata) { return this->order[strata]; }

    arma::uvec getCensorVec() { return this->censorVec; }

    arma::uvec getStrata() { return this->strata; }
    
    std::vector<arma::uvec> getCensor() { return this->censor; }

    arma::uvec getCensor(int strata) { return this->censor[strata]; }

    std::vector<arma::uvec> getH() { return this->H; }

    arma::uvec getH(int strata) { return this->H[strata]; }

    std::vector<arma::uvec> getIndex() { return this->index; }

    arma::uvec getIndex(int strata) { return this->index[strata]; }

    arma::uword getNumStrata() { return numStrata; }

    // arma::uvec getCC() { return this->CC; }

    // arma::uword getCC(int i) { return this->CC[i]; }

    arma::uword getNEvents() { return nEvents; }

    // void resetCoxModel(CoxIRLSRegression& coxRegression) { this->coxRegression = coxRegression; }

    void setNeighbors(std::vector<std::string>& neighbors) { this->neighbors = neighbors; }
    
    std::vector<std::string> getNeighbors() { return neighbors; }
    
    void setWZ(arma::vec& WZ) { this->WZ = WZ; }
    
    arma::vec getWZ() { return WZ; }

    bool isMissingValue(const double& val) { return std::isnan(val); }

    bool isMissingValue(const std::string& val) { return val==MISSING_VALUE_STRING; }

    // bool checkValue(const double& val) { return true; }

    bool checkValue(const std::string& val);
  
    CensoredVariable* like(const std::string& name) { return new CensoredVariable(name); }

    friend void test_censored(const std::string& name, arma::vec& values, arma::uvec& censor, arma::uvec& strata);
    friend std::ostream& operator<<(std::ostream& os, CensoredVariable& cv);
};

#endif /* CENSOREDVARIABLE_HPP_ */
