#ifndef INDEPENDENCETESTRANDOM_HPP_
#define INDEPENDENCETESTRANDOM_HPP_

#include "IndependenceTest.hpp"
#include <RcppArmadillo.h>

/**
 * A class that exists purely to test algorithms which require
 * an IndependenceTest. The result of the test is random,
 * but seeded with the variable names, so that it will give the
 * same result every time.
 * 
 * Not every function is implemented
 */ 

class IndependenceTestRandom : public IndependenceTest {

private:
    std::vector<Variable*> variables;

public:
    IndependenceTestRandom(std::vector<Variable*>& variables) { this->variables = variables; }

    /**
     * @return true about 25% of the time
     */
    bool isIndependent(Variable* x, Variable* y, std::vector<Variable*>& z);

    bool isDependent(Variable* x, Variable* y, std::vector<Variable*>& z) { return !isIndependent(x, y, z); }

    double getPValue() { return 0; }

    std::vector<Variable*> getVariables() { return variables; }

    Variable* getVariable(std::string name) { return NULL; }

    std::vector<std::string*> getVariableNames() { return {}; }

    bool determines(std::vector<Variable*>& z, Variable* y) { return false; }

    double getAlpha() { return 0; }

    void setAlpha(double alpha) {;}

    std::vector<DataSet*> getDataSets() { return {}; }

    int getSampleSize() { return 0; }
    
    double getScore() { return 0; }

};

#endif /* INDEPENDENCETESTRANDOM_HPP_ */