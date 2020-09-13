#ifndef INDEPENDENCETASK_HPP
#define INDEPENDENCETASK_HPP

#include "Variable.hpp"

/**
 * Represents a task testing if x _||_ y | z
 */ 
class IndependenceTask {

public:
    IndependenceTask(Variable* _x, Variable* _y, std::vector<Variable*>& _z) : x(_x), y(_y), z(_z) {}
    Variable* x;
    Variable* y;
    std::vector<Variable*> z;

};

#endif /* INDEPENDENCETASK_HPP */