#ifndef TRIPLE_HPP_
#define TRIPLE_HPP_

#include "Variable.hpp"

class Triple {
private:
    Variable* x;
    Variable* y;
    Variable* z;

public:

    Triple(Variable* x, Variable* y, Variable* z);

    Variable* getX() { return x; }
    Variable* getY() { return y; }
    Variable* getZ() { return z; }

    friend std::ostream& operator<<(std::ostream& os, Triple& triple);
    friend bool operator==(const Triple& t1, const Triple& t2);
    friend struct std::hash<Triple>;

};

template<> struct std::hash<Triple> {
    std::size_t operator()(const Triple& t) const {
        return 17 + (19 * std::hash<Variable*>()(t.x) * std::hash<Variable*>()(t.z)) + 
                (23 * std::hash<Variable*>()(t.y));
    }
};


#endif /* TRIPLE_HPP_ */