#include "Triple.hpp"

std::ostream& operator<<(std::ostream& os, const Triple& triple) {
    os << "<" << triple.x->getName() << "," << triple.y->getName() << "," << triple.z->getName() << ">";
    return os;
}

bool operator==(const Triple& t1, const Triple& t2) {
    return (t1.x == t2.x && t1.y == t2.y && t1.z == t2.z) ||
           (t1.x == t2.z && t1.y == t2.y && t1.z == t2.x);
}