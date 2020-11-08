#include "Triple.hpp"

Triple::Triple(Variable* x, Variable* y, Variable* z) {
    if (x == NULL || y == NULL || z == NULL)
        throw std::invalid_argument("Nodes in Triple cannot be NULL");

    this->x = x;
    this->y = y;
    this->z = z;
}

std::ostream& operator<<(std::ostream& os, Triple& triple) {
    os << "<" << triple.x->getName() << "," << triple.y->getName() << "," << triple.z->getName() << ">";
    return os;
}

bool operator==(const Triple& t1, const Triple& t2) {
    return (t1.x == t2.x && t1.y == t2.y && t1.z == t2.z) ||
           (t1.x == t2.z && t1.y == t2.y && t1.z == t2.x);
}