#include "Triple.hpp"

Triple::Triple(const Node& x, const Node& y, const Node& z) {
    this->y = y;
    if (x < z) {
	this->x = x;
	this->z = z;
    } else {
	this->x = z;
	this->z = x;
    }
}

std::string Triple::toString() const {
    std::ostringstream result;

    result << *this;

    return result.str();
}

bool Triple::contains(const Node& n) {
    return (n == x) || (n == y) || (n == z);
}

std::ostream& operator<<(std::ostream& os, const Triple& triple) {
    if (triple.x.getName() < triple.z.getName())
	os << "<" << triple.x.getName() << "," << triple.y.getName() << "," << triple.z.getName() << ">";
    else
	os << "<" << triple.z.getName() << "," << triple.y.getName() << "," << triple.x.getName() << ">";
    return os;
}

bool operator<(const Triple& t1, const Triple& t2) {
    Triple trip1(t1);
    Triple trip2(t2);
    // Rcpp::Rcout << trip1 << " < " << trip2 << " = " << (trip1.toString() < trip2.toString()) << std::endl;
    return trip1.toString() < trip2.toString();
}

bool operator==(const Triple& t1, const Triple& t2) {
    return (t1.x == t2.x && t1.y == t2.y && t1.z == t2.z) ||
           (t1.x == t2.z && t1.y == t2.y && t1.z == t2.x);
}
