#include "Node.hpp"

bool operator==(const Node& n1, const Node& n2) {
    return n1.name == n2.name;
}

bool operator!=(const Node& n1, const Node& n2) {
    return !(n1 == n2);
}

bool operator< (const Node& n1, const Node& n2) {
    return n1.name < n2.name;
}

bool operator<=(const Node& n1, const Node& n2) {
    return n1 < n2 || n1 == n2;
}

bool operator> (const Node& n1, const Node& n2) {
    return n1.name > n2.name;
}

bool operator>=(const Node& n1, const Node& n2) {
    return n1 > n2 || n1 == n2;
}

std::ostream& operator<<(std::ostream& os, const Node& n) {
    os << n.getName();
    return os;
}

std::size_t hash_value(const Node& n) {
    std::size_t res = 17;
    res = res * 31 + std::hash<std::string>()(n.name);
//    res = res * 43 + std::hash<DataType>()(n.type);
    return res;
}
