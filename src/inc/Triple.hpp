#ifndef TRIPLE_HPP_
#define TRIPLE_HPP_

#include "Node.hpp"

class Triple {
private:
    

public:

    Node x;
    Node y;
    Node z;

    Triple() {}
    Triple(const Node& x, const Node& y, const Node& z);
    Triple(const Triple& other) = default;
    Triple& operator=(const Triple& other) = default;
    Triple(Triple&& other) = default;
    Triple& operator=(Triple&& other) = default;
    ~Triple() = default;

    Node getX() { return x; }
    Node getY() { return y; }
    Node getZ() { return z; }

    std::string toString() const;

    bool contains(const Node& n);

    friend std::ostream& operator<<(std::ostream& os, const Triple& triple);
    friend bool operator<(const Triple& t1, const Triple& t2);
    friend bool operator==(const Triple& t1, const Triple& t2);
    friend struct std::hash<Triple>;

};

template<> struct std::hash<Triple> {
    std::size_t operator()(const Triple& t) const {
        return 17 + (19 * std::hash<Node>()(t.x) * std::hash<Node>()(t.z)) + 
                (23 * std::hash<Node>()(t.y));
    }
};


#endif /* TRIPLE_HPP_ */
