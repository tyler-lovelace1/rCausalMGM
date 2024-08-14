#include "SepsetMap.hpp"


// SepsetMap::SepsetMap(SepsetMap& map) { 
//     this->sepsets = map.sepsets; 
//     this->pValues = map.pValues;  
//     this->parents = map.parents;
//     // this->correlations = map.correlations;
//     this->returnEmptyIfNotSet = map.returnEmptyIfNotSet;
// }

/**
 * Sets the sepset for {x, y} to be z. Note that {x, y} is unordered.
 */
void SepsetMap::set(const Node& x, const Node& y, std::vector<Node>& z) {
    NodePair pair = std::minmax(x, y);
    sepsets[pair] = z;
}

/**
 * Sets the sepset for {x, y} to be z. If {x, y} is already in the SepsetMap, 
 * then z is only updated if p is greater than the previously recorded sepset.
 * Note that {x, y} is unordered.
 */
void SepsetMap::set(const Node& x, const Node& y, std::vector<Node>& z, double p) {
    NodePair pair = std::minmax(x, y);
    if (sepsets.count(pair) > 0) {
	if (p > pValues[pair]) {
	    sepsets[pair] = z;
	    pValues[pair] = p;
	}
    } else {
	sepsets[pair] = z;
	pValues[pair] = p;
    }
}

/** 
 * Removes the list associated with the pair
 */
void SepsetMap::remove(const Node& x, const Node& y) {
    NodePair pair = std::minmax(x, y);
    sepsets.erase(pair);
}

void SepsetMap::setPValue(const Node& x, const Node& y, double p) {
    NodePair pair = std::minmax(x, y);
    pValues[pair] = p;
}


/**
 * Retrieves the sepset previously set for {a, b}, or NULL if no such set was previously set.
 */
std::vector<Node> SepsetMap::get(const Node& a, const Node& b) {
    NodePair pair = std::minmax(a, b);

    // if (correlations != boost::none && correlations.get().count(pair) == 0) {
    //     return &emptyList;
    // }

    // If the pair is not set
    if (sepsets.count(pair) == 0) {
        if (returnEmptyIfNotSet) return emptyList;
        else {
	  std::vector<Node> nullVec = { Node() };
	  return nullVec;
	}
    }

    return sepsets[pair];
}

bool SepsetMap::isInSepsetMap(const Node& a, const Node& b) {
    NodePair pair = std::minmax(a, b);
    return sepsets.count(pair) > 0;
}

double SepsetMap::getPValue(const Node& x, const Node& y) {
    NodePair pair = std::minmax(x, y);
    return pValues[pair];
}

void SepsetMap::set(const Node& x, std::unordered_set<Node>& z) {
    if (parents.count(x) != 0) {
        parents[x].insert(z.begin(), z.end());
    } else {
        parents[x] = z;
    }
}

std::unordered_set<Node> SepsetMap::get(const Node& x) {
    if (parents.count(x) == 0) return std::unordered_set<Node>();
    else return parents[x];
}

void SepsetMap::addAll(SepsetMap& newSepsets) {
    sepsets.insert(newSepsets.sepsets.begin(), newSepsets.sepsets.end());
}

int SepsetMap::size() {
    return sepsets.size();
}

std::ostream& operator<<(std::ostream& os, SepsetMap& ssm) {
    // Sepsets
    os << "Sepsets:\n";
    for (std::pair<NodePair, std::vector<Node>> element : ssm.sepsets) {
        NodePair varPair = element.first;
        os << "{" << varPair.first << ", " << varPair.second << "} -> [";
        for (const Node& var : element.second) {
            os << var << ", ";
        }
        os << "]\n";
    }

    // Parents
    os << "Parents:\n";
    for (std::pair<Node, std::unordered_set<Node>> element : ssm.parents) {
        os << element.first << " -> {";
        for (const Node& par : element.second) {
            os << par << ", ";
        }
        os << "}\n";
    }

    return os;
}

bool operator==(const SepsetMap& ssm1, const SepsetMap& ssm2) {
    return ssm1.sepsets == ssm2.sepsets;
}
