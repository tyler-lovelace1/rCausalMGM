#include "SepsetMap.hpp"


/**
 * Sets the sepset for {x, y} to be z. Note that {x, y} is unordered.
 */
void SepsetMap::set(Variable* x, Variable* y, std::vector<Variable*>& z) {
    VariablePair pair = std::minmax(x, y);
    sepsets[pair] = z;
}

/** 
 * Removes the list associated with the pair
 */
void SepsetMap::remove(Variable* x, Variable* y) {
    VariablePair pair = std::minmax(x, y);
    sepsets.erase(pair);
}

void SepsetMap::setPValue(Variable* x, Variable* y, double p) {
    VariablePair pair = std::minmax(x, y);
    pValues[pair] = p;
}


/**
 * Retrieves the sepset previously set for {a, b}, or null if no such set was previously set.
 */
boost::optional<std::vector<Variable*>> SepsetMap::get(Variable* a, Variable* b) {
    VariablePair pair = std::minmax(a, b);

    if (correlations == boost::none && correlations.get().count(pair) == 0) {
        return std::vector<Variable*>();
    }

    // If the pair is not set
    if (sepsets.count(pair) == 0) {
        if (returnEmptyIfNotSet) return std::vector<Variable*>();
        else return boost::none;
    }

    return sepsets[pair];
}

double SepsetMap::getPValue(Variable* x, Variable* y) {
    VariablePair pair = std::minmax(x, y);
    return pValues[pair];
}

void SepsetMap::set(Variable* x, std::unordered_set<Variable*>& z) {
    if (parents.count(x) != 0) {
        parents[x].insert(z.begin(), z.end());
    } else {
        parents[x] = z;
    }
}

std::unordered_set<Variable*> SepsetMap::get(Variable* x) {
    if (parents.count(x) == 0) return std::unordered_set<Variable*>();
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
    for (std::pair<VariablePair, std::vector<Variable*>> element : ssm.sepsets) {
        VariablePair varPair = element.first;
        os << "{" << varPair.first->getName() << ", " << varPair.second->getName() << "} -> [";
        for (Variable* var : element.second) {
            os << var->getName() << ", ";
        }
        os << "]\n";
    }

    // Parents
    os << "Parents:\n";
    for (std::pair<Variable*, std::unordered_set<Variable*>> element : ssm.parents) {
        os << element.first->getName() << " -> {";
        for (Variable* par : element.second) {
            os << par->getName() << ", ";
        }
        os << "}\n";
    }

    return os;
}

bool operator==(const SepsetMap& ssm1, const SepsetMap& ssm2) {
    return ssm1.sepsets == ssm2.sepsets;
}
