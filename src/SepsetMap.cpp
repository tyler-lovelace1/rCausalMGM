#include "SepsetMap.hpp"


/**
 * Sets the sepset for {x, y} to be z. Note that {x, y} is unordered.
 */
void SepsetMap::set(Variable* x, Variable* y, boost::optional<std::vector<Variable*>>& z) {
    VariablePair pair = std::minmax(x, y);
    if (z == boost::none) {
        sepsets.erase(pair);
    } else {
        sepsets[pair] = z;
    }
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
        parents[x].get().insert(z.begin(), z.end());
    } else {
        parents[x].get() = z;
    }
}

std::unordered_set<Variable*> SepsetMap::get(Variable* x) {
    if (parents.count(x) == 0) return std::unordered_set<Variable*>();
    else return parents[x].get();
}

void SepsetMap::addAll(SepsetMap newSepsets) {
    sepsets.insert(newSepsets.sepsets.begin(), newSepsets.sepsets.end());
}

int SepsetMap::size() {
    return sepsets.size();
}