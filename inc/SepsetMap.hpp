#ifndef SEPSETMAP_HPP_
#define SEPSETMAP_HPP_

#include "boost/optional/optional.hpp"
#include "Variable.hpp"

/**
 * <p>Stores a map from pairs of variables to separating sets--that is, for each unordered pair of variables {variable1, variable2} in a
 * graph, stores a set of variables conditional on which variable1 and variable2 are independent 
 * or stores boost::none if the pair was not judged to be independent. (Note that if a sepset is non-boost::none and empty,
 * that should means that the compared variables were found to be independent conditional on the empty set, whereas if a
 * sepset is boost::none, that should mean that no set was found yet conditional on which the compared variables are independent.
 * So at the end of the search, a boost::none sepset carries different information from an empty sepset.)</p> 
 *
 * @author Joseph Ramsey
 * @author Max Dudek (C++ conversion)
 */
class SepsetMap {

//TODO - serialization?

private:

    std::unordered_map<std::unordered_set<Variable*>, boost::optional<std::vector<Variable*>>> sepsets;
    std::unordered_map<std::unordered_set<Variable*>, double> pValues;

    std::unordered_map<Variable*, boost::optional<std::unordered_set<Variable*>>> parents;
    boost::optional<std::unordered_set<std::unordered_set<Variable*>>> correlations = boost::none;
    bool returnEmptyIfNotSet = false;

public:

    // TODO
    SepsetMap() {}
    SepsetMap(SepsetMap& map) { this->sepsets = map.sepsets; this->pValues = map.pValues; }


    void setCorrelations(boost::optional<std::unordered_set<std::unordered_set<Variable*>>>& pairs) { this->correlations = pairs; }
    bool isReturnEmptyIfNotSet() { return returnEmptyIfNotSet; }
    void setReturnEmptyIfNotSet(bool returnEmptyIfNotSet) { this->returnEmptyIfNotSet = returnEmptyIfNotSet; }

    /**
     * Sets the sepset for {x, y} to be z. Note that {x, y} is unordered.
     */
    void set(Variable* x, Variable* y, boost::optional<std::vector<Variable*>>& z);

    void setPValue(Variable* x, Variable* y, double p);

    /**
     * Retrieves the sepset previously set for {a, b}, or boost::none if no such set was previously set.
     */
    boost::optional<std::vector<Variable*>> get(Variable* a, Variable* b);

    double getPValue(Variable* x, Variable* y);

    void set(Variable* x, std::unordered_set<Variable*>& z);

    std::unordered_set<Variable*> get(Variable* x);

    //TODO: overload equals operator

    //TODO: overload << operator

    void addAll(SepsetMap newSepsets);

    int size();

};

#endif /* #define SEPSETMAP_HPP_ */