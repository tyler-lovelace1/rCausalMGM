#ifndef SEPSETMAP_HPP_
#define SEPSETMAP_HPP_

#include <boost/optional/optional.hpp>
#include <boost/functional/hash.hpp>
#include "Variable.hpp"

typedef std::pair<Variable*, Variable*> VariablePair;

/**
 * <p>Stores a map from pairs of variables to separating sets--that is, for each unordered pair of variables {variable1, variable2} in a
 * graph, stores a set of variables conditional on which variable1 and variable2 are independent 
 * or stores NULL if the pair was not judged to be independent. (Note that if a sepset is non-NULL and empty,
 * that should means that the compared variables were found to be independent conditional on the empty set, whereas if a
 * sepset is NULL, that should mean that no set was found yet conditional on which the compared variables are independent.
 * So at the end of the search, a NULL sepset carries different information from an empty sepset.)</p> 
 *
 * @author Joseph Ramsey
 * @author Max Dudek (C++ conversion)
 */
class SepsetMap {

//TODO - serialization?

private:

    std::unordered_map<VariablePair, std::vector<Variable*>, boost::hash<VariablePair>> sepsets;
    std::unordered_map<VariablePair, double, boost::hash<VariablePair>> pValues;

    std::unordered_map<Variable*, std::unordered_set<Variable*>> parents;
    boost::optional<std::unordered_set<VariablePair, boost::hash<VariablePair>>> correlations = boost::none;
    bool returnEmptyIfNotSet = false;
    std::vector<Variable*> emptyList; // stays empty

public:

    SepsetMap() {}
    SepsetMap(SepsetMap& map);

    void setCorrelations(boost::optional<std::unordered_set<VariablePair, boost::hash<VariablePair>>>& pairs) { this->correlations = pairs; }
    bool isReturnEmptyIfNotSet() { return returnEmptyIfNotSet; }
    void setReturnEmptyIfNotSet(bool returnEmptyIfNotSet) { this->returnEmptyIfNotSet = returnEmptyIfNotSet; }

    /**
     * Sets the sepset for {x, y} to be z. Note that {x, y} is unordered.
     */
    void set(Variable* x, Variable* y, std::vector<Variable*>& z);

    /** 
     * Removes the list associated with the pair
     */
    void remove(Variable* x, Variable* y);

    void setPValue(Variable* x, Variable* y, double p);

    /**
     * Retrieves the sepset previously set for {a, b}, or NULL if no such set was previously set.
     */
    std::vector<Variable*>* get(Variable* a, Variable* b);

    double getPValue(Variable* x, Variable* y);

    void set(Variable* x, std::unordered_set<Variable*>& z);

    std::unordered_set<Variable*> get(Variable* x);

    void addAll(SepsetMap& newSepsets);

    int size();

    friend std::ostream& operator<<(std::ostream& os, SepsetMap& ssm);
    friend bool operator==(const SepsetMap& ssm1, const SepsetMap& ssm2);

};

#endif /* #define SEPSETMAP_HPP_ */