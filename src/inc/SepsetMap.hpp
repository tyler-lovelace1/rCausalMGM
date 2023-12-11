#ifndef SEPSETMAP_HPP_
#define SEPSETMAP_HPP_

// [[Rcpp::depends(BH)]]

// #include <boost/optional/optional.hpp>
#include <boost/functional/hash.hpp>
#include "Node.hpp"

typedef std::pair<Node, Node> NodePair;

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

    std::unordered_map<NodePair, std::vector<Node>, boost::hash<NodePair>> sepsets;
    std::unordered_map<NodePair, double, boost::hash<NodePair>> pValues;

    std::unordered_map<Node, std::unordered_set<Node>> parents;
    // boost::optional<std::unordered_set<NodePair, boost::hash<NodePair>>> correlations = boost::none;
    bool returnEmptyIfNotSet = false;
    std::vector<Node> emptyList; // stays empty

public:

    SepsetMap() {}
    SepsetMap(const SepsetMap& other) = default;
    SepsetMap(SepsetMap&& other) = default;
    SepsetMap& operator=(const SepsetMap& other) = default;
    SepsetMap& operator=(SepsetMap&& other) = default;
    ~SepsetMap() = default;

    
    // void setCorrelations(boost::optional<std::unordered_set<NodePair, boost::hash<NodePair>>>& pairs) { this->correlations = pairs; }
    bool isReturnEmptyIfNotSet() { return returnEmptyIfNotSet; }
    void setReturnEmptyIfNotSet(bool returnEmptyIfNotSet) { this->returnEmptyIfNotSet = returnEmptyIfNotSet; }

    /**
     * Sets the sepset for {x, y} to be z. Note that {x, y} is unordered.
     */
    void set(const Node& x, const Node& y, std::vector<Node>& z);

    /**
     * Sets the sepset for {x, y} to be z. If {x, y} is already in the SepsetMap, 
     * then z is only updated if p is greater than the previously recorded sepset.
     * Note that {x, y} is unordered.
     */
    void set(const Node& x, const Node& y, std::vector<Node>& z, double p);

    /** 
     * Removes the list associated with the pair
     */
    void remove(const Node& x, const Node& y);

    void setPValue(const Node& x, const Node& y, double p);

    /**
     * Retrieves the sepset previously set for {a, b}, or NULL if no such set was previously set.
     */
    std::vector<Node> get(const Node& a, const Node& b);

    bool isInSepsetMap(const Node& a, const Node& b);

    double getPValue(const Node& x, const Node& y);

    void set(const Node& x, std::unordered_set<Node>& z);

    std::unordered_set<Node> get(const Node& x);

    void addAll(SepsetMap& newSepsets);

    int size();

    friend std::ostream& operator<<(std::ostream& os, SepsetMap& ssm);
    friend bool operator==(const SepsetMap& ssm1, const SepsetMap& ssm2);

};

#endif /* #define SEPSETMAP_HPP_ */
