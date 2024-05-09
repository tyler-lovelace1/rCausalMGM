#ifndef GROWSHRINK_HPP_
#define GROWSHRINK_HPP_

#include "Score.hpp"

class GrowShrink {
private:
    Score* scorer;
    bool verbose = false;
    // double scoreVal;
    
public:

    GrowShrink() {}

    GrowShrink(Score* scorer);

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setPenalty(bool penalty) { this->scorer->setPenalty(penalty); }

    std::vector<Node> grow(const Node& target, std::vector<Node> regressors,
			   double* scoreReturn = NULL);

    std::vector<Node> shrink(const Node& target, std::vector<Node> active,
			     double score, double* scoreReturn = NULL);

    std::vector<Node> search(const Node& target, std::vector<Node> regressors,
			     double* scoreReturn = NULL);

    std::vector<Node> search(const Node& target, double* scoreReturn = NULL);

    friend Rcpp::StringVector GrowShrinkSubSetTest(const Rcpp::DataFrame &df, std::string target, int numSub);

};

#endif /* GROWSHRINK_HPP_ */
