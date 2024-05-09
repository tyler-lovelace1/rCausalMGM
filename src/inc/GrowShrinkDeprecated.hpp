#ifndef GROWSHRINKDEPRECATED_HPP_
#define GROWSHRINKDEPRECATED_HPP_

#include "Score.hpp"
#include "EdgeListGraph.hpp"

class GrowShrinkDeprecated {
private:
    Score* scorer;
    bool verbose = false;
    EdgeListGraph graph;
    
    // double scoreVal;
    
public:

    GrowShrinkDeprecated() {}

    GrowShrinkDeprecated(Score* scorer);

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setPenalty(bool penalty) { this->scorer->setPenalty(penalty); }

    void setGraph(EdgeListGraph graph) { this->graph = graph; }

    std::vector<Node> grow(const Node& target, std::vector<Node> regressors,
			   double* scoreReturn = NULL);

    std::vector<Node> shrink(const Node& target, std::vector<Node> active,
			     double score, double* scoreReturn = NULL);

    std::vector<Node> search(const Node& target, std::vector<Node> regressors,
			     double* scoreReturn = NULL);

    std::vector<Node> search(const Node& target, double* scoreReturn = NULL);

    friend Rcpp::StringVector GrowShrinkDeprecatedSubSetTest(const Rcpp::DataFrame &df, std::string target, int numSub);

};

#endif /* GROWSHRINKDEPRECATED_HPP_ */
