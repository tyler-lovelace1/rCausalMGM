#ifndef GROWSHRINKTREE_HPP_
#define GROWSHRINKTREE_HPP_

#include "Score.hpp"

struct GrowShrinkNode {
    Node node;
    
    double score = 0.0;
    mutable double shrinkScore = 0.0;
    mutable bool shrinkFlag = true;
    std::vector<Node> condSet;
    mutable std::vector<Node> shrinkCondSet;
    
    // std::set<Node> childSet;
    // std::set<GrowShrinkNode> children;

    GrowShrinkNode() {}

    GrowShrinkNode(Node node, std::vector<Node> condSet) : node(node),
							   condSet(condSet) {
	this->condSet.push_back(node);
	shrinkFlag = false;
    }

    void shrink(Score* scorer, Node target, bool verbose) {
	if (!shrinkFlag) {    
	    double oldScore = 1e20, curScore = score, tempScore;
	    bool changeFlag = false;

	    if (verbose) {
		Rcpp::Rcout << "    Score = " << score << " : { ";
		for (Node n : condSet) {
		    Rcpp::Rcout << n << " ";
		}
		Rcpp::Rcout << "}\n";
	    }
    
	    std::list<Node> activeList(condSet.begin(), condSet.end());

	    std::list<Node>::iterator removeIt;

	    while (curScore < oldScore) {
		changeFlag = false;
		oldScore = curScore;
		for (auto it = activeList.begin(); it != activeList.end(); it++) {
		    std::vector<Node> tempActive;
		    for (auto jt = activeList.begin(); jt != activeList.end(); jt++) {
			if (jt != it) {
			    tempActive.push_back(*jt);
			}
		    }
	    
		    tempScore = scorer->localScore(target, tempActive);
	    
		    if (tempScore < curScore) {
			curScore = tempScore;
			removeIt = it;
			changeFlag = true;
		    }	
		}
		if (changeFlag) {
		    activeList.erase(removeIt);
		    if (verbose) {
			Rcpp::Rcout << "    Score = " << curScore << " : { ";
			for (Node n : activeList) {
			    Rcpp::Rcout << n << " ";
			}
			Rcpp::Rcout << "}\n";
		    }
		}
	    }

	    shrinkFlag = true;
	    shrinkCondSet = std::vector<Node>(activeList.begin(), activeList.end());
	    shrinkScore = curScore;
	}
    }

    bool operator==(const GrowShrinkNode& rhs) const {
	return score == rhs.score && node == rhs.node;
    }
    
    bool operator<(const GrowShrinkNode& rhs) const {
	if (condSet.size() == rhs.condSet.size()) {
	    if (score == rhs.score) {
		return node < rhs.node;
	    }
	    return score < rhs.score;
	}
	return condSet.size() < rhs.condSet.size();
    }
};

class GrowShrinkTree {
private:
    Score* scorer;
    bool verbose = false;
    Node target;
    GrowShrinkNode root;
    std::map<GrowShrinkNode, std::set<Node>> childSet;
    std::map<GrowShrinkNode, std::set<GrowShrinkNode>> children;
    
public:

    GrowShrinkTree() {}

    GrowShrinkTree(Score* scorer, Node target);

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setPenalty(bool penalty) { this->scorer->setPenalty(penalty); }

    GrowShrinkNode grow(std::vector<Node> candidates, GrowShrinkNode& parent);

    // std::vector<Node> shrink(std::vector<Node> active, double score, double* scoreReturn = NULL);

    std::vector<Node> search(std::vector<Node> candidates, double* scoreReturn = NULL);

    std::vector<Node> search(double* scoreReturn = NULL);

    friend Rcpp::StringVector GrowShrinkTreeTest(const Rcpp::DataFrame &df, std::string target);
    friend Rcpp::StringVector GrowShrinkTreeSubSetTest(const Rcpp::DataFrame &df, std::string target, int numSub);

};

#endif /* GROWSHRINKTREE_HPP_ */
