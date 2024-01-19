#include "GrowShrink.hpp"

GrowShrink::GrowShrink(Score* scorer) {
    if (scorer == NULL) 
        throw std::invalid_argument("scorer may not be NULL.");

    this->scorer = scorer;
}

std::vector<Node> GrowShrink::grow(const Node& target,
				   std::vector<Node> regressors,
				   double* scoreReturn) {
    std::vector<Node> active;
    double oldScore = 1e20, curScore = 0.0, score;
    Node bestNode;
    bool changeFlag = false;

    if (verbose) Rcpp::Rcout << "    Score = " << curScore << " : { }\n";

    while (curScore < oldScore) {
	changeFlag = false;
	oldScore = curScore;
	for (auto it = regressors.begin(); it != regressors.end(); it++) {
	    if (std::find(active.begin(), active.end(), *it)==active.end()) {
		active.push_back(*it);
		score = scorer->localScore(target, active);
		if (score < curScore) {
		    curScore = score;
		    bestNode = *it;
		    changeFlag = true;
		}
		active.pop_back();
	    }
	}
	if (changeFlag) {
	    active.push_back(bestNode);
	    if (verbose) {
		Rcpp::Rcout << "    Score = " << curScore << " : { ";
		for (Node n : active) {
		    Rcpp::Rcout << n << " ";
		}
		Rcpp::Rcout << "}\n";
	    }
	}
    }
    if (scoreReturn != NULL) {
	*scoreReturn = curScore;
    }
    return active;
}


std::vector<Node> GrowShrink::shrink(const Node& target, std::vector<Node> active,
				     double curScore, double* scoreReturn) {
    double oldScore = 1e20, score;
    // Node bestNode;
    bool changeFlag = false;

    if (verbose) {
	Rcpp::Rcout << "    Score = " << curScore << " : { ";
	for (Node n : active) {
	    Rcpp::Rcout << n << " ";
	}
	Rcpp::Rcout << "}\n";
    }
    
    std::list<Node> activeList(active.begin(), active.end());

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
	    
	    score = scorer->localScore(target, tempActive);
	    
	    if (score < curScore) {
		curScore = score;
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
    if (scoreReturn != NULL) {
	*scoreReturn = curScore;
    }
    return std::vector<Node>(activeList.begin(), activeList.end());
}

std::vector<Node> GrowShrink::search(const Node& target,
				     std::vector<Node> regressors,
				     double* scoreReturn) {
    
    auto it = std::remove(regressors.begin(), regressors.end(), target);
    regressors.erase(it, regressors.end());

    double score = 0.0;
    
    if (verbose) {
	Rcpp::Rcout << "Searching for Markov Boundary of " << target.getName() << "...\n";
	Rcpp::Rcout << "  Growing...\n";
    }
    
    std::vector<Node> active = grow(target, regressors, &score);
    if (verbose) Rcpp::Rcout << "  Shrinking...\n";
    active = shrink(target, active, score, &score);

    if (scoreReturn != NULL) {
	*scoreReturn = score;
    }

    if (verbose) RcppThread::Rcout << "Finished. \n";

    return active;
}

std::vector<Node> GrowShrink::search(const Node& target, double* scoreReturn) {
    return search(target, scorer->getVariables(), scoreReturn);
}
