#include "GrowShrinkTree.hpp"
#include "DegenerateGaussianScore.hpp"
#include "RegressionBicScore.hpp"

GrowShrinkTree::GrowShrinkTree(Score* scorer, Node target) : target(target) {
    if (scorer == NULL) 
        throw std::invalid_argument("scorer may not be NULL.");

    this->scorer = scorer;

    children[root] = {};
    childSet[root] = {};
}

GrowShrinkNode GrowShrinkTree::grow(std::vector<Node> candidates,
				    GrowShrinkNode& parent) {
    // std::vector<Node> active = root.condSet;
    double oldScore = 1e20, curScore = root.score, score;
    GrowShrinkNode bestNode = root;
    bool changeFlag = false;

    std::set<Node> candidateSet(candidates.begin(), candidates.end());
    // std::set<Node> activeCandidateSet(candidates.begin(), candidates.end());

    if (verbose) Rcpp::Rcout << "    Score = " << bestNode.score << " : { }\n";

    while (bestNode.score < oldScore) {
	changeFlag = false;
	oldScore = bestNode.score;

	if (childSet.count(bestNode)==0) {
	    childSet[bestNode] = {};
	    children[bestNode] = {};
	}

	for (Node n : candidateSet) {
	    if (childSet[bestNode].count(n)==0) {
		childSet[bestNode].insert(n);

		GrowShrinkNode child(n, bestNode.condSet);

		// if (verbose) {
		//     Rcpp::Rcout << "      Scoring child " << child.node << " = { ";

		//     for (Node n2 : child.condSet) {
		// 	Rcpp::Rcout << n2 << " ";
		//     }
		//     Rcpp::Rcout << "}\n";
		// }
		
		child.score = scorer->localScore(target, child.condSet);

		if (child.score < bestNode.score) {	
		    children[bestNode].insert(child);
		}

		// if (verbose) Rcpp::Rcout << "      Inserting child = { " << child.node << ", " << child.score << " }\n";
	    }
	}

	for (auto it = children[bestNode].begin(); it != children[bestNode].end(); it++) {
	    // if (verbose) Rcpp::Rcout << "      Checking child = { " << it->node << ", " << it->score << " }\n";
	    if (candidateSet.count(it->node)>0) {
		// if (verbose) Rcpp::Rcout << "        In candidate set\n";
		if (it->score < bestNode.score) {
		    // if (verbose) Rcpp::Rcout << "        Child has a better score\n";
		    changeFlag = true;
		    parent = bestNode;
		    bestNode = *it;
		    candidateSet.erase(it->node);
		}
		break;
	    }
	}
	
	if (changeFlag) {
	    if (verbose) {
		Rcpp::Rcout << "    Score = " << bestNode.score << " : { ";
		for (Node n : bestNode.condSet) {
		    Rcpp::Rcout << n << " ";
		}
		Rcpp::Rcout << "}\n";
	    }
	}
    }
    
    return bestNode;
}


// void GrowShrinkNode::shrink(Score* scorer, Node target) {

//     if (!shrinkFlag) {    
// 	double oldScore = 1e20, curScore = score, tempScore;
// 	bool changeFlag = false;

// 	if (verbose) {
// 	    Rcpp::Rcout << "    Score = " << score << " : { ";
// 	    for (Node n : condSet) {
// 		Rcpp::Rcout << n << " ";
// 	    }
// 	    Rcpp::Rcout << "}\n";
// 	}
    
// 	std::list<Node> activeList(condSet.begin(), condSet.end());

// 	std::list<Node>::iterator removeIt;

// 	while (curScore < oldScore) {
// 	    changeFlag = false;
// 	    oldScore = curScore;
// 	    for (auto it = activeList.begin(); it != activeList.end(); it++) {
// 		std::vector<Node> tempActive;
// 		for (auto jt = activeList.begin(); jt != activeList.end(); jt++) {
// 		    if (jt != it) {
// 			tempActive.push_back(*jt);
// 		    }
// 		}
	    
// 		tempScore = scorer->localScore(target, tempActive);
	    
// 		if (tempScore < curScore) {
// 		    curScore = tempScore;
// 		    removeIt = it;
// 		    changeFlag = true;
// 		}	
// 	    }
// 	    if (changeFlag) {
// 		activeList.erase(removeIt);
// 		if (verbose) {
// 		    Rcpp::Rcout << "    Score = " << curScore << " : { ";
// 		    for (Node n : activeList) {
// 			Rcpp::Rcout << n << " ";
// 		    }
// 		    Rcpp::Rcout << "}\n";
// 		}
// 	    }
// 	}

// 	shrinkFlag = true;
// 	shrinkCondSet = std::vector<Node>(activeList.begin(), activeList.end());
// 	shrinkScore = curScore;
//     }
// }

std::vector<Node> GrowShrinkTree::search(std::vector<Node> candidates,
					 double* scoreReturn) {
    
    auto it = std::remove(candidates.begin(), candidates.end(), target);

    if (it != candidates.end()) {
	candidates.erase(it, candidates.end());
    }

    it = std::remove_if(candidates.begin(), candidates.end(),
			[] (const Node& n) {
			    return n.isCensored();
			});

    if (it != candidates.end()) {
	candidates.erase(it, candidates.end());
    }
        
    if (verbose) {
	Rcpp::Rcout << "Searching for Markov Boundary of " << target.getName() << "...\n";
	Rcpp::Rcout << "  Growing...\n";
    }

    GrowShrinkNode parent(root);
    GrowShrinkNode gsNode = grow(candidates, parent);
    
    if (verbose) Rcpp::Rcout << "  Shrinking...\n";

    bool shrinkFlag = gsNode.shrinkFlag;

    gsNode.shrink(scorer, target, verbose);

    if (!shrinkFlag) {
	children[parent].erase(gsNode);
	children[parent].insert(gsNode);
    }

    if (scoreReturn != NULL) {
	*scoreReturn = gsNode.shrinkScore;
    }

    if (verbose) RcppThread::Rcout << "Finished. \n";

    return gsNode.shrinkCondSet;
}

std::vector<Node> GrowShrinkTree::search(double* scoreReturn) {
    // std::vector<Node> potentialMB;
    // for (Node n : scorer->getVariables()) {
    // 	if (n==target) continue;
    // 	if (n.isCensored()) continue;
    // 	potentialMB.push_back(n);
    // }
    return search(scorer->getVariables(), scoreReturn);
}

//[[Rcpp::export]]
Rcpp::StringVector GrowShrinkTreeTest(const Rcpp::DataFrame &df, std::string target) {
    DataSet ds = DataSet(df);
    ds.dropMissing();

    Node targetNode = ds.getVariable(target);

    std::vector<Node> regressors(ds.getVariables());
    auto it = std::remove(regressors.begin(), regressors.end(), targetNode);

    regressors.erase(it, regressors.end());

    std::vector<Node> mb;
    double score;

    if (!ds.isCensored()) {
	DegenerateGaussianScore scorer(ds, 1.0);
	GrowShrinkTree gst(&scorer, targetNode);
	gst.setVerbose(true);
	
	mb = gst.search(regressors, &score);
    } else {
        RegressionBicScore scorer(ds, 1.0);
	GrowShrinkTree gst(&scorer, targetNode);
	gst.setVerbose(true);
	
	mb = gst.search(regressors, &score);
    }
    
    RcppThread::checkUserInterrupt();

    Rcpp::StringVector _mb;

    for (Node n : mb) {
	_mb.push_back(n.getName());
    }

    _mb.attr("Score") = Rcpp::wrap(score);
    
    return _mb;
}


//[[Rcpp::export]]
Rcpp::StringVector GrowShrinkTreeSubSetTest(const Rcpp::DataFrame &df, std::string target, int numSub) {
    DataSet ds = DataSet(df);
    ds.dropMissing();

    Node targetNode = ds.getVariable(target);

    std::vector<Node> regressors(ds.getVariables());
    auto it = std::remove(regressors.begin(), regressors.end(), targetNode);

    regressors.erase(it, regressors.end());

    std::vector<Node> mb;
    double score;

    if (!ds.isCensored()) {
	DegenerateGaussianScore scorer(ds, 1.0);
	GrowShrinkTree gst(&scorer, targetNode);
	gst.setVerbose(false);

	for (int i = 0; i < numSub; i++) {
	    std::vector<Node> candidates(regressors);
	    std::random_shuffle(candidates.begin(), candidates.end(), randWrapper);
	    candidates.erase(candidates.begin() + candidates.size() / 2,
			     candidates.end());
	    mb = gst.search(candidates, &score);
	}
    } else {
        RegressionBicScore scorer(ds, 1.0);
	GrowShrinkTree gst(&scorer, targetNode);
	gst.setVerbose(false);
	
	for (int i = 0; i < numSub; i++) {
	    std::vector<Node> candidates(regressors);
	    std::random_shuffle(candidates.begin(), candidates.end(), randWrapper);
	    candidates.erase(candidates.begin() + candidates.size() / 2,
			     candidates.end());
	    mb = gst.search(candidates, &score);
	}
    }
    
    RcppThread::checkUserInterrupt();

    Rcpp::StringVector _mb;

    for (Node n : mb) {
	_mb.push_back(n.getName());
    }

    _mb.attr("Score") = Rcpp::wrap(score);
    
    return _mb;
}
