#include "GrowShrinkTree.hpp"
#include "DegenerateGaussianScore.hpp"
#include "RegressionBicScore.hpp"

GrowShrinkTree::GrowShrinkTree(Score* scorer, Node target) : target(target) {
    if (scorer == NULL) 
        throw std::invalid_argument("scorer may not be NULL.");

    this->scorer = scorer;

    int idx = 0;
    for (Node n : scorer->getVariables()) {
	if (n == target) continue;
	node2idx[n] = idx;
	idx++;
    }

    // RcppThread::Rcout << scorer->getVariables().size()-1 << std::endl;

    // RcppThread::Rcout << idx << std::endl;

    numNodes = idx;

    // root = new GrowShrinkNode(numNodes);
    root.reset(new GrowShrinkNode(numNodes));
}

std::vector<Node> GrowShrinkTree::growShrink(std::vector<Node>& candidates, double* scoreReturn) {

    if (verbose) RcppThread::Rcout << "  Growing...\n";

    double oldScore = 1e20, curScore = root->score, score;
    GrowShrinkNode* bestNode = root.get();
    bool changeFlag = false;

    std::set<Node> candidateSet(candidates.begin(), candidates.end());
    // std::set<Node> activeCandidateSet(candidates.begin(), candidates.end());

    // std::set<std::shared_ptr<GrowShrinkNode>,GrowShrinkNode::cmp_ptr> candidateChildren;

    if (verbose) RcppThread::Rcout << "    Score = " << bestNode->score << " : { }\n";

    while (bestNode->score < oldScore) {
	changeFlag = false;
	oldScore = bestNode->score;

	// if (childSet.count(bestNode)==0) {
	//     childSet[bestNode] = {};
	//     children[bestNode] = {};
	// }

	// newChildren[bestNode] = {};

	for (Node n : candidateSet) {
	    if (!bestNode->visited[node2idx[n]].load()) {		
		GrowShrinkNode* child = new GrowShrinkNode(n, bestNode->condSet, numNodes);

		// if (verbose) {
		//     RcppThread::Rcout << "      Scoring child " << child->node << " = { ";

		//     for (Node n2 : child->condSet) {
		// 	RcppThread::Rcout << n2 << " ";
		//     }
		//     RcppThread::Rcout << "}\n";
		// }
		
		child->score = scorer->localScore(target, child->condSet);

		// childSet[bestNode].insert(n);
		if (child->score < bestNode->score) {
		    // std::lock_guard<std::mutex> treeLock(treeMutex);
		    // newChildren[bestNode].insert(child);
		    bestNode->addChild(child);
		    // if (verbose) RcppThread::Rcout << "      Inserting child = { " << child->node << ", " << child->score << " }\n";
		} else {
		    // child.reset();
		    delete child;
		}

		bestNode->visited[node2idx[n]].store(true);
		// bestNode->visited[node2idx[n]] = true;
		// 
	    }
	}

	// std::set<Node> candidateChildren;
	// {
	//     std::lock_guard<std::mutex> treeLock(treeMutex);
	//     candidateChildren = std::move(std::set<Node>(children[bestNode].begin(), children[bestNode].end()));
	// }

	// candidateChildren.insert(newChildren[bestNode].begin(), newChildren[bestNode].end());
	// {
	//     std::lock_guard<std::mutex> treeLock(treeMutex);
	// candidateChildren = bestNode->getChildren();
	// }

	for (auto it = bestNode->children.begin(); it != bestNode->children.end(); it++) {
	    // if (verbose) RcppThread::Rcout << "      Checking child = { " << (*it)->node << ", " << (*it)->score << " }\n";
	    if (candidateSet.count((*it)->node)>0) {
		// if (verbose) RcppThread::Rcout << "        In candidate set\n";
		// if ((*it)->score < bestNode->score) {
		//     if (verbose) RcppThread::Rcout << "        Child has a better score\n";
		changeFlag = true;
		// parent = bestNode;
		bestNode = *it;
		candidateSet.erase((*it)->node);
		// }
		break;
	    }
	}
	
	if (changeFlag) {
	    if (verbose) {
		RcppThread::Rcout << "    Score = " << bestNode->score << " : { ";
		for (Node n : bestNode->condSet) {
		    RcppThread::Rcout << n << " ";
		}
		RcppThread::Rcout << "}\n";
	    }
	}
    }

    if (verbose) RcppThread::Rcout << "  Shrinking...\n";

    bestNode->shrink(scorer, target, verbose);

    if (scoreReturn != NULL) {
	*scoreReturn = bestNode->shrinkScore.load();
    }
    
    return bestNode->shrinkCondSet;
}


// void GrowShrinkNode::shrink(Score* scorer, Node target) {

//     if (!shrinkFlag) {    
// 	double oldScore = 1e20, curScore = score, tempScore;
// 	bool changeFlag = false;

// 	if (verbose) {
// 	    RcppThread::Rcout << "    Score = " << score << " : { ";
// 	    for (Node n : condSet) {
// 		RcppThread::Rcout << n << " ";
// 	    }
// 	    RcppThread::Rcout << "}\n";
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
// 		    RcppThread::Rcout << "    Score = " << curScore << " : { ";
// 		    for (Node n : activeList) {
// 			RcppThread::Rcout << n << " ";
// 		    }
// 		    RcppThread::Rcout << "}\n";
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
        
    if (verbose)
	RcppThread::Rcout << "Searching for Markov Boundary of " << target << "...\n";

    // std::set<GrowShrinkNode> newChildren;
    // std::shared_ptr<GrowShrinkNode> gsNode = grow(candidates);
    
    // if (verbose) RcppThread::Rcout << "  Shrinking...\n";

    std::vector<Node> condSet = growShrink(candidates, scoreReturn);    

    // {
    // 	std::lock_guard<std::mutex> treeLock(treeMutex);
	
    // 	for (auto it = newChildren.begin(); it != newChildren.end(); it++) {
    // 	    children[it->first].insert(it->second.begin(), it->second.end());
    // 	}

    // 	auto bestIt = children[growOut[0]].find(growOut[1]);
    // 	bestIt->shrink(scorer, target, verbose);
    // 	if (scoreReturn != NULL) {
    // 	    *scoreReturn = bestIt->shrinkScore;
    // 	}
    // 	condSet = bestIt->shrinkCondSet;
    // }

    // bool shrinkFlag = gsNode.shrinkFlag;

    // gsNode->shrink(scorer, target, verbose);

    // if (!shrinkFlag) {
    // 	children[parent].erase(gsNode);
    // 	children[parent].insert(gsNode);
    // }

    // if (scoreReturn != NULL) {
    // 	*scoreReturn = gsNode->shrinkScore.load();
    // }

    // condSet = gsNode->shrinkCondSet;

    if (verbose) RcppThread::Rcout << "Finished. \n";

    return condSet;
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

// no export [[Rcpp::export]]
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


// no export  [[Rcpp::export]]
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
	    // if ((i+1) % 50 == 0) {
	    // 	gst.reset();
	    // }
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


//  no export [[Rcpp::export]]
std::vector<double> GrowShrinkTreeParallelSubSetTest(const Rcpp::DataFrame &df, std::string target, int numSub, int threads) {
    DataSet ds = DataSet(df);
    ds.dropMissing();

    Node targetNode = ds.getVariable(target);

    std::vector<Node> regressors(ds.getVariables());
    auto it = std::remove(regressors.begin(), regressors.end(), targetNode);

    regressors.erase(it, regressors.end());

    std::vector<Node> mb;
    double score;

    RcppThread::ThreadPool pool(std::max(1, threads));

    DegenerateGaussianScore scorer(ds, 1.0);
    GrowShrinkTree* gst = new GrowShrinkTree(&scorer, targetNode);
    gst->setVerbose(true);

    auto gsTask = [&](std::vector<Node> candidates) {
	double score;
	RcppThread::checkUserInterrupt();
	gst->search(candidates, &score);
	return score;
    };

    std::vector<std::future<double>> futures(numSub);
    std::vector<double> scores;

    for (int i = 0; i < numSub; i++) {
	std::vector<Node> candidates(regressors);
	std::random_shuffle(candidates.begin(), candidates.end(), randWrapper);
	candidates.erase(candidates.begin() + candidates.size() / 2,
			 candidates.end());
	futures[i] = pool.pushReturn(gsTask, candidates);
    }

    for (int i = 0; i < numSub; i++) {
	scores[i] = futures[i].get();
    }

    pool.join();
    
    delete gst;
    
    return scores;
}
