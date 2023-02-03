// [[Rcpp::depends(BH,RcppThread)]]

#include "Grasp.hpp"

bool Grasp::OrderGraph::operator==(const OrderGraph& rhs) const {
    if (score == rhs.score) {
	return order == rhs.order;
    }
    return false;
}

bool Grasp::OrderGraph::operator<(const OrderGraph& rhs) const {
    if (score == rhs.score) {
	return order < rhs.order;
    }
    return score < rhs.score;
}

bool Grasp::OrderGraph::operator!=(const OrderGraph& rhs) const {
    return !(*this == rhs);
}
	
bool Grasp::OrderGraph::operator> (const OrderGraph& rhs) const {
    return rhs < *this;
}
	
bool Grasp::OrderGraph::operator<=(const OrderGraph& rhs) const {
    return !(*this > rhs);
}
	
bool Grasp::OrderGraph::operator>=(const OrderGraph& rhs) const {
    return !(*this < rhs);
}

Grasp::Grasp(DataSet& data, int threads) : ds(data), growShrink(data, threads) {
    // growShrink = GrowShrink(data, threads);

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 2;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }

    parallelism = std::max(parallelism, 1);
    
    growShrink.setPenalty(penalty);
    graph = EdgeListGraph(data.getVariables());
    // std::vector<Node> _nodes(growShrink.getVariables());
    // std::random_shuffle(_nodes.begin(), _nodes.end(), randWrapper);
    // nodes = std::list<Node>(_nodes.begin(), _nodes.end());
    // pi = OrderGraph(nodes);

    // Rcpp::Rcout << "Scoring... " << std::endl;

    // Rcpp::Rcout << "First element: " << *pi.start << std::endl;
    // Rcpp::Rcout << "Last element: " << *std::prev(pi.stop) << std::endl;
    // Rcpp::Rcout << "pi.stop == pi.order.end(): " << (pi.stop == pi.order.end()) << std::endl;
    
    // score = update(pi);

    // Rcpp::Rcout << "Initial score: " << score << std::endl;
}

void Grasp::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

void Grasp::setNumStarts(int numStarts) {
    if (numStarts < 1)
        throw std::invalid_argument("Number of starts must be >= 1.");

    if (numStarts > 100)
        throw std::invalid_argument("Number of starts must be <= 100.");

    this->numStarts = numStarts;
}

double Grasp::update(OrderGraph& tau) {
    // double score = 0;
    // tau.parentMap[tau.order.front()] = std::unordered_set<Node>();
    // tau.scoreMap[tau.order.front()] = 0;
    // auto start = tau.order.begin();
    // start++;

    // double oldScore = 0;
    // for (const Node& n : tau.order) oldScore += tau.scoreMap[n];

    double score = tau.score;
    double bic = tau.bic;
    
    for (auto it = tau.start; it != tau.stop; it++) {
	score -= tau.scoreMap[*it];
	bic -= tau.bicMap[*it];
	// if (it == tau.order.begin()) {
	//     tau.parentMap[*it] = std::unordered_set<Node>();
	//     tau.scoreMap[*it] = 0;
	// } else {
	std::vector<Node> prefix;

	if (initialGraph!=NULL) {
	    std::vector<Node> neighbors = initialGraph->getAdjacentNodes(*it);
	    std::vector<Node> possPrefix(tau.order.begin(), it);

	    std::sort(neighbors.begin(), neighbors.end());
	    std::sort(possPrefix.begin(), possPrefix.end());

	    std::set_intersection(possPrefix.begin(), possPrefix.end(),
				  neighbors.begin(), neighbors.end(),
				  std::back_inserter(prefix));
	    
	} else {
	    prefix = std::vector<Node>(tau.order.begin(), it);
	}
	
	// Rcpp::Rcout << "Target: " << it->getName() << std::endl;
	// Rcpp::Rcout << "Prefix: ";
	// for (Node n : prefix) Rcpp::Rcout << n.getName() << " ";
	// Rcpp::Rcout << std::endl;
	double localBic = 1e20;
	// std::list<Node> parents = growShrink.search(*it, prefix, &bic);
	std::list<Node> parents = growShrink.search(*it, prefix, &localBic);
	// Rcpp::Rcout << "Parents: ";
	// for (Node n : parents) Rcpp::Rcout << n.getName() << " ";
	// Rcpp::Rcout << std::endl;
	tau.parentMap[*it] = std::unordered_set<Node>(parents.begin(), parents.end());
	tau.scoreMap[*it] = parents.size();
	score += parents.size();
	tau.bicMap[*it] = localBic;
	bic += localBic;
	// Rcpp::Rcout << "Score: " << score << std::endl << std::endl;
	// }
    }
    // double newScore = 0;
    // for (const Node& n : tau.order) newScore += tau.scoreMap[n];
    // Rcpp::Rcout << "Old Score: " << oldScore << std::endl;
    // Rcpp::Rcout << "Full Score: " << newScore << std::endl;
    // Rcpp::Rcout << "Score: " << score << std::endl << std::endl;
    tau.score = score;
    tau.bic = bic;
    return score;
}


// double Grasp::update(OrderGraph& tau) {
//     // double score = 0;
//     // tau.parentMap[tau.order.front()] = std::unordered_set<Node>();
//     // tau.scoreMap[tau.order.front()] = 0;
//     // auto start = tau.order.begin();
//     // start++;

//     // double oldScore = 0;
//     // for (const Node& n : tau.order) oldScore += tau.scoreMap[n];

//     double score = tau.score;
//     double bic = tau.bic;

//     RcppThread::ThreadPool pool(parallelism);

//     auto updateTask = [&] (const std::list<Node>::iterator& it) {
// 			  {
// 			      std::lock_guard<std::mutex> scoreLock(scoreMutex);
// 			      score -= tau.scoreMap[*it];
// 			      bic -= tau.bicMap[*it];
// 			  }
// 			  // if (it == tau.order.begin()) {
// 			  //     tau.parentMap[*it] = std::unordered_set<Node>();
// 			  //     tau.scoreMap[*it] = 0;
// 			  // } else {
// 			  std::vector<Node> prefix;

// 			  if (initialGraph!=NULL) {
// 			      std::vector<Node> neighbors(initialGraph->getAdjacentNodes(*it));
			      
// 			      std::vector<Node> possPrefix(tau.order.begin(), it);

// 			      std::sort(neighbors.begin(), neighbors.end());
// 			      std::sort(possPrefix.begin(), possPrefix.end());

// 			      std::set_intersection(possPrefix.begin(), possPrefix.end(),
// 						    neighbors.begin(), neighbors.end(),
// 						    std::back_inserter(prefix));
	    
// 			  } else {
// 			      prefix = std::vector<Node>(tau.order.begin(), it);
// 			  }
	
// 			  // RcppThread::Rcout << "Target: " << it->getName() << std::endl;
// 			  // RcppThread::Rcout << "Prefix: ";
// 			  // for (Node n : prefix) RcppThread::Rcout << n.getName() << " ";
// 			  // RcppThread::Rcout << std::endl;
// 			  double localBic = 1e20;

// 			  // GrowShrink localGrowShrink(growShrink);
// 			  // std::list<Node> parents = growShrink.search(*it, prefix, &bic);
// 			  std::list<Node> parents = growShrink.search(*it, prefix,
// 								      &localBic);
// 			  // RcppThread::Rcout << "Parents: ";
// 			  // for (Node n : parents) RcppThread::Rcout << n.getName() << " ";
// 			  // RcppThread::Rcout << std::endl;
// 			  {
// 			      std::lock_guard<std::mutex> scoreLock(scoreMutex);
// 			      tau.parentMap[*it] = std::unordered_set<Node>(parents.begin(),
// 									    parents.end());
// 			      tau.scoreMap[*it] = parents.size();
// 			      score += parents.size();
// 			      tau.bicMap[*it] = localBic;
// 			      bic += localBic;
// 			      // RcppThread::Rcout << "Score: " << score << std::endl << std::endl;
// 			  }
			  
// 			  // }
// 		      };
    
//     for (auto it = tau.start; it != tau.stop; it++) {
// 	// RcppThread::Rcout << "Loading task for: " << it->getName() << std::endl;
// 	pool.push(updateTask, it);
//     }
    
//     pool.join();

//     // double newScore = 0;
//     // for (const Node& n : tau.order) newScore += tau.scoreMap[n];
//     // Rcpp::Rcout << "Old Score: " << oldScore << std::endl;
//     // Rcpp::Rcout << "Full Score: " << newScore << std::endl;
//     // Rcpp::Rcout << "Score: " << score << std::endl << std::endl;
//     tau.score = score;
//     tau.bic = bic;
//     return score;
// }

// double Grasp::update(OrderGraph& tau,
// 		     std::list<Node>::iterator start,
// 		     std::list<Node>::iterator stop) {
//     double score = tau.score;
    
//     for (auto it = start; it != stop; it++) {
// 	if (tau.scoreMap.count(*it)) score -= tau.scoreMap[*it];
// 	std::vector<Node> prefix(tau.order.begin(), it);
// 	// Rcpp::Rcout << "Target: " << it->getName() << std::endl;
// 	// Rcpp::Rcout << "Prefix: ";
// 	// for (Node n : prefix) Rcpp::Rcout << n.getName() << " ";
// 	// Rcpp::Rcout << std::endl;
// 	std::list<Node> parents = growShrink.search(*it, prefix);
// 	// Rcpp::Rcout << "Parents: ";
// 	// for (Node n : parents) Rcpp::Rcout << n.getName() << " ";
// 	// Rcpp::Rcout << std::endl;
// 	tau.parentMap[*it] = std::unordered_set<Node>(parents.begin(), parents.end());
// 	tau.scoreMap[*it] = parents.size();
// 	score += parents.size();
// 	// Rcpp::Rcout << "Score: " << score << std::endl;
//     }
//     // Rcpp::Rcout << "Score: " << score << std::endl;
//     tau.score = score;
//     return score;
// }


bool Grasp::graspDFS(int tier) {
    OrderGraph tau = pi;
    OrderGraph oldTau;
    double oldScore = tau.score;
    int curDepth = 0;
    std::set<Grasp::OrderGraph> visited = {};
    std::set<std::set<NodePair>> visitedTuckSets = {};
    std::set<NodePair> tucked = {};
    int N = ds.getNumRows();
    double buffer = std::log(N);

    std::stack<Grasp::OrderGraph> graphStack;
    std::stack<std::set<NodePair>> tuckStack;
    std::stack<int> depthStack;
    
    bool first = true;

    graphStack.push(tau);
    tuckStack.push(tucked);
    depthStack.push(0);

    while (!graphStack.empty()) {
	tau = graphStack.top();
	graphStack.pop();
	
	tucked = tuckStack.top();
	tuckStack.pop();

	curDepth = depthStack.top();
	depthStack.pop();

	// if (tier==0) curDepth /= 4;

	if (visitedTuckSets.count(tucked) || visited.count(tau) || curDepth > depth)
	    continue;

	// if (visitedTuckSets.count(tucked) || visited.count(tau) ||
	//     (tier > 0 && curDepth > depth) || (tier == 0 && curDepth > 2 * depth))
	//     continue;

	// if (tier==0) curDepth *= 4;
	
	if (!first) {
	    update(tau);

	    // Rcpp::Rcout << "    Old Score = " << oldTau.score << std::endl;
	    // Rcpp::Rcout << "    Old BIC = " << oldTau.bic << std::endl;

	    if (tau.score < pi.score || (tau.score==pi.score && tau.bic < pi.bic)) {
		pi = tau;
		bestGraphs.clear();
		bestGraphs.insert(pi);
		// growShrink.clearHistory();
		if (verbose) Rcpp::Rcout << "    Score = " << pi.score << ", BIC = " << pi.bic << "\r";
		return true;
	    }

	    if (tau.score > pi.score) continue;

	    // double dBic = tau.bic - pi.bic;
	    // double postProbTau = std::exp(-0.5 * dBic);
	    // postProbTau = postProbTau / (1 + postProbTau);

	    // // if (verbose) Rcpp::Rcout << "    Score = " << tau.score << ", BIC = " << tau.bic << ", p(tau|data) = " << postProbTau <<  "\n";

	    // if (tau.bic != pi.bic && postProbTau < 0.05) continue;

	    if (verbose) Rcpp::Rcout << "    Score = " << tau.score << ", BIC = " << tau.bic << "\r";

	    if (tau.bic == pi.bic) bestGraphs.insert(tau);
	
	    // Rcpp::Rcout << "    Current Score = " << tau.score << std::endl;
	    // Rcpp::Rcout << "    Current BIC = " << tau.bic << std::endl;
	    
	    // if (tau.score == oldTau.score && tau.bic < oldTau.bic) {
	    // 	pi = tau;
	    // 	Rcpp::Rcout << "    New BIC = " << pi.bic << std::endl;
	    // 	return true;
	    // }
	}

	first = false;

	// Rcpp::Rcout << "  Depth: " << curDepth << std::endl;

	// Rcpp::Rcout << "    Old Score = " << oldTau.score << std::endl;

	visited.insert(tau);

	visitedTuckSets.insert(tucked);

	oldTau = tau;

	// bool dfsBreak = false;
	
	for (const Node& y : tau.getShuffledNodes()) {
	    if (tau.parentMap[y].size()==0) continue;
	    std::unordered_set<Node> parents(tau.parentMap[y]);
	    std::unordered_set<Node> ancestors;
	    tau.getAncestors(y, ancestors);
	    // Rcpp::Rcout << "Node " << y.getName() << "\n  Ancestors: ";
	    // for (Node n : ancestors) Rcpp::Rcout << n.getName() << " ";
	    // Rcpp::Rcout << "Node " << y.getName() << std::endl;
	    
	    for (const Node& x : parents) {
		// Rcpp::Rcout << "Parent " << x.getName() << "\n";
	    
		bool covered = tau.isCovered(x, y);
		bool singular = true;

		NodePair tuck = std::minmax(x, y);

		// Rcpp::Rcout << "  isCovered? " << covered << std::endl;
		// Rcpp::Rcout << "  already tucked? " << tucked.count(tuck) << std::endl;

		if (tier==0 && !covered) continue;
		if (covered && tucked.count(tuck)) continue;

		// if (tier == 0) {
		//     Rcpp::Rcout << "  Tucking covered edge " << x.getName() << " --> "
		// 		<< y.getName() << std::endl;
		// }

		// if (tier == 2) {
		//     Rcpp::Rcout << "  Tucking edge " << x.getName() << " --> "
		// 		<< y.getName() << std::endl;
		// }

		auto start = std::find(tau.order.begin(), tau.order.end(), x);
		auto stop = std::find(start, tau.order.end(), y);

		// std::list<Node> newOrder(tau.order.begin(), start);
		// std::list<Node> gamma, gammac;
		// std::list<Node> newOrderEnd(std::next(stop), tau.order.end());

		// start++;

		// for (auto it = start; it != stop; it++) {
		//     if (ancestors.count(*it)) {
		// 	if (tau.parentMap[*it].count(x)) {
		// 	    // Rcpp::Rcout << "    nonsingular\n";
		// 	    singular = false;
		// 	    if (tier==1) break;
		// 	}
		// 	gamma.push_back(*it);
		//     } else {
		// 	gammac.push_back(*it);
		//     }
		// }

		// newOrder.splice(newOrder.end(), gamma);
		// newOrder.push_back(y);
		// newOrder.push_back(x);
		// newOrder.splice(newOrder.end(), gammac);
		// newOrder.splice(newOrder.end(), newOrderEnd);

		// tau.order = newOrder;
		// tau.start = gamma.empty() ? std::find(newOrder.begin(), newOrder.end(), y) : std::find(newOrder.begin(), newOrder.end(), gamma.front());
		// tau.stop = newOrderEnd.empty() ? newOrder.end() : std::find(newOrder.begin(), newOrder.end(), newOrderEnd.front());

		start = tau.order.insert(start, y);
		stop = tau.order.erase(stop);
		auto jt = std::next(start, 2);
		auto insertLoc = start;
		// start--;

		// Rcpp::Rcout << "  Ancestors: ";
		// for (Node n : ancestors) Rcpp::Rcout << n.getName() << " ";
		// Rcpp::Rcout << std::endl;

		bool first = true;
		
		while (jt != stop) {
		    // Rcpp::Rcout << "    Order: ";
		    // for (auto it = start; it != stop; it++) Rcpp::Rcout << it->getName() << " ";
		    // Rcpp::Rcout << std::endl;
		    // Rcpp::Rcout << "Current: " << jt->getName() << std::endl;
		    if (ancestors.count(*jt)) {
			// Rcpp::Rcout << "  in ancestors\n";
			if (tau.parentMap[*jt].count(x)) {
			    // Rcpp::Rcout << "    nonsingular\n";
			    singular = false;
			    if (tier==1) break;
			}
			
			if (first) {
			    start = tau.order.insert(insertLoc, *jt);
			    first = false;
			} else {
			    tau.order.insert(insertLoc, *jt);
			}
			
			jt = tau.order.erase(jt);
			// Rcpp::Rcout << "  order updated\n";
		    } else {
			jt++;
		    }
		}

		if (tier==1 && !singular) {
		    tau.order = oldTau.order;
		    continue;
		}

		// if (tier == 1) {
		//     Rcpp::Rcout << "  Tucking singular edge " << x.getName() << " --> "
		// 		<< y.getName() << std::endl;
		// }

		tau.start = start;
		tau.stop = stop;
		
		// update(tau, start, stop);

		// if (tau.score < oldTau.score) {
		//     pi = tau;
		//     Rcpp::Rcout << "  Tucked " << x.getName() << " & " << y.getName()
		// 		<< ": New Score = " << pi.score << std::endl;
		//     return true;
		// }

		// if (tau.score == oldTau.score && curDepth < depth) {
		//     Rcpp::Rcout << "  Tucked " << x.getName() << " & " << y.getName()
		// 		<< ": Equal Score = " << tau.score << std::endl;
		    
		    // visited.insert(tau);
		tucked.insert(tuck);
		graphStack.push(tau);
		tuckStack.push(tucked);
		depthStack.push(curDepth+1);
		tucked.erase(tuck);
		    // dfsBreak = true;
		    // break;
		    // }

		tau = oldTau;
		
	    }
	    // if (dfsBreak) {
	    // 	Rcpp::Rcout << "At break, stack size = " << graphStack.size() << std::endl;
	    // 	break;
	    // }
	}

	// Rcpp::Rcout << "At DFS iteration, stack size = " << graphStack.size() << std::endl;
    }

    return false;
}

// bool Grasp::graspDFSr(OrderGraph& tau, int curDepth, int tier,
// 		      std::set<NodePair>& tucks,
// 		      std::set<std::set<NodePair>>& dfsHistory) {

//     OrderGraph oldTau = tau;

//     Rcpp::Rcout << "Depth = " << curDepth << std::endl;

//     for (const Node& y : tau.getShuffledNodes()) {
// 	if (tau.parentMap[y].size()==0) continue;
// 	std::unordered_set<Node> parents(tau.parentMap[y]);
// 	std::unordered_set<Node> ancestors;
// 	tau.getAncestors(y, ancestors);
// 	// Rcpp::Rcout << "Node " << y.getName() << "\n  Ancestors: ";
// 	// for (Node n : ancestors) Rcpp::Rcout << n.getName() << " ";
// 	// Rcpp::Rcout << "Node " << y.getName() << std::endl;
	    
// 	for (const Node& x : parents) {
// 	    // Rcpp::Rcout << "Parent " << x.getName() << "\n";
	    
// 	    bool covered = tau.isCovered(x, y);
// 	    bool singular = true;

// 	    NodePair tuck = std::minmax(x, y);

// 	    // Rcpp::Rcout << "  isCovered? " << covered << std::endl;
// 	    // Rcpp::Rcout << "  already tucked? " << tucked.count(tuck) << std::endl;

// 	    if (tier==0 && !covered) continue;
// 	    if (covered && tucks.count(tuck)) continue;

// 	    // Rcpp::Rcout << "  Tucking " << x.getName() << " & " << y.getName()
// 	    // 		<< std::endl;

// 	    if (tier == 0) {
// 		Rcpp::Rcout << "  Tucking covered edge " << x.getName() << " --> "
// 			    << y.getName() << std::endl;
// 	    }

// 	    if (tier == 2) {
// 		Rcpp::Rcout << "  Tucking edge " << x.getName() << " --> "
// 			    << y.getName() << std::endl;
// 	    }
		    

// 	    auto start = std::find(tau.order.begin(), tau.order.end(), x);
// 	    auto stop = std::find(start, tau.order.end(), y);

// 	    start = tau.order.insert(start, y);
// 	    stop = tau.order.erase(stop);
// 	    auto jt = start;
// 	    std::advance(jt, 2);
// 	    auto insertLoc = start;
// 	    // start--;

// 	    // Rcpp::Rcout << "  Ancestors: ";
// 	    // for (Node n : ancestors) Rcpp::Rcout << n.getName() << " ";
// 	    // Rcpp::Rcout << std::endl;

// 	    bool first = true;
		
// 	    while (jt != stop) {
// 		// Rcpp::Rcout << "    Order: ";
// 		// for (auto it = start; it != stop; it++) Rcpp::Rcout << it->getName() << " ";
// 		// Rcpp::Rcout << std::endl;
// 		// Rcpp::Rcout << "Current: " << jt->getName() << std::endl;
// 		if (ancestors.count(*jt)) {
// 		    // Rcpp::Rcout << "  in ancestors\n";
// 		    if (tau.parentMap[*jt].count(x)) {
// 			// Rcpp::Rcout << "    nonsingular\n";
// 			singular = false;
// 			if (tier==1) break;
// 		    }
			
// 		    if (first) {
// 			start = tau.order.insert(insertLoc, *jt);
// 			first = false;
// 		    } else {
// 			tau.order.insert(insertLoc, *jt);
// 		    }
			
// 		    jt = tau.order.erase(jt);
// 		    // Rcpp::Rcout << "  order updated\n";
// 		} else {
// 		    jt++;
// 		}
// 	    }

// 	    if (tier==1 && !singular) {
// 		tau.order = oldTau.order;
// 		continue;
// 	    }

// 	    if (tier == 1) {
// 		Rcpp::Rcout << "  Tucking singular edge " << x.getName() << " --> "
// 			    << y.getName() << std::endl;
// 	    }
		
// 	    update(tau, start, stop);

// 	    if (tau.score < oldTau.score) {
// 		pi = tau;
// 		Rcpp::Rcout << "    New Score = " << pi.score << std::endl;
// 		return true;
// 	    }

// 	    if (tau.score == oldTau.score && curDepth < depth) {
// 		Rcpp::Rcout << "    Equal Score = " << tau.score << std::endl;
		    
// 		// visited.insert(tau);
// 		tucks.insert(tuck);
// 		if (dfsHistory.count(tucks)==0) {
// 		    dfsHistory.insert(tucks);
// 		    if (graspDFSr(tau, curDepth+1, tier, tucks, dfsHistory))
// 			return true;
// 		    Rcpp::Rcout << "Depth = " << curDepth << std::endl;
// 		}
// 		tucks.erase(tuck);
// 		// dfsBreak = true;
// 		// break;
// 	    }

// 	    tau = oldTau;
	    
// 	}
//     }

//     return false;
// }

std::map<EdgeListGraph, std::pair<int, double>> Grasp::search() {

    int inputDepth = depth;

    std::vector<Node> _nodes(growShrink.getVariables());
    std::random_shuffle(_nodes.begin(), _nodes.end(), randWrapper);
    nodes = std::list<Node>(_nodes.begin(), _nodes.end());
    pi = OrderGraph(nodes);
    score = update(pi);

    OrderGraph bestOrder = pi;

    std::map<EdgeListGraph, std::pair<int, double>> cpdagMap;

    for (int i = 0; i < numStarts; i++) {
	if (verbose) Rcpp::Rcout << "Run " << i+1 << "...\n";
	
	if (i > 0) {
	    std::vector<Node> _nodes(growShrink.getVariables());
	    std::random_shuffle(_nodes.begin(), _nodes.end(), randWrapper);
	    nodes = std::list<Node>(_nodes.begin(), _nodes.end());
	    pi = OrderGraph(nodes);
	    update(pi);
	}
	
	for (int tier = 0; tier < 3; tier++) {
	    double oldScore = nodes.size() * nodes.size();
	    double curScore = pi.score;
	    double oldBIC = 1e20;
	    double curBIC = pi.bic;
	    // growShrink.clearHistory();

	    if (tier==0)  depth = std::max(inputDepth, 4);
	    else          depth = inputDepth;

	    if (verbose) Rcpp::Rcout << "  Running GRaSP" << tier << "...\n";
	    while (curScore < oldScore || curBIC < oldBIC) {
		oldScore = curScore;
		oldBIC = curBIC;
		// OrderGraph tau = pi;
		// std::set<NodePair> tucks;
		// std::set<std::set<NodePair>> dfsHistory;
		// dfsHistory.insert(tucks);
		// graspDFSr(tau, 0, tier, tucks, dfsHistory);
		graspDFS(tier);
		curScore = pi.score;
		curBIC = pi.bic;
	    }
	    // double oldScore = std::numeric_limits<double>::max();
	    if (verbose) Rcpp::Rcout << "\n";
	}

	if (pi.score < bestOrder.score) {
	    bestOrder = pi;

	    cpdagMap.clear();

	    // std::set<EdgeListGraph> runCpdagSet;

	    for (const OrderGraph order : bestGraphs) {
	    
		graph = order.toGraph();

		std::ostringstream alg;
		if (initialGraph==NULL) {
		    alg << "GRaSP";
		} else {
		    alg << initialGraph->getAlgorithm() << "-" << "GRaSP";
		}
		graph.setAlgorithm(alg.str());

		if (cpdagMap.count(graph)>0) {
		    cpdagMap[graph].first++;
		    cpdagMap[graph].second = std::min(cpdagMap[graph].second, order.bic);
		} else {
		    cpdagMap[graph] = std::pair<int, double>(1, order.bic);
		}
	    }
	} else if (pi.score == bestOrder.score) {
	    for (const OrderGraph order : bestGraphs) {
	    
		graph = order.toGraph();

		std::ostringstream alg;
		if (initialGraph==NULL) {
		    alg << "GRaSP";
		} else {
		    alg << initialGraph->getAlgorithm() << "-" << "GRaSP";
		}
		graph.setAlgorithm(alg.str());

		if (cpdagMap.count(graph)>0) {
		    cpdagMap[graph].first++;
		    cpdagMap[graph].second = std::min(cpdagMap[graph].second, order.bic);
		} else {
		    cpdagMap[graph] = std::pair<int, double>(1, order.bic);
		}
	    }
	}

	bestGraphs.clear();
	// else if (pi.score == bestOrder.score && pi.bic < bestOrder.bic) bestOrder = pi;
    }

    // graph = bestOrder.toGraph();

    // std::ostringstream alg;
    // if (initialGraph==NULL) {
    // 	alg << "GRaSP";
    // } else {
    // 	alg << initialGraph->getAlgorithm() << "-" << "GRaSP";
    // }
    // graph.setAlgorithm(alg.str());
        
    return cpdagMap;
}
