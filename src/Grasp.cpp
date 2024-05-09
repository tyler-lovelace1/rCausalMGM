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


std::list<Node> Grasp::initializeRandom() {
    std::vector<std::set<Node>> tiers = knowledge.getTiers();

    if (tiers.empty()) {
	std::set<Node> _nodes(nodes.begin(), nodes.end());
	tiers.push_back(_nodes);
    }

    std::list<Node> nodeList;

    // int tierIdx=0;
    for (const std::set<Node>& tier : tiers) {
	// tierIdx++;
	std::vector<Node> _tier(tier.begin(), tier.end());

	// Rcpp::Rcout << "  Tier " << tierIdx << ":\n";
	// for (const Node& n: _tier) Rcpp::Rcout << n << " ";
	// Rcpp::Rcout << std::endl;
	
	std::random_shuffle(_tier.begin(), _tier.end(), randWrapper);
	nodeList.insert(nodeList.end(), _tier.begin(), _tier.end());
    }

    std::set<Node> set1(this->nodes.begin(), this->nodes.end());
    std::set<Node> set2(nodeList.begin(), nodeList.end());

    if (set1 != set2)
	throw std::runtime_error("GRASP random initialization does not contain every node.");

    return nodeList;
}

std::list<Node> Grasp::initializeBOSS() {

    std::list<Node> nodeList = initializeRandom();

    OrderGraph tau(nodeList);

    updateParallel(tau);

    double oldBic = 1e20;
    double oldScore = nodeList.size() * nodeList.size();
    double bic = tau.bic;
    double score = tau.score;

    if (verbose) Rcpp::Rcout << "  Initializing with BOSS...\n";

    int iter = 0;
    if (verbose) Rcpp::Rcout << "\r    Iter:  " << iter << "    Edges = " << tau.score << ", Score = " << tau.bic << "      ";
	

    while (oldBic - bic > 1e-6) {
	oldBic = bic;
	oldScore = score;
	nodeList = std::list<Node>(tau.order.begin(), tau.order.end());
	for (Node n : nodeList) {
	    // if (verbose) Rcpp::Rcout << "\n      Node " << n << "...";
	    tau = bestMove(n, tau);
	    // if (verbose) Rcpp::Rcout << "\r    Edges = " << tau.score << ", Score = " << tau.bic;
	    RcppThread::checkUserInterrupt();
	}
	score = tau.score;
	bic = tau.bic;
	iter++;

	if (verbose) Rcpp::Rcout << "\r    Iter:  " << iter << "    Edges = " << tau.score << ", Score = " << tau.bic << "      ";
	// break;
    }

    if (verbose) Rcpp::Rcout << "\n";
    
    return tau.order;
}


// Grasp::OrderGraph Grasp::internalBOSS(Node n, OrderGraph tau) {
//     // OrderGraph oldTau = tau;

//     auto oldIt = std::find(tau.order.begin(), tau.order.end(), n);

//     // RcppThread::Rcout << "      Node n = " << *oldIt << "\n";

//     RcppThread::ThreadPool pool(parallelism);

//     auto bossTask = [this, &n, &tau] (int index) {

// 	OrderGraph newTau(tau);

// 	// RcppThread::Rcout << "      Node n = " << n << " index = " << index << "\n";

// 	std::list<Node>::iterator oldIt = std::find(newTau.order.begin(), newTau.order.end(), n);
// 	std::list<Node>::iterator newIt = newTau.order.begin();

// 	// RcppThread::Rcout << "Tau order: \n";
// 	// for (auto node : newTau.order) {
// 	//     RcppThread::Rcout << node << " ";
// 	// }	
// 	// RcppThread::Rcout << "\n";

// 	// RcppThread::Rcout << "      Node n = " << *oldIt << "\n";
		
// 	bool beforeFlag = false;
// 	for (int i = 0; i < index; i++) {
// 	    if (newIt == oldIt) 
// 		beforeFlag = true;
// 	    newIt++;
// 	}
	
//  	// Node n2 = *newIt;

// 	// RcppThread::Rcout << "      Node at new position = " << *newIt << "\n";
	
// 	// RcppThread::Rcout << "      beforeFlag = " << beforeFlag << "\n";
	
// 	if (oldIt == newIt) {
// 	    return newTau;
// 	} else if (beforeFlag) {
// 	    newTau.start = newTau.order.erase(oldIt);
// 	    // newIt = std::find(newTau.order.begin(), newTau.order.end(), n2);
// 	    // newIt++;
// 	    newTau.stop = newTau.order.insert(newIt, n);
// 	    newTau.stop++;
// 	    // newTau.stop--;
// 	} else {
// 	    newTau.stop = newTau.order.erase(oldIt);
// 	    newTau.start = newTau.order.insert(newIt, n);
// 	    // newTau.start++;
// 	}

// 	// RcppThread::Rcout << "New tau order: \n";
// 	// for (auto node : newTau.order) {
// 	//     RcppThread::Rcout << node << " ";
// 	// }	
// 	// RcppThread::Rcout << "\n";

// 	// RcppThread::Rcout << "      Node at start position = " << *newTau.start << "\n";

// 	// RcppThread::Rcout << "      Node at stop position = " << *newTau.stop << "\n";

// 	// newTau.start = newTau.order.begin();
// 	// newTau.stop = newTau.order.end();

// 	// RcppThread::Rcout << "      Old Score = " << newTau.bic << "\n";

// 	update(newTau);

// 	// RcppThread::Rcout << "      New Score = " << newTau.bic << "\n";
	
// 	return newTau;
//     };

//     // Rcpp::Rcout << "    Scoring new orders...\n";
    
//     std::vector<std::future<OrderGraph>> futures(tau.order.size());
//     std::vector<OrderGraph> results(tau.order.size());
//     for (int i = 0; i < tau.order.size(); ++i)
// 	futures[i] = pool.pushReturn(bossTask, i);
//     for (int i = 0; i < tau.order.size(); ++i)
// 	results[i] = futures[i].get();
//     pool.join();

//     // for (int i = 0; i < tau.order.size(); ++i)
//     //  	results[i] = bossTask(i);

//     // Rcpp::Rcout << "    Scoring of new orders complete\n";

//     double minBic = tau.bic;
//     int minIdx = 0;
//     auto minIt = tau.order.begin();
//     for (int i = 0; i < tau.order.size(); i++) {
// 	if (minIt == oldIt)
// 	    break;
// 	minIdx++;
// 	minIt++;
//     }
    
//     for (int i = 0; i < tau.order.size(); ++i) {
// 	if (results[i].bic < minBic) {
// 	    minIdx = i;
// 	    minBic = results[i].bic;
// 	}
//     }

//     return results.at(minIdx);
// }


Grasp::OrderGraph Grasp::bestMove(Node n, OrderGraph tau) {
    // OrderGraph oldTau = tau;

    auto oldIt = std::find(tau.order.begin(), tau.order.end(), n);

    // RcppThread::Rcout << "\n      Node n = " << *oldIt << "\n";

    RcppThread::ThreadPool pool(parallelism);

    auto scoreTask = [this] (Node node, std::vector<Node> prefix) {
	double localBic = 1e20;
	std::vector<Node> parents;

	// RcppThread::Rcout << "Scoring " << node << "\n";

	RcppThread::checkUserInterrupt();

	{
	    std::lock_guard<std::mutex> gstLock(gstMutexMap[node]);
	    parents = gstMap[node]->search(prefix, &localBic);
	}
	
	return localBic;
    };
	
    // Rcpp::Rcout << "    Scoring new orders...\n";

    std::vector<std::future<double>> futureWith(tau.order.size());
    std::vector<std::future<double>> futureWithout(tau.order.size());
    std::vector<std::future<double>> futureScores(tau.order.size());
	
    std::vector<double> with(tau.order.size());
    std::vector<double> without(tau.order.size());
    std::vector<double> scores(tau.order.size());

    int currIdx = 0;
    
    std::vector<Node> prefix;
    auto jt = tau.order.begin();

    // futureScores[0] = pool.pushReturn(scoreTask, n, prefix);
    // Rcpp::Rcout << "Order: ";
    for (int i = 0; i < tau.order.size(); ++i) {

	// Rcpp::Rcout << *jt << " ";
	
	if (*jt == n) {
	    currIdx = i;
	    jt++;
	    continue;
	}
	
	futureWithout[i] = pool.pushReturn(scoreTask, *jt, prefix);
	futureScores[i] = pool.pushReturn(scoreTask, n, prefix);

	prefix.push_back(n);
	futureWith[i] = pool.pushReturn(scoreTask, *jt, prefix);
	prefix.pop_back();

	prefix.push_back(*jt);
	
	jt++;
    }

    // Rcpp::Rcout << "\n";

    for (int i = 0; i < tau.order.size(); ++i) {
	if (i == currIdx) continue;
	with[i] = futureWith[i].get();
	without[i] = futureWithout[i].get();
	scores[i] = futureScores[i].get();
    }

    pool.join();

    double runningScore = 0.0;
    
    // Rcpp::Rcout << "  currIdx = " << currIdx << "\n";

    // Rcpp::Rcout << "Original Order: ";
    // for (Node n : tau.order) {
    // 	Rcpp::Rcout << n << " ";
    // }
    // Rcpp::Rcout << "\n";
    
    for (int i = tau.order.size() - 1; i >= 0; i--) {
	if (i == currIdx) continue;
	scores[i] += runningScore;
	runningScore += with[i];
    }

    runningScore = 0.0;
    for (int i = 0; i < tau.order.size(); i++) {
	if (i == currIdx) continue;
	runningScore += without[i];
	scores[i] += runningScore;
    }

    scores[currIdx] = tau.bic;

    // Rcpp::Rcout << "    Scoring of new orders complete\n";

    double minBic = tau.bic;
    int minIdx = currIdx;
    
    for (int i = 0; i < tau.order.size(); ++i) {
	if (scores[i] + 1e-6 < minBic) {
	    minIdx = i;
	    minBic = scores[i];
	}
    }

    if (currIdx == minIdx) return tau;

    std::list<Node>::iterator newIt = std::next(tau.order.begin(), minIdx);

    // Rcpp::Rcout << "\n    Best Score = " << scores[minIdx] << "\n";
    // Rcpp::Rcout << "    Best Score Index = " << minIdx << "\n";
    // Rcpp::Rcout << "    Node at Index = " << *newIt << "\n";

    if (currIdx < minIdx) {
	// Rcpp::Rcout << "currIdx < minIdx" << std::endl;
	tau.start = tau.order.erase(oldIt);
	newIt++;
	tau.stop = tau.order.insert(newIt, n);
	tau.stop++;
    } else {
	// Rcpp::Rcout << "currIdx >= minIdx" << std::endl;
	tau.stop = tau.order.erase(oldIt);
	tau.start = tau.order.insert(newIt, n);
    }

    // Rcpp::Rcout << "New Order: ";
    // for (Node n : tau.order) {
    // 	Rcpp::Rcout << n << " ";
    // }

    // Rcpp::Rcout << "\nstart = " << *tau.start << std::endl;

    // if (tau.stop == tau.order.end()) {
    // 	Rcpp::Rcout << "stop = END" << std::endl;
    // } else {
    // 	Rcpp::Rcout << "stop = " << *tau.stop << std::endl;
    // }
    
    updateParallel(tau);

    // Rcpp::Rcout << "    Best Score After Update = " << tau.bic << "\n";

    return tau;
}



std::list<Node> Grasp::initializeMinDegree() {
    std::vector<std::set<Node>> tiers = knowledge.getTiers();

    if (tiers.empty()) {
	std::set<Node> _nodes(nodes.begin(), nodes.end());
	tiers.push_back(_nodes);
    }

    std::list<Node> nodeList;

    if (initialGraph==NULL) {

	// int tierIdx=0;
	for (const std::set<Node>& tier : tiers) {
	    // tierIdx++;
	    std::vector<Node> _tier(tier.begin(), tier.end());

	    EdgeListGraph tierUndir(_tier);

	    for (int i = 0; i < _tier.size(); i++) {
		std::vector<Node> parents = gstMap[_tier.at(i)]->search(_tier);
		for (Node n : parents) {
		    tierUndir.addUndirectedEdge(_tier.at(i), n);
		}
	    }

	    Rcpp::Rcout << "GrowShrink undirected graph:\n" << tierUndir;

	    initialGraph = new EdgeListGraph(tierUndir);

	    std::list<Node> mdOrder = minDegreeAlgorithm(tierUndir);
	
	    // std::random_shuffle(_tier.begin(), _tier.end(), randWrapper);
	
	    nodeList.insert(nodeList.end(), mdOrder.begin(), mdOrder.end());

	    
	}

    } else {

	std::list<Node> mdOrder = minDegreeAlgorithm(*initialGraph);
		
	nodeList.insert(nodeList.end(), mdOrder.begin(), mdOrder.end());

    }

    std::set<Node> set1(this->nodes.begin(), this->nodes.end());
    std::set<Node> set2(nodeList.begin(), nodeList.end());

    if (set1 != set2)
	throw std::runtime_error("GRASP random initialization does not contain every node.");

    return nodeList;
}

std::list<Node> Grasp::minDegreeAlgorithm(EdgeListGraph graph) {

    int p = graph.getNumNodes();
    int minDegree;
    int degree;
    Node mdNode;
    std::set<Node> active;
    std::list<Node> mdOrder;

    for (Node n : graph.getNodes()) {
	active.insert(n);
    }

    while (!active.empty()) {
	minDegree = 2 * p;
	for (Node n : active) {
	    degree = graph.getDegree(n);
	    if (degree < minDegree) {
		minDegree = degree;
		mdNode = n;
	    } else if (degree==minDegree && R::runif(0,1) < 0.5) {
		mdNode = n;
	    }
	}

	Rcpp::Rcout << "Min Degree = " << minDegree << " : " << mdNode << std::endl;

	std::vector<Node> adjNodes = graph.getAdjacentNodes(mdNode);

	if (minDegree > 1) {
	    for (int i = 1; i < adjNodes.size(); i++) {
		for (int j = 0; j < i; j++) {
		    graph.addUndirectedEdge(adjNodes.at(i), adjNodes.at(j));
		}

		// std::vector<Node> possibleAdj(graph.getAdjacentNodes(adjNodes.at(i)));
		// possibleAdj.insert(possibleAdj.end(), adjNodes.begin(), adjNodes.end());
		// std::vector<Node> mb = growShrink.search(adjNodes.at(i), possibleAdj);
		
		// for (Node n : mb) {
		//     graph.addUndirectedEdge(adjNodes.at(i), n);
		// }
	    }
	}

	graph.removeNode(mdNode);

	active.erase(mdNode);

	mdOrder.push_front(mdNode);

	// Rcpp::Rcout << "\nUpdated Graph:\n\n" << graph << "\n\n";
    }

    Rcpp::Rcout << "Minimum degree order:\n";
    for (Node n : mdOrder) Rcpp::Rcout << n << " ";
    Rcpp::Rcout << "\n";

    return mdOrder;
}

Grasp::Grasp(Score *scorer, int threads) {
    // growShrink = GrowShrink(data, threads);

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
    }

    parallelism = std::max(parallelism, 1);

    std::vector<Node> _nodes = scorer->getVariables();

    this->penalty = scorer->getPenalty();
    
    // growShrink.setPenalty(penalty);
    graph = EdgeListGraph(_nodes);
    nodes = std::list<Node>(_nodes.begin(), _nodes.end());

    for (Node node : nodes) {
	gstMap[node] = std::make_unique<GrowShrinkTree>(scorer, node);
	gstMutexMap[node];
    }
    
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
    double score = tau.score;
    double bic = tau.bic;

    // std::vector<RcppThread::Thread> threads;
    
    for (auto it = tau.start; it != tau.stop; it++) {
	score -= tau.scoreMap[*it];
	bic -= tau.bicMap[*it];
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
	
	double localBic = 1e20;
	std::vector<Node> parents;

	{
	    std::lock_guard<std::mutex> gstLock(gstMutexMap[*it]);
	    parents = gstMap[*it]->search(prefix, &localBic);
	}
	
	tau.parentMap[*it] = std::unordered_set<Node>(parents.begin(), parents.end());
	tau.scoreMap[*it] = parents.size();
	score += parents.size();
	tau.bicMap[*it] = localBic;
	bic += localBic;
    }
    tau.score = score;
    tau.bic = bic;
    return score;
}

double Grasp::updateParallel(OrderGraph& tau) {
    // double score = 0;
    // tau.parentMap[tau.order.front()] = std::unordered_set<Node>();
    // tau.scoreMap[tau.order.front()] = 0;
    // auto start = tau.order.begin();
    // start++;

    // double oldScore = 0;
    // for (const Node& n : tau.order) oldScore += tau.scoreMap[n];

    double score = tau.score;
    double bic = tau.bic;

    RcppThread::ThreadPool pool(parallelism);

    auto updateTask = [&] (std::list<Node>::iterator it) {
			  {
			      std::lock_guard<std::mutex> scoreLock(scoreMutex);
			      score -= tau.scoreMap[*it];
			      bic -= tau.bicMap[*it];
			  }
			  
			  std::vector<Node> prefix;

			  if (initialGraph!=NULL) {
			      std::vector<Node> neighbors(initialGraph->getAdjacentNodes(*it));
			      
			      std::vector<Node> possPrefix(tau.order.begin(), it);

			      std::sort(neighbors.begin(), neighbors.end());
			      std::sort(possPrefix.begin(), possPrefix.end());

			      std::set_intersection(possPrefix.begin(), possPrefix.end(),
						    neighbors.begin(), neighbors.end(),
						    std::back_inserter(prefix));
	    
			  } else {
			      prefix = std::vector<Node>(tau.order.begin(), it);
			  }
	

			  double localBic = 1e20;

			  std::vector<Node> parents;

			  
			  // std::lock_guard<std::mutex> gstLock(gstMutexMap[*it]);
			  parents = gstMap[*it]->search(prefix, &localBic);
			  
			  {
			      std::lock_guard<std::mutex> scoreLock(scoreMutex);
			      tau.parentMap[*it] = std::unordered_set<Node>(parents.begin(),
									    parents.end());
			      tau.scoreMap[*it] = parents.size();
			      score += parents.size();
			      tau.bicMap[*it] = localBic;
			      bic += localBic;
			  }
			  
		      };
    
    for (auto it = tau.start; it != tau.stop; it++) {
	pool.push(updateTask, it);
    }
    
    pool.join();

    tau.score = score;
    tau.bic = bic;
    return score;
}

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
    double oldBic = tau.bic;
    int curDepth = 0;
    double bicThresh = std::log(10.0); // Bayes Factor K = sqrt(10),
				           // substantially worse
				           // prune models that are
				           // not substantively the
				           // same as the current best
				           // model
    // std::set<std::list<Node>> visited = {};
    // std::set<EdgeListGraph> visitedDAGs = {};
    std::set<std::set<NodePair>> visitedTuckSets = {};
    std::set<NodePair> tucked = {};

    std::stack<Grasp::OrderGraph> graphStack;
    std::stack<std::set<NodePair>> tuckStack;
    std::stack<int> depthStack;
    
    bool first = false;

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

	// if ((curDepth > depth) || visited.count(tau.order)>0 || visitedTuckSets.count(tucked)>0 || !tau.isValid(knowledge))
	//     continue;

	if (curDepth > depth || !tau.isValid(knowledge))
	    continue;

	RcppThread::checkUserInterrupt();
	
	if (!first) {
	    
	    updateParallel(tau);

	    // if (scoreEdges) {

	    // 	if (tau.score < pi.score || (tau.score == pi.score && tau.bic < pi.bic)) {
	    // 	    pi = tau;
	    // 	    // bestGraphs.clear();
	    // 	    // bestGraphs.insert(pi);
		
	    // 	    // if (verbose) Rcpp::Rcout << "\n    Depth = " << curDepth << ", Score = " << pi.score << ", Score = " << pi.bic;
	    // 	    return true;
	    // 	}
	    // 	// 
	    // 	if (tau.score > pi.score || (tau.score == pi.score && tau.bic - pi.bic > bicThresh)) continue;
	    // } else {
		
	    if (tau.bic - pi.bic < -1e-6 || (tau.bic - pi.bic < 1e-6 && tau.score - pi.score < -0.5)) {
		pi = tau;
		// bestGraphs.clear();
		// bestGraphs.insert(pi);
		
		// if (verbose) Rcpp::Rcout << "\n    Depth = " << curDepth << ", Score = " << pi.score << ", Score = " << pi.bic;
		return true;
	    }

	    if (tau.bic - pi.bic > 1e-6) continue;

	    // double dBic = tau.bic - pi.bic; //  + 1e-6

	    // double acceptProb = std::exp(-dBic/2.0);

	    // // // acceptProb /= (1 + acceptProb);

	    // if (acceptProb < Rcpp::runif(1)[0]) {
	    //     // Rcpp::Rcout << "\n    Branch Pruned:  Depth = " << curDepth << ", Score = " << tau.score << ", Score = " << tau.bic << ", deltaScore = " << dBic;
	    //     continue;
	    // }
	    // }

	    // if (verbose) Rcpp::Rcout << "\n    Depth = " << curDepth << ", Score = " << tau.score << ", Score = " << tau.bic;

	    // if (tau.bic == pi.bic) bestGraphs.insert(tau);
	
	}

	first = false;

	// visited.insert(tau.order);

	// visitedTuckSets.insert(tucked);

	oldTau = tau;

	// Rcpp::Rcout << "\n    Original Node Order:\n";

	// for (Node n : tau.order) {
	//     Rcpp::Rcout << n << " ";
	// }
	// Rcpp::Rcout << std::endl;

	std::vector<Node> nodeList(tau.getShuffledNodes());
	
	for (const Node& y : nodeList) {
	    if (tau.parentMap[y].size()==0) continue;
	    std::unordered_set<Node> parents(tau.parentMap[y]);
	    std::unordered_set<Node> ancestors;
	    tau.getAncestors(y, ancestors);
	    
	    for (const Node& x : parents) {
	    
		bool covered = tau.isCovered(x, y);
		bool singular = true;

		NodePair tuck = std::minmax(x, y);

		if (tier==0 && !covered) continue;
		if (covered && tucked.count(tuck)) continue;


		auto start = std::find(tau.order.begin(), tau.order.end(), x);
		auto stop = std::find(start, tau.order.end(), y);

		start = tau.order.insert(start, y);
		stop = tau.order.erase(stop);
		auto jt = std::next(start, 2);
		auto insertLoc = start;

		bool first = true;
		
		while (jt != stop) {
		    if (ancestors.count(*jt)) {
			if (tau.parentMap[*jt].count(x)) {
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
		    } else {
			jt++;
		    }
		}

		if (tier==1 && !singular) {
		    tau.order = oldTau.order;
		    continue;
		}

		tau.start = start;
		tau.stop = stop;
		
		tucked.insert(tuck);
		
		if (visitedTuckSets.count(tucked)==0) {
		    visitedTuckSets.insert(tucked);
		    graphStack.push(tau);
		    tuckStack.push(tucked);
		    depthStack.push(curDepth+1);
		    // Rcpp::Rcout << "    Performed tuck : { " << x << ", " << y << " }\n";

		    // Rcpp::Rcout << "    New Node Order:\n";

		    // for (Node n : tau.order) {
		    //     Rcpp::Rcout << n << " ";
		    // }
		    // Rcpp::Rcout << std::endl;

		    // Rcpp::Rcout << "    start = " << *tau.start << std::endl;
		    // if (tau.stop != tau.order.end()) {
		    // 	Rcpp::Rcout << "    stop = " << *tau.stop << std::endl;
		    // } else {
		    // 	Rcpp::Rcout << "    stop = end\n";
		    // }
		}

		tucked.erase(tuck);

		// Rcpp::Rcout << "    Performed tuck : { " << x << ", " << y << " }\n";

		// Rcpp::Rcout << "    Old graph:\n" << oldTau.toDAG() << std::endl;

		// Rcpp::Rcout << "    New Node Order:\n";

		// for (Node n : tau.order) {
		//     Rcpp::Rcout << n << " ";
		// }
		// Rcpp::Rcout << std::endl;

		tau = oldTau;
		
	    }
	}
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

EdgeListGraph Grasp::search() {

    int inputDepth = depth;

    // std::vector<Node> _nodes(growShrink.getVariables());
    // std::random_shuffle(_nodes.begin(), _nodes.end(), randWrapper);
    // nodes = initializeRandom();
    // pi = OrderGraph(nodes);
    // score = update(pi);

    double bicThresh = std::log(10.0);

    OrderGraph bestOrder;

    // std::map<EdgeListGraph, std::pair<int, double>> cpdagMap;

    for (int i = 0; i < numStarts; i++) {
	if (verbose) Rcpp::Rcout << "Run " << i+1 << "...\n";
	
	if (bossInit) {
	    nodes = initializeBOSS();
	} else {
	  nodes = initializeRandom();
	}
	pi = OrderGraph(nodes);
	score = updateParallel(pi);
	
	for (int tier = 0; tier < 3; tier++) {
	    
	    if (tier==0)  depth = std::max(inputDepth, 4);
	    else          depth = inputDepth;

	    if (verbose) Rcpp::Rcout << "  Running GRaSP" << tier << "...\n";
	    bool improved = true;
	    while (improved) {

		if (verbose) Rcpp::Rcout << "\r    Edges = " << pi.score << ", Score = " << pi.bic << "      ";
		improved = graspDFS(tier);
	    }
	    if (verbose) Rcpp::Rcout << "\n";
	}

	if (pi.bic - bestOrder.bic < -1e-6 || (pi.bic - bestOrder.bic < 1e-6 && pi.score - bestOrder.score < -0.5)) {
	// if (pi.score < bestOrder.score || (pi.score == bestOrder.score && pi.bic < bestOrder.bic)) {
	    bestOrder = pi;

	//     cpdagMap.clear();

	//     // std::set<EdgeListGraph> runCpdagSet;

	//     for (const OrderGraph order : bestGraphs) {

	// 	// Rcpp::Rcout << "Best Order:\n";
	// 	// for (const Node& n : order.order) Rcpp::Rcout << n << " ";
	// 	// Rcpp::Rcout << "\n\n";
	    
	// 	graph = order.toGraph();

	// 	std::ostringstream alg;
	// 	if (initialGraph==NULL) {
	// 	    alg << "GRaSP";
	// 	} else {
	// 	    alg << initialGraph->getAlgorithm() << "-" << "GRaSP";
	// 	}
	// 	graph.setAlgorithm(alg.str());

	// 	if (cpdagMap.count(graph)>0) {
	// 	    cpdagMap[graph].first++;
	// 	    cpdagMap[graph].second = std::min(cpdagMap[graph].second, order.bic);
	// 	} else {
	// 	    cpdagMap[graph] = std::pair<int, double>(1, order.bic);
	// 	}
	//     }
	// } else if (pi.bic == bestOrder.bic) {
	//     for (const OrderGraph order : bestGraphs) {

	// 	// Rcpp::Rcout << "Best Order:\n";
	// 	// for (const Node& n : order.order) Rcpp::Rcout << n << " ";
	// 	// Rcpp::Rcout << "\n\n";
		
	// 	graph = order.toGraph();

	// 	std::ostringstream alg;
	// 	if (initialGraph==NULL) {
	// 	    alg << "GRaSP";
	// 	} else {
	// 	    alg << initialGraph->getAlgorithm() << "-" << "GRaSP";
	// 	}
	// 	graph.setAlgorithm(alg.str());

	// 	if (cpdagMap.count(graph)>0) {
	// 	    cpdagMap[graph].first++;
	// 	    cpdagMap[graph].second = std::min(cpdagMap[graph].second, order.bic);
	// 	} else {
	// 	    cpdagMap[graph] = std::pair<int, double>(1, order.bic);
	// 	}
	//     }
	}

	// bestGraphs.clear();
	resetTrees();
	// else if (pi.score == bestOrder.score && pi.bic < bestOrder.bic) bestOrder = pi;
    }

    graph = bestOrder.toGraph();

    std::ostringstream alg;
    if (initialGraph==NULL) {
	alg << "GRaSP";
    } else {
	alg << initialGraph->getAlgorithm() << "-" << "GRaSP";
    }
    graph.setAlgorithm(alg.str());

    graph.setScore(bestOrder.bic);

    graph.setHyperParam("penalty", arma::vec({ penalty }));

    // graph = bestOrder.toGraph();

    // std::ostringstream alg;
    // if (initialGraph==NULL) {
    // 	alg << "GRaSP";
    // } else {
    // 	alg << initialGraph->getAlgorithm() << "-" << "GRaSP";
    // }
    // graph.setAlgorithm(alg.str());
        
    return graph;
}
