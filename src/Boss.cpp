// [[Rcpp::depends(BH,RcppThread)]]

#include "Boss.hpp"

bool Boss::OrderGraph::operator==(const OrderGraph& rhs) const {
    if (score == rhs.score) {
	return order == rhs.order;
    }
    return false;
}

bool Boss::OrderGraph::operator<(const OrderGraph& rhs) const {
    if (score == rhs.score) {
	return order < rhs.order;
    }
    return score < rhs.score;
}

bool Boss::OrderGraph::operator!=(const OrderGraph& rhs) const {
    return !(*this == rhs);
}
	
bool Boss::OrderGraph::operator> (const OrderGraph& rhs) const {
    return rhs < *this;
}
	
bool Boss::OrderGraph::operator<=(const OrderGraph& rhs) const {
    return !(*this > rhs);
}
	
bool Boss::OrderGraph::operator>=(const OrderGraph& rhs) const {
    return !(*this < rhs);
}


std::list<Node> Boss::initializeRandom() {
    std::vector<std::set<Node>> tiers = knowledge.getTiers();

    if (tiers.empty()) {
	std::set<Node> _nodes(nodes.begin(), nodes.end());
	tiers.push_back(_nodes);
    }

    std::list<Node> nodeList;

    for (const std::set<Node>& tier : tiers) {
	std::vector<Node> _tier(tier.begin(), tier.end());

	std::random_shuffle(_tier.begin(), _tier.end(), randWrapper);
	nodeList.insert(nodeList.end(), _tier.begin(), _tier.end());
    }

    std::set<Node> set1(this->nodes.begin(), this->nodes.end());
    std::set<Node> set2(nodeList.begin(), nodeList.end());

    if (set1 != set2)
	throw std::runtime_error("BOSS random initialization does not contain every node.");

    return nodeList;
}


Boss::Boss(Score *scorer, int threads) : taskQueue(MAX_QUEUE_SIZE) {

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
    this->scorer = scorer;
    
    graph = EdgeListGraph(_nodes);
    nodes = std::list<Node>(_nodes.begin(), _nodes.end());

    for (Node node : nodes) {
	gstMap[node] = std::make_unique<GrowShrinkTree>(scorer, node);
	gstMutexMap[node];
    }
    
}

void Boss::setNumStarts(int numStarts) {
    if (numStarts < 1)
        throw std::invalid_argument("Number of starts must be >= 1.");

    if (numStarts > 100)
        throw std::invalid_argument("Number of starts must be <= 100.");

    this->numStarts = numStarts;
}

double Boss::update(OrderGraph& tau) {
    double score = tau.score;
    double bic = tau.bic;
    
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
	
	tau.parentMap[*it] = std::set<Node>(parents.begin(), parents.end());
	tau.scoreMap[*it] = parents.size();
	score += parents.size();
	tau.bicMap[*it] = localBic;
	bic += localBic;
    }
    tau.score = score;
    tau.bic = bic;
    return score;
}

double Boss::updateParallel(OrderGraph& tau) {

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

			  // {
			  //     std::lock_guard<std::mutex> gstLock(gstMutexMap[*it]);
			  parents = gstMap[*it]->search(prefix, &localBic);
			  // }


			  {
			      std::lock_guard<std::mutex> scoreLock(scoreMutex);
			      tau.parentMap[*it] = std::set<Node>(parents.begin(),
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

// Boss::OrderGraph Boss::bestMove(Node n, OrderGraph tau) {
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

// 	// RcppThread::Rcout << "      Old BIC = " << newTau.bic << "\n";

// 	update(newTau);

// 	// RcppThread::Rcout << "      New BIC = " << newTau.bic << "\n";
	
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

Boss::OrderGraph Boss::bestMove(Node n, OrderGraph tau) {
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

bool Boss::bossForward() {

    OrderGraph tau(pi);

    double oldBic = 1e20;
    // double oldScore = tau.order.size() * tau.order.size();
    double bic = tau.bic;
    // double score = tau.score;

    // if (verbose) Rcpp::Rcout << "  Running BOSS...\n";

    int iter = 0;

    if (verbose) Rcpp::Rcout << "\r      Iter:  " << iter << "    Edges = " << tau.score << ", BIC = " << tau.bic << "      ";

    while (oldBic - bic > 1e-6) {
	oldBic = bic;
	// oldScore = score;
	std::vector<Node> nodeList(tau.getShuffledNodes());
	for (Node n : nodeList) {
	    tau = bestMove(n, tau);
	    // if (verbose) Rcpp::Rcout << "\r    Score = " << tau.score << ", BIC = " << tau.bic;
	    RcppThread::checkUserInterrupt();
	}
	// score = tau.score;
	bic = tau.bic;
	iter++;

	if (verbose) Rcpp::Rcout << "\r      Iter:  " << iter << "    Edges = " << tau.score << ", BIC = " << tau.bic << "      ";
    }

    if (verbose) Rcpp::Rcout << "\n";

    bool improved = tau.bic < pi.bic + 1e-6;

    if (improved)
	pi = tau;
    
    return improved;
}

std::set<Node> Boss::getNaXY(const Node& x, const Node& y, EdgeListGraph& graph) {
    std::set<Node> naXY;
    std::vector<Node> adjy = graph.getAdjacentNodes(y);

    for (const Node& z : adjy) {
	Edge edge = graph.getEdge(y, z);
	if (!edge.isUndirected()) continue;
	if (!graph.isAdjacentTo(x,z)) continue;
	naXY.insert(z);
    }

    return naXY;
}

bool Boss::isClique(const std::set<Node>& nodes, EdgeListGraph& graph) {
    if (nodes.size()<=1) return true;
    
    std::vector<Node> nodeList(nodes.begin(), nodes.end());

    for (int i = 0; i < nodes.size()-1; i++) {
	for (int j = i + 1; j < nodes.size(); j++) {
	    if (!graph.isAdjacentTo(nodeList[i], nodeList[j])) {
		return false;
	    }
	}
    }

    // for (Node n1 : nodes) {
    // 	for (Node n2 : nodes) {
    // 	    if (n2 >= n1) break;
    // 	    if (!graph.isAdjacentTo(n1, n2)) {
    // 		return false;
    // 	    }
    // 	}
    // }
    
    return true;
}

bool Boss::validDelete(EdgeListGraph& graph, DeleteTask& task) {
    std::set<Node> complement(task.naXY.begin(), task.naXY.end());
    for (Node h : task.H) {
	complement.erase(h);
    }

    if (!isClique(complement, graph)) {
	return false;
    }

    return true;
}

void Boss::deleteProducer(EdgeListGraph& graph) {
    // int idx = 0;
    for (Edge edge : graph.getEdges()) {
	// idx++;
	// RcppThread::Rcout << "      Edge " << idx << ": " << edge << std::endl;
	Node x;
	Node y;
	std::set<Node> naXY, parents, H;
	std::vector<Node> naXYvec, paVec, Hvec;
	int depth;
	std::vector<int>* choice;
	if (edge.isDirected()) {
	    if (edge.getEndpoint1() == ENDPOINT_TAIL) {
		x = edge.getNode1();
		y = edge.getNode2();
	    } else {
		y = edge.getNode1();
		x = edge.getNode2();
	    }

	    naXY = getNaXY(x, y, graph);
	    naXYvec = std::vector<Node>(naXY.begin(), naXY.end());
	    paVec = graph.getParents(y);
	    parents = std::set<Node>(paVec.begin(), paVec.end());
	    parents.erase(x);

	    depth = naXY.size();
	    DepthChoiceGenerator cg(naXY.size(), depth);

	    for (choice = cg.next(); choice != NULL; choice = cg.next()) {
	        Hvec = GraphUtils::asList(*choice, naXYvec);
		H = std::set<Node>(Hvec.begin(), Hvec.end());

		DeleteTask task(x, y, naXY, H, parents);

		// if (validDelete(graph, task)) {
		// tasks.push_back(task);
		taskQueue.push(task);
		// }
	    }
	    
	} else {
	    x = edge.getNode1();
	    y = edge.getNode2();

	    naXY = getNaXY(x, y, graph);
	    naXYvec = std::vector<Node>(naXY.begin(), naXY.end());
	    paVec = graph.getParents(y);
	    parents = std::set<Node>(paVec.begin(), paVec.end());
	    parents.erase(x);

	    depth = naXY.size();
	    DepthChoiceGenerator cg(naXY.size(), depth);

	    for (choice = cg.next(); choice != NULL; choice = cg.next()) {
	        Hvec = GraphUtils::asList(*choice, naXYvec);
		H = std::set<Node>(Hvec.begin(), Hvec.end());

		DeleteTask task(x, y, naXY, H, parents);

		// if (validDelete(graph, task)) {
		// tasks.push_back(task);
		taskQueue.push(task);
		// }
	    }
	    
	    x = edge.getNode2();
	    y = edge.getNode1();

	    naXY = getNaXY(x, y, graph);
	    naXYvec = std::vector<Node>(naXY.begin(), naXY.end());
	    paVec = graph.getParents(y);
	    parents = std::set<Node>(paVec.begin(), paVec.end());
	    parents.erase(x);

	    depth = naXY.size();
	    cg = DepthChoiceGenerator(naXY.size(), depth);

	    for (choice = cg.next(); choice != NULL; choice = cg.next()) {
	        Hvec = GraphUtils::asList(*choice, naXYvec);
		H = std::set<Node>(Hvec.begin(), Hvec.end());

		DeleteTask task(x, y, naXY, H, parents);

		// if (validDelete(graph, task)) {
		// tasks.push_back(task);
		taskQueue.push(task);
		// }
	    }
	}
	if (RcppThread::isInterrupted()) {
	    break;
	}
    }

    DeleteTask poisonPill;
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(poisonPill);
    }
}

void Boss::deleteConsumer(EdgeListGraph& graph) {
    while(true) {
	
        DeleteTask task = taskQueue.pop();

        //If poison, return
        if (task.x.isNull() && task.y.isNull()) return;

	if (RcppThread::isInterrupted()) {
	    break;
	}

	std::set<Node> complement(task.naXY.begin(), task.naXY.end());
	for (Node h : task.H) {
	    complement.erase(h);
	}

	// RcppThread::Rcout << "    " << task.x << " --> " << task.y << " naXY size = " << task.naXY.size() << "\n";
	// RcppThread::Rcout << "    " << task.x << " --> " << task.y << " H size = " << task.H.size() << "\n";
	// RcppThread::Rcout << "    " << task.x << " --> " << task.y << " complement size = " << complement.size() << "\n";

	if (!isClique(complement, graph)) continue;

	// RcppThread::Rcout << "    " << task.x << " --> " << task.y << " is valid\n";

	// std::set<Node> complement(task.naXY.begin(), task.naXY.end());
	// for (Node h : task.H) {
	//     complement.erase(h);
	// }

	std::vector<Node> condSet(task.parents.begin(), task.parents.end());
	for (Node c : complement) {
	    condSet.push_back(c);
	}

	// RcppThread::Rcout << "    " << task.x << " --> " << task.y << " condSet size = " << condSet.size() << "\n";


	double newScore = scorer->localScore(task.y, condSet);

	condSet.push_back(task.x);

	double oldScore = scorer->localScore(task.y, condSet);

	double scoreChange = newScore - oldScore;

	// RcppThread::Rcout << "    " << task.x << " --> " << task.y << " : " << scoreChange << "\n";

	if (scoreChange < -1e-6) {
	    std::lock_guard<std::mutex> deleteLock(deleteMutex);
	    // RcppThread::Rcout << "      Deletion: " << task.x
	    // 		      << " --> " << task.y << " : " << scoreChange << std::endl;
	    if (deleteScores.count(task)==0) {
		deleteScores[task] = scoreChange;
	    } else if (scoreChange - deleteScores[task] < -1e-6) {
		deleteScores.erase(task);
		deleteScores[task] = scoreChange;
	    }
	}
    }
}


bool Boss::bestDelete(EdgeListGraph& graph) {

    deleteScores.clear();

    std::vector<RcppThread::Thread> threads;

    threads.push_back(RcppThread::Thread( [&] { deleteProducer(graph); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( [&] { deleteConsumer(graph); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
	if (threads[i].joinable()) {
	    threads[i].join();
	} else {
	    Rcpp::Rcout << "#### THREAD " << i << " NOT JOINABLE ####\n";
	}
    }

    if (deleteScores.size()==0)
	return false;

    DeleteTask bestTask;
    double bestScoreChange = 0.0;
    for (auto it = deleteScores.begin(); it != deleteScores.end(); it++) {
	if (it->second < bestScoreChange) {
	    // Rcpp::Rcout << "      Current Best Deletion: " << it->first.x
	    // 		<< " --> " << it->first.y << " : " << it->second << std::endl;
	    bestTask = it->first;
	    bestScoreChange = it->second;
	}
    }

    if (verbose) {
	Rcpp::Rcout << "      Deletion: " << bestTask.x
		    << " --> " << bestTask.y << " : Score = "
		    << graph.getScore() + bestScoreChange << std::endl;
    }
	    

    graph.removeEdge(bestTask.x, bestTask.y);

    for (Node h : bestTask.H) {
	graph.removeEdge(h, bestTask.y);
	graph.addDirectedEdge(bestTask.y, h);

	if (graph.isAdjacentTo(h, bestTask.x)) {
	    if (graph.getEdge(h, bestTask.x).isUndirected()) {
		graph.removeEdge(h, bestTask.x);
		graph.addDirectedEdge(bestTask.x, h);
	    }
	}
    }

    double oldScore = graph.getScore();
    
    graph = graph.getCPDAG();

    graph.setScore(oldScore + bestScoreChange);

    return true;
}

bool Boss::bossBES() {
    piGraph = pi.toGraph();
    double oldBic = pi.bic;
    piGraph.setScore(pi.bic);

    bool improved = true;

    while (improved) {
	improved = bestDelete(piGraph);
    }

    pi = OrderGraph(piGraph.getCausalOrdering());

    updateParallel(pi);

    if (verbose) Rcpp::Rcout << "      Edges = " << pi.score << ", Score = " << pi.bic << "\n";
    
    return (piGraph.getScore() - oldBic < -1e-6) || (pi.bic - oldBic < -1e-6); // std::abs(pi.bic - oldBic) > 1e-6; //  || (pi.score - piGraph.getNumEdges() < -0.5)
}


EdgeListGraph Boss::search() {

    OrderGraph bestOrder;

    for (int i = 0; i < numStarts; i++) {
	if (verbose) Rcpp::Rcout << "Run " << i+1 << "...\n";
	
	nodes = initializeRandom();
	pi = OrderGraph(nodes);
	score = updateParallel(pi);

	if (verbose) Rcpp::Rcout << "  Running BOSS...\n";

	bool improved = true;
	bool first = true;
	EdgeListGraph oldPiGraph = pi.toGraph();

	while (improved) {
	    if (verbose) Rcpp::Rcout << "    Forward Best Order Score Search...\n";
	    improved = bossForward();
	    if (!first && !improved)
		break;
	    if (verbose) Rcpp::Rcout << "    Backward Equivalence Search...\n";
	    improved = bossBES();
	    first = false;
	    if (piGraph == oldPiGraph)
		break;
	    oldPiGraph = piGraph;
	}

	// EdgeListGraph piGraph2 = pi.toGraph();

	// Rcpp::Rcout << "    pi order: { ";
	// for (Node n : pi.order) {
	//     Rcpp::Rcout << n << " ";
	// }
	// Rcpp::Rcout << "}\n    cpdag causal order 1: { ";
	// for (Node n : piGraph.getCausalOrdering()) {
	//     Rcpp::Rcout << n << " ";
	// }
	// Rcpp::Rcout << "}\n    cpdag causal order 2: { ";
	// for (Node n : piGraph2.getCausalOrdering()) {
	//     Rcpp::Rcout << n << " ";
	// }
	// Rcpp::Rcout << "}\n";	
		
	if (verbose) Rcpp::Rcout << "\n";
        
	if (piGraph.getScore() - graph.getScore() < -1e-6 || (piGraph.getScore() - graph.getScore() < 1e-6 && piGraph.getNumEdges() - graph.getNumEdges() < -0.5)) {
	    bestOrder = pi;
	    graph = piGraph;
	    
	}

	resetTrees();
    }

    // graph = bestOrder.toGraph();

    std::ostringstream alg;
    if (initialGraph==NULL) {
	alg << "BOSS";
    } else {
	alg << initialGraph->getAlgorithm() << "-" << "BOSS";
    }
    graph.setAlgorithm(alg.str());

    // graph.setScore(bestOrder.bic);

    graph.setHyperParam("penalty", arma::vec({ penalty }));
        
    return graph;
}
