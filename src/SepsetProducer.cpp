#include "SepsetProducer.hpp"
#include "GraphUtils.hpp"

void SepsetProducer::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

void SepsetProducer::producer() {
    for (auto&& b : graph.getNodes()) {
        
        std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);
	Node nullNode;
	std::vector<Node> possParents = possibleParents(b, adjacentNodes, nullNode);

        if (possParents.size() < 2) {
            continue;
        }

	// std::sort(possParents.begin(),
	// 	  possParents.end(),
	// 	  [] (const Node& a, const Node& b) {return a < b; }
	//     );

        ChoiceGenerator cg(possParents.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Node a = possParents[(*choice)[0]];
            Node c = possParents[(*choice)[1]];

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c) || (a.isCensored() && c.isCensored())) {
                continue;
            }

	    std::vector<Node> possibleDsep;
	    if (sepsets.isInSepsetMap(a, c))
		possibleDsep = sepsets.get(a, c);
	    
	    if (!possibleDsep.empty()) {
		taskQueue.push(IndependenceTask(a, b, c, possibleDsep));
	    }

	    std::vector<Node> adja = graph.getAdjacentNodes(a);
	    std::vector<Node> adjc = graph.getAdjacentNodes(c);

	    std::vector<Node> ppa = possibleParents(a, adja, c);
	    std::vector<Node> ppc = possibleParents(c, adjc, a);

	    DepthChoiceGenerator cg1(ppa.size(), depth);
	    std::vector<int> *comb2;
	    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb2, ppa);
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }

	    DepthChoiceGenerator cg2(ppc.size(), depth);
	    std::vector<int> *comb3;
	    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb3, ppc);
		if (s.empty()) continue;
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }
        }

	if (RcppThread::isInterrupted()) {
	    break;
	}
    }

    // Poison pill
    for (int i = 0; i < parallelism; i++) {
	taskQueue.push(IndependenceTask());
    }
}

void SepsetProducer::consumer() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a.isNull() || it.b.isNull() || it.c.isNull()) break;

	double score = 0.0;

        bool indep = test->isIndependent(it.a, it.c, it.s, &score);

	if (indep) {
	    
	    std::lock_guard<std::mutex> mapLock(mapMutex);
	    
	    Triple triple(it.a, it.b, it.c);

	    if (sepsetCount.count(triple) == 0) {
		sepsetCount[triple] = {0, 0};
		maxP[triple] = 0.0;
		maxPCollider[triple] = false;
	    }
		
	    if (std::find(it.s.begin(), it.s.end(), it.b) != it.s.end()) {
		std::pair<int, int> current = sepsetCount[triple];
		sepsetCount[triple] = {current.first + 1, current.second};
		if (score > maxP[triple]) {
		    maxP[triple] = score;
		    maxPCollider[triple] = false;
		}
	    } else {
		std::pair<int, int> current = sepsetCount[triple];
		sepsetCount[triple] = {current.first, current.second + 1};
		if (score > maxP[triple]) {
		    maxP[triple] = score;
		    maxPCollider[triple] = true;
		}
	    }
	}
    }

}


void SepsetProducer::producerSepsetMap() {
    for (auto&& b : graph.getNodes()) {
        
        std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);
	Node nullNode;
	std::vector<Node> possParents = possibleParents(b, adjacentNodes, nullNode);

        if (possParents.size() < 2) {
            continue;
        }

	// std::sort(possParents.begin(),
	// 	  possParents.end(),
	// 	  [] (const Node& a, const Node& b) {return a < b; }
	//     );

        ChoiceGenerator cg(possParents.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Node a = possParents[(*choice)[0]];
            Node c = possParents[(*choice)[1]];

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c) || (a.isCensored() && c.isCensored())) {
                continue;
            }

	    if (sepsets.isInSepsetMap(a, c))
		continue;
	    
	    std::vector<Node> adja = graph.getAdjacentNodes(a);
	    std::vector<Node> adjc = graph.getAdjacentNodes(c);

	    std::vector<Node> ppa = possibleParents(a, adja, c);
	    std::vector<Node> ppc = possibleParents(c, adjc, a);

	    DepthChoiceGenerator cg1(ppa.size(), depth);
	    std::vector<int> *comb2;
	    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb2, ppa);
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }

	    DepthChoiceGenerator cg2(ppc.size(), depth);
	    std::vector<int> *comb3;
	    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb3, ppc);
		if (s.empty()) continue;
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }
        }

	if (RcppThread::isInterrupted()) {
	    break;
	}
    }

    // Poison pill
    for (int i = 0; i < parallelism; i++) {
	taskQueue.push(IndependenceTask());
    }
}

void SepsetProducer::consumerSepsetMap() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a.isNull() || it.b.isNull() || it.c.isNull()) break;

	if (sepsets.isInSepsetMap(it.a, it.c))
	    continue;

	double score = 0.0;

        bool indep = test->isIndependent(it.a, it.c, it.s, &score);

	if (indep) {
	    std::lock_guard<std::mutex> mapLock(mapMutex);
	    sepsets.set(it.a, it.c, it.s, score);
	}
    }
}


void SepsetProducer::fillMap() {  
    if (!mapFilled) {
	sepsetCount.clear();
	maxP.clear();
	maxPCollider.clear();

	std::vector<RcppThread::Thread> threads;
    
	if (rule == ORIENT_SEPSETS) {
	    threads.push_back(RcppThread::Thread( [this] { producerSepsetMap(); } ));

	    for (int i = 0; i < parallelism; i++) {
		threads.push_back(RcppThread::Thread( [this] { consumerSepsetMap(); } ));
	    }

	    for (int i = 0; i < threads.size(); i++) {
		threads[i].join();
	    }
	} else {

	    threads.push_back(RcppThread::Thread( [this] { producer(); } ));

	    for (int i = 0; i < parallelism; i++) {
		threads.push_back(RcppThread::Thread( [this] { consumer(); } ));
	    }

	    for (int i = 0; i < threads.size(); i++) {
		threads[i].join();
	    }
	}
    }

    mapFilled = true;

    if (rule == ORIENT_MAJORITY || rule == ORIENT_CONSERVATIVE) {
	for (auto element : sepsetCount) {
	    Triple t = element.first;

	    if (!isCollider(t.x, t.y, t.z) && !isNoncollider(t.x, t.y, t.z)) {
		ambiguous.insert(t);
	    }
	}
    }

}

std::vector<Node> SepsetProducer::getSepset(const Node& a, const Node& b) {

    if (rule==ORIENT_SEPSETS) {
	if (sepsets.isInSepsetMap(a, b))
	    return sepsets.get(a, b);
    }
    
    std::vector<Node> sepset;
    double pval = 0;

    if (a.isCensored() && b.isCensored()) {
	sepsets.set(a, b, sepset, 1.0);
	return sepset;
    }
    
    std::vector<Node> possibleDsep;
    if (sepsets.isInSepsetMap(a, b))
	possibleDsep = sepsets.get(a, b);
    
    if (!possibleDsep.empty()) {
	double score;
	bool indep = test->isIndependent(a, b, possibleDsep, &score);
	pval = score;
	sepset = possibleDsep;
    }

    std::vector<Node> adja = graph.getAdjacentNodes(a);
    std::vector<Node> adjb = graph.getAdjacentNodes(b);

    std::vector<Node> ppa = possibleParents(a, adja, b);
    std::vector<Node> ppb = possibleParents(b, adjb, a);

    int maxDepth = std::max(ppa.size(), ppb.size());
    maxDepth = std::min(maxDepth, std::max(depth, 1000));

    for (int d = 0; d <= maxDepth; d++) {

	if (d <= ppa.size()) {
	    ChoiceGenerator cg1(ppa.size(), d);
	    std::vector<int> *comb2;
	    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb2, ppa);
		double score;
		bool indep = test->isIndependent(a, b, s, &score);
		if (indep) {
		    pval = score;
		    sepset = s;
		    // if (rule==ORIENT_SEPSETS) break;
		}
	    }
	}

	if (d <= ppb.size()) {
	    ChoiceGenerator cg2(ppb.size(), d);
	    std::vector<int> *comb3;
	    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb3, ppb);
		double score;
		bool indep = test->isIndependent(a, b, s, &score);
		if (indep && (score > pval)) {
		    pval = score;
		    sepset = s;
		    // if (rule==ORIENT_SEPSETS) break;
		}
	    }
	}

	if (rule==ORIENT_SEPSETS && pval > 0) break;
    }

    sepsets.set(a, b, sepset, pval);
    
    return sepset;
}


bool SepsetProducer::isCollider(const Node& i, const Node& j, const Node& k) {
    if (rule == ORIENT_SEPSETS) {
	std::vector<Node> sepset;
	if (sepsets.isInSepsetMap(i, k)) {
	    sepset = sepsets.get(i, k);
	    
	} else {
	    sepset = getSepset(i, k);
	}

	// Rcpp::Rcout << "   <" << i << ", " << k << "> : { ";
	// for (auto n : sepset) {
	//     Rcpp::Rcout << n << " ";
	// }
	// Rcpp::Rcout << "}\n";
	
	
	return std::count(sepset.begin(), sepset.end(), j) == 0;
    }
    
    if (sepsetCount.count(Triple(i,j,k)) == 0) {
    	// std::vector<Node> sepset = getSepset(i,k);
    	// return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
	return false;
    }
    
    if (rule == ORIENT_MAJORITY || rule == ORIENT_CONSERVATIVE) {
	std::pair<int, int> count = sepsetCount[Triple(i,j,k)];
	return count.second > 0 && count.second > count.first &&
	    (rule == ORIENT_MAJORITY || count.first == 0);
    } else if (rule == ORIENT_MAXP) {
	return maxPCollider[Triple(i,j,k)];
    } else {
	throw std::runtime_error("Orientation rule is not set");
    }
    return false;
}

bool SepsetProducer::isNoncollider(const Node& i, const Node& j, const Node& k) {
    if (rule == ORIENT_SEPSETS) {
        std::vector<Node> sepset;
	if (sepsets.isInSepsetMap(i, k)) {
	    sepset = sepsets.get(i, k);
	} else {
	    sepset = getSepset(i, k);
	}
	
	return (std::count(sepset.begin(), sepset.end(), j) >= 1);
    }
    
    if (sepsetCount.count(Triple(i,j,k)) == 0) {
    	// std::vector<Node> sepset = getSepset(i,k);
    	// return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
	return false;
    }
    
    if (rule == ORIENT_MAJORITY || rule == ORIENT_CONSERVATIVE) {
	std::pair<int, int> count = sepsetCount[Triple(i,j,k)];
	return count.first > 0 && count.first > count.second &&
	    (rule == ORIENT_MAJORITY || count.second == 0);
    } else if (rule == ORIENT_MAXP) {
	return maxPCollider[Triple(i,j,k)];
    } else {
	throw std::runtime_error("Orientation rule is not set");
    }
    return false;
}

std::vector<Node> SepsetProducer::possibleParents(const Node& x,
						  const std::vector<Node>& adjx,
						  const Node& y) {
    std::vector<Node> possParents;

    for (const Node& z : adjx) {
	if (z==x) continue;
	if (z==y) continue;

	if (possibleParentOf(x, z)) {
	    possParents.push_back(z);
	}
    }

    return possParents;
}

bool SepsetProducer::possibleParentOf(const Node& x, const Node& z) {
    return !knowledge.isForbidden(z, x) && !knowledge.isRequired(x, z);
}


std::vector<Triple> SepsetProducer::getOrderedColliders() {
    std::vector<Triple> orderedColliders;
    std::map<Triple, double> score;

    // if (!mapFilled) {
    // 	throw std::runtime_error("fillMap not called");
    // }
    
    for (auto&& b : graph.getNodes()) {
        std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);
	Node nullNode;
	std::vector<Node> possParents = possibleParents(b, adjacentNodes, nullNode);

        if (possParents.size() < 2) {
            continue;
        }

	// std::sort(possParents.begin(),
	// 	  possParents.end(),
	// 	  [] (const Node& a, const Node& b) {return a < b; }
	//     );

        ChoiceGenerator cg(possParents.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Node a = possParents[(*choice)[0]];
            Node c = possParents[(*choice)[1]];

	    // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c) || (a.isCensored() && c.isCensored())) {
                continue;
            }

	    // Rcpp::Rcout << Triple(a,b,c) << std::endl;

	    // Rcpp::Rcout << "{ " << a << ", " << c << " }: p value = " << sepsets.getPValue(a, c) << std::endl;

	    if (isCollider(a, b, c)) {
		Triple t(a,b,c);
		// Rcpp::Rcout << t << " is a collider" << std::endl;
		orderedColliders.push_back(t);
		score[t] = 0.0;
		if (rule == ORIENT_SEPSETS) {
		    // Rcpp::Rcout << "{ " << a << ", " << c << " }: p value = " << sepsets.getPValue(a, c) << std::endl;
		    score[t] = sepsets.getPValue(a, c);
		} else if (rule == ORIENT_MAXP) {
		    // Rcpp::Rcout << "Max P: " << maxP[t] << std::endl;
		    // Rcpp::Rcout << "{ " << a << ", " << c << " }: p value = " << maxP[t] << std::endl;
		    score[t] = maxP[t];
		} else if (rule == ORIENT_MAJORITY) {
		    std::pair<int, int> count = sepsetCount[t];
		    // Rcpp::Rcout << "Count: " << count.first << ", " << count.second << std::endl;
		    // Rcpp::Rcout << "Ratio: " << count.second / (double) std::max(count.first + count.second, 1) << std::endl;
		    score[t] = count.second / (double) std::max(count.first + count.second, 1);
		    if (score[t] == 1)
			score[t] = (double) count.second;
		    
		} else if (rule == ORIENT_CONSERVATIVE) {
		    std::pair<int, int> count = sepsetCount[t];
		    // Rcpp::Rcout << "Count: " << count.first << ", " << count.second << std::endl;
		    score[t] = (double) count.second;
		} else {
		    throw std::invalid_argument("Invalid orientation rule detected");
		}
		// Rcpp::Rcout << t << ": " << score[t] << std::endl;
	    }
	}
    }

    // Rcpp::Rcout << "Score map filled:\n";

    std::sort(orderedColliders.begin(), orderedColliders.end(),
	      [&](const Triple& t1, const Triple& t2) {
		  if (score[t1] == score[t2])
		      return t1 < t2;
		  return score[t1] < score[t2];
	      });

    // for (auto t : orderedColliders) {
    // 	Rcpp::Rcout << "  " << t << ": " << score[t] << std::endl;
    // }

    return orderedColliders;
}
