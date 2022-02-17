#include "SepsetProducerMaxP.hpp"
#include "GraphUtils.hpp"

void SepsetProducerMaxP::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

void SepsetProducerMaxP::producer() {
    for (auto&& b : graph.getNodes()) {
        
        std::vector<Node> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            continue;
        }

	std::sort(adjacentNodes.begin(),
		  adjacentNodes.end(),
		  [] (const Node& a, const Node& b) {return a < b; }
	    );

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Node a = adjacentNodes[(*choice)[0]];
            Node c = adjacentNodes[(*choice)[1]];

            if (graph.isAdjacentTo(a, c)) {
                continue;
            }

	    std::vector<Node> possibleDsep = *extraSepsets.get(a, c);
	    if (!possibleDsep.empty()) {
		taskQueue.push(IndependenceTask(a, b, c, possibleDsep));
	    }

	    std::vector<Node> adja = graph.getAdjacentNodes(a);
	    std::vector<Node> adjc = graph.getAdjacentNodes(c);

	    DepthChoiceGenerator cg1(adja.size(), depth);
	    std::vector<int> *comb2;
	    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb2, adja);
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }

	    DepthChoiceGenerator cg2(adjc.size(), depth);
	    std::vector<int> *comb3;
	    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
		std::vector<Node> s = GraphUtils::asList(*comb3, adjc);
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

void SepsetProducerMaxP::consumer() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
	if (it.a.isNull() || it.b.isNull() || it.c.isNull()) break;
	
        double score;
        bool indep = test->isIndependent(it.a, it.c, it.s, &score);

	{
	    std::lock_guard<std::mutex> mapLock(mapMutex);
	    
	    if (scores.count(Triple(it.a, it.b, it.c))) {
		if (indep && (score > scores[Triple(it.a, it.b, it.c)])) {
		    scores.at(Triple(it.a, it.b, it.c)) = score;
		    colliders.at(Triple(it.a, it.b, it.c)) = (std::find(it.s.begin(),
									it.s.end(),
									it.b) == it.s.end());
		}
	    } else {
		if (indep) {
		    scores[Triple(it.a, it.b, it.c)] = score;
		    colliders[Triple(it.a, it.b, it.c)] = (std::find(it.s.begin(),
								     it.s.end(),
								     it.b) == it.s.end());
		}
		else {
		    scores[Triple(it.a, it.b, it.c)] = 0;
		    colliders[Triple(it.a, it.b, it.c)] = 0;
		}
	    }
	}
    }

}

void SepsetProducerMaxP::fillMap() {

    colliders.clear();
    scores.clear();

    std::vector<RcppThread::Thread> threads;

    threads.push_back(RcppThread::Thread( [this] { producer(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( [this] { consumer(); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

}

std::vector<Node>  SepsetProducerMaxP::getSepset(const Node& a, const Node& b) {
    std::vector<Node> sepset;
    double pval = 0;

    // if (extraSepsets != null) {
    std::vector<Node> possibleDsep = *extraSepsets.get(a, b);
    if (!possibleDsep.empty()) {
	double score;
	bool indep = test->isIndependent(a, b, possibleDsep, &score);
	pval = score;
	sepset = possibleDsep;
    }
	// }

    std::vector<Node> adja = graph.getAdjacentNodes(a);
    std::vector<Node> adjb = graph.getAdjacentNodes(b);

    DepthChoiceGenerator cg1(adja.size(), depth);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
	std::vector<Node> s = GraphUtils::asList(*comb2, adja);
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep && (score > pval)) {
	    pval = score;
	    sepset = s;
	}
    }

    DepthChoiceGenerator cg2(adjb.size(), depth);
    std::vector<int> *comb3;
    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
	std::vector<Node> s = GraphUtils::asList(*comb3, adjb);
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep && (score > pval)) {
	    pval = score;
	    sepset = s;
	}
    }
    
    return sepset;
}

bool SepsetProducerMaxP::isCollider(const Node& i, const Node& j, const Node& k) {
    if (colliders.count(Triple(i,j,k)) == 0) {
	std::vector<Node> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
    }
    return colliders[Triple(i,j,k)];
}

bool SepsetProducerMaxP::isNoncollider(const Node& i, const Node& j, const Node& k) {
    if (colliders.count(Triple(i,j,k)) == 0) {
	std::vector<Node> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) != sepset.end());
    }
    return !colliders[Triple(i,j,k)];
}
