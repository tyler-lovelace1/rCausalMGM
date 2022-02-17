#include "SepsetProducerConservative.hpp"
#include "GraphUtils.hpp"

void SepsetProducerConservative::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

void SepsetProducerConservative::producer() {
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

            // Skip triples that are shielded.
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

void SepsetProducerConservative::consumer() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a.isNull() || it.b.isNull() || it.c.isNull()) break;

        bool indep = test->isIndependent(it.a, it.c, it.s);

	{
	    std::lock_guard<std::mutex> mapLock(mapMutex);

	    if (sepsetCount.count(Triple(it.a, it.b, it.c)) == 0) {
		sepsetCount[Triple(it.a, it.b, it.c)] = {0, 0};
	    }
	    
	    if (indep) {
		if (std::find(it.s.begin(), it.s.end(), it.b) != it.s.end()) {
		    std::pair<int, int> current = sepsetCount[Triple(it.a, it.b, it.c)];
		    sepsetCount[Triple(it.a, it.b, it.c)] = {current.first + 1, current.second};
		} else {
		    std::pair<int, int> current = sepsetCount[Triple(it.a, it.b, it.c)];
		    sepsetCount[Triple(it.a, it.b, it.c)] = {current.first, current.second + 1};
		}
	    }
	}
    }

}

void SepsetProducerConservative::fillMap() {

    sepsetCount.clear();

    std::vector<RcppThread::Thread> threads;

    threads.push_back(RcppThread::Thread( [this] { producer(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( [this] { consumer(); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    for (auto element : sepsetCount) {
        Triple t = element.first;

        if (!isCollider(t.x, t.y, t.z) && !isNoncollider(t.x, t.y, t.z)) {
            ambiguous.insert(t);
        }
    }

}

std::vector<Node>  SepsetProducerConservative::getSepset(const Node& a, const Node& b) {
    std::vector<Node> sepset;
    double pval = 0;

    std::vector<Node> possibleDsep = *extraSepsets.get(a, b);
    if (!possibleDsep.empty()) {
	double score;
	bool indep = test->isIndependent(a, b, possibleDsep, &score);
	pval = score;
	sepset = possibleDsep;
    }

    std::vector<Node> adja = graph.getAdjacentNodes(a);
    std::vector<Node> adjb = graph.getAdjacentNodes(b);

    DepthChoiceGenerator cg1(adja.size(), depth);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
	std::vector<Node> s = GraphUtils::asList(*comb2, adja);
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep) {
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

bool SepsetProducerConservative::isCollider(const Node& i, const Node& j, const Node& k) {
    if (sepsetCount.count(Triple(i,j,k)) == 0) {
	std::vector<Node> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
    }
    std::pair<int, int> count = sepsetCount[Triple(i,j,k)];
    return (count.first == 0) && (count.second > 0);
}

bool SepsetProducerConservative::isNoncollider(const Node& i, const Node& j, const Node& k) {
    if (sepsetCount.count(Triple(i,j,k)) == 0) {
	std::vector<Node> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) != sepset.end());
    }
    std::pair<int, int> count = sepsetCount[Triple(i,j,k)];
    return (count.first > 0) && (count.second == 0);
}
