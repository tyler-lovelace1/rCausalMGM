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
    for (Variable* b : graph.getNodes()) {
        
        std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            continue;
        }

	std::sort(adjacentNodes.begin(),
		  adjacentNodes.end(),
		  [] (Variable* a, Variable* b) {return a->getName() < b->getName(); }
	    );

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Variable* a = adjacentNodes[(*choice)[0]];
            Variable* c = adjacentNodes[(*choice)[1]];

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c)) {
                continue;
            }

	    std::vector<Variable*> possibleDsep = *extraSepsets.get(a, c);
	    if (!possibleDsep.empty()) {
		taskQueue.push(IndependenceTask(a, b, c, possibleDsep));
	    }

	    // std::vector<Variable*> adja = graph.getAdjacentNodes(a);
	    // std::vector<Variable*> adjc = graph.getAdjacentNodes(c);

	    std::vector<Variable*> adja = graph.getPossibleParents(a);
	    std::vector<Variable*> adjc = graph.getPossibleParents(c);

	    DepthChoiceGenerator cg1(adja.size(), depth);
	    std::vector<int> *comb2;
	    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
		std::vector<Variable*> s = GraphUtils::asList(*comb2, adja);
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }

	    DepthChoiceGenerator cg2(adjc.size(), depth);
	    std::vector<int> *comb3;
	    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
		std::vector<Variable*> s = GraphUtils::asList(*comb3, adjc);
		if (s.empty()) continue;
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }
        }
    }

    // Poison pill
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(IndependenceTask(NULL, NULL, NULL, {}));
    }
}

void SepsetProducerConservative::consumer() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a == NULL && it.b == NULL && it.c == NULL) break;

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

    std::vector<std::thread> threads;

    threads.push_back(std::thread( [this] { producer(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( [this] { consumer(); } ));
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

std::vector<Variable*>  SepsetProducerConservative::getSepset(Variable* a, Variable* b) {
    std::vector<Variable*> sepset;
    double pval = 0;

    std::vector<Variable*> possibleDsep = *extraSepsets.get(a, b);
    if (!possibleDsep.empty()) {
	double score;
	bool indep = test->isIndependent(a, b, possibleDsep, &score);
	pval = score;
	sepset = possibleDsep;
    }

    std::vector<Variable*> adja = graph.getAdjacentNodes(a);
    std::vector<Variable*> adjb = graph.getAdjacentNodes(b);

    DepthChoiceGenerator cg1(adja.size(), depth);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
	std::vector<Variable*> s = GraphUtils::asList(*comb2, adja);
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
	std::vector<Variable*> s = GraphUtils::asList(*comb3, adjb);
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep && (score > pval)) {
	    pval = score;
	    sepset = s;
	}
    }
    
    return sepset;
}

bool SepsetProducerConservative::isCollider(Variable* i, Variable* j, Variable* k) {
    if (sepsetCount.count(Triple(i,j,k)) == 0) {
	std::vector<Variable*> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
    }
    std::pair<int, int> count = sepsetCount[Triple(i,j,k)];
    return (count.first == 0) && (count.second > 0);
}

bool SepsetProducerConservative::isNoncollider(Variable* i, Variable* j, Variable* k) {
    if (sepsetCount.count(Triple(i,j,k)) == 0) {
	std::vector<Variable*> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) != sepset.end());
    }
    std::pair<int, int> count = sepsetCount[Triple(i,j,k)];
    return (count.first > 0) && (count.second == 0);
}
