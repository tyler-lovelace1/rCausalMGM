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
    for (Variable* b : graph.getNodes()) {
        
        std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            continue;
        }

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Variable* a = adjacentNodes[(*choice)[0]];
            Variable* c = adjacentNodes[(*choice)[1]];

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c)) {
                continue;
            }

	    std::vector<Variable*> adja = graph.getAdjacentNodes(a);
	    std::vector<Variable*> adjc = graph.getAdjacentNodes(c);

	    {
		std::unique_lock<std::mutex> mapLock(mapMutex);
		scores[Triple(a, b, c)] = 0;
	    }

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
		taskQueue.push(IndependenceTask(a, b, c, s));
	    }
        }
    }

    // Poison pill
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(IndependenceTask(NULL, NULL, NULL, {}));
    }
}

void SepsetProducerMaxP::consumer() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a == NULL && it.b == NULL && it.c == NULL) break;

        double score;
        bool indep = test->isIndependent(it.a, it.c, it.s, &score);

	if (indep) {
	    std::unique_lock<std::mutex> mapLock(mapMutex);
	    // if (verbose) {
	    // 	Rcpp::Rcout << it.a->getName() << " _||_ " << it.c->getName() << " | [";
	    // 	for (Variable* n : it.s) Rcpp::Rcout << n->getName() << ",";
	    // 	Rcpp::Rcout << "]\n";
	    // }
	    if (score > scores[Triple(it.a, it.b, it.c)]) {
		scores[Triple(it.a, it.b, it.c)] = score;
		colliders[Triple(it.a, it.b, it.c)] = (std::find(it.s.begin(), it.s.end(), it.b) == it.s.end());
	    }
	}
        
    }

}

void SepsetProducerMaxP::fillMap() {

    colliders.clear();
    scores.clear();

    std::vector<std::thread> threads;

    threads.push_back(std::thread( [this] { producer(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( [this] { consumer(); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    if (verbose) Rcpp::Rcout << "Map filled\n";

}

std::vector<Variable*>  SepsetProducerMaxP::getSepset(Variable* a, Variable* b) {
    std::vector<Variable*> sepset;
    double pval = 0;

    std::vector<Variable*> adja = graph.getAdjacentNodes(a);
    std::vector<Variable*> adjb = graph.getAdjacentNodes(b);

    DepthChoiceGenerator cg1(adja.size(), depth);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
	std::vector<Variable*> s = GraphUtils::asList(*comb2, adja);
	if (verbose) {
	    Rcpp::Rcout << a->getName() << " ? " << b->getName() << " [";
	    for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
	    Rcpp::Rcout << "]" << std::endl;
	}
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep && (score > pval)) {
	    if (verbose) {
		Rcpp::Rcout << a->getName() << " _||_ " << b->getName() << " [";
		for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
		Rcpp::Rcout << "]\t" << score << "\t" << indep << std::endl;
	    }
	    pval = score;
	    sepset = s;
	}
    }

    DepthChoiceGenerator cg2(adjb.size(), depth);
    std::vector<int> *comb3;
    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
	std::vector<Variable*> s = GraphUtils::asList(*comb3, adjb);
	if (verbose) {
	    Rcpp::Rcout << a->getName() << " ? " << b->getName() << " [";
	    for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
	    Rcpp::Rcout << "]" << std::endl;
	}
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep && (score > pval)) {
	    if (verbose) {
		Rcpp::Rcout << a->getName() << " _||_ " << b->getName() << " [";
		for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
		Rcpp::Rcout << "] " << score << std::endl;
	    }
	    pval = score;
	    sepset = s;
	}
    }
    
    return sepset;
}

bool SepsetProducerMaxP::isCollider(Variable* i, Variable* j, Variable* k) {
    if (verbose) Rcpp::Rcout << "Checking " << Triple(i,j,k) << "\n";
    if (colliders.find(Triple(i,j,k)) == colliders.end()) {
	if (verbose) Rcpp::Rcout << "getSepset called\n";
	std::vector<Variable*> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
    }
    if (verbose) Rcpp::Rcout << Triple(i,j,k) << " mapped: " << colliders[Triple(i,j,k)] << "\n";
    return colliders[Triple(i,j,k)];
}

bool SepsetProducerMaxP::isNoncollider(Variable* i, Variable* j, Variable* k) {
    if (colliders.find(Triple(i,j,k)) == colliders.end()) {
	if (verbose) Rcpp::Rcout << "getSepset called\n";
	std::vector<Variable*> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) != sepset.end());
    }
    return !colliders[Triple(i,j,k)];
}
