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

	std::sort(adjacentNodes.begin(),
		  adjacentNodes.end(),
		  [] (Variable* a, Variable* b) {return a->getName() < b->getName(); }
	    );

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Variable* a = adjacentNodes[(*choice)[0]];
            Variable* c = adjacentNodes[(*choice)[1]];
	    // if (adjacentNodes[(*choice)[0]]->getName()
	    // 	< adjacentNodes[(*choice)[1]]->getName()) {
	    // 	a = adjacentNodes[(*choice)[0]];
	    // 	c = adjacentNodes[(*choice)[1]];
	    // } else {
	    // 	a = adjacentNodes[(*choice)[1]];
	    // 	c = adjacentNodes[(*choice)[0]];
	    // }

            // Skip triples that are shielded.
            if (graph.isAdjacentTo(a, c)) {
                continue;
            }

	    std::vector<Variable*> possibleDsep = *extraSepsets.get(a, c);
	    if (!possibleDsep.empty()) {
		taskQueue.push(IndependenceTask(a, b, c, possibleDsep));
	    }

	    std::vector<Variable*> adja = graph.getAdjacentNodes(a);
	    std::vector<Variable*> adjc = graph.getAdjacentNodes(c);

	    // {
	    // 	std::unique_lock<std::mutex> mapLock(mapMutex);
	    //     mapCondition.wait(mapLock, [this] { return !mapModifying; });
	    // 	mapModifying = true;
	    // 	scores[Triple(a, b, c)] = 0;
	    // 	colliders[Triple(a, b, c)] = 0;
	    // 	mapModifying = false;
	    // }
	    // mapCondition.notify_one();

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

void SepsetProducerMaxP::consumer() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a == NULL && it.b == NULL && it.c == NULL) break;

        double score;
        bool indep = test->isIndependent(it.a, it.c, it.s, &score);

	{
	    std::unique_lock<std::mutex> mapLock(mapMutex);
	    mapCondition.wait(mapLock,
			      [this] {
				  if (!mapModifying) {
				      return mapModifying = true;
				  }
				  return false;
			      });
	    // mapModifying = true;
	    // if (verbose) {
	    // 	Rcpp::Rcout << it.a->getName() << " _||_ " << it.c->getName() << " | [";
	    // 	for (Variable* n : it.s) Rcpp::Rcout << n->getName() << ",";
	    // 	Rcpp::Rcout << "]\n";
	    // }
	    if (scores.count(Triple(it.a, it.b, it.c))) {
		if (indep && (score > scores[Triple(it.a, it.b, it.c)])) {
		    scores.at(Triple(it.a, it.b, it.c)) = score;
		    colliders.at(Triple(it.a, it.b, it.c)) = (std::find(it.s.begin(), it.s.end(), it.b) == it.s.end());
		}
	    } else {
		if (indep) {
		    scores[Triple(it.a, it.b, it.c)] = score;
		    colliders[Triple(it.a, it.b, it.c)] = (std::find(it.s.begin(), it.s.end(), it.b) == it.s.end());
		}
		else {
		    scores[Triple(it.a, it.b, it.c)] = 0;
		    colliders[Triple(it.a, it.b, it.c)] = 0;
		}
	    }
	    // if (score > scores[Triple(it.a, it.b, it.c)]) {
	    // 	scores[Triple(it.a, it.b, it.c)] = score;
	    // 	colliders[Triple(it.a, it.b, it.c)] = (std::find(it.s.begin(), it.s.end(), it.b) == it.s.end());
	    // }
	    mapModifying = false;
	}
	mapCondition.notify_one();
    }

}

void SepsetProducerMaxP::fillMap() {

    if (verbose) Rcpp::Rcout << "Filling Triple map...\n";

    colliders.clear();
    scores.clear();

    mapModifying = false;

    // for (Variable* b : graph.getNodes()) {
        
    //     std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(b);

    //     if (adjacentNodes.size() < 2) {
    //         continue;
    //     }

    //     ChoiceGenerator cg(adjacentNodes.size(), 2);
    //     std::vector<int> *choice;

    //     for (choice = cg.next(); choice != NULL; choice = cg.next()) {
    //         Variable* a = adjacentNodes[(*choice)[0]];
    //         Variable* c = adjacentNodes[(*choice)[1]];

    //         // Skip triples that are shielded.
    //         if (graph.isAdjacentTo(a, c)) {
    //             continue;
    //         }

    // 	    std::vector<Variable*> adja = graph.getAdjacentNodes(a);
    // 	    std::vector<Variable*> adjc = graph.getAdjacentNodes(c);

    // 	    // {
    // 	    // 	std::unique_lock<std::mutex> mapLock(mapMutex);
    // 	    //     mapCondition.wait(mapLock, [this] { return !mapModifying; });
    // 	    // 	mapModifying = true;
    // 	    // 	scores[Triple(a, b, c)] = 0;
    // 	    // 	colliders[Triple(a, b, c)] = 0;
    // 	    // 	mapModifying = false;
    // 	    // }
    // 	    // mapCondition.notify_one();

    // 	    DepthChoiceGenerator cg1(adja.size(), depth);
    // 	    std::vector<int> *comb2;
    // 	    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
    // 		std::vector<Variable*> s = GraphUtils::asList(*comb2, adja);

    // 		double score;
    // 		bool indep = test->isIndependent(a, c, s, &score);
        
    // 		// if (verbose) {
    // 		// 	Rcpp::Rcout << a->getName() << " _||_ " << c->getName() << " | [";
    // 		// 	for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
    // 		// 	Rcpp::Rcout << "]\n";
    // 		// }
    // 		if (scores.count(Triple(a, b, c))) {
    // 		    if (indep && (score > scores[Triple(a, b, c)])) {
    // 			scores.at(Triple(a, b, c)) = score;
    // 			colliders.at(Triple(a, b, c)) = (std::find(s.begin(), s.end(), b) == s.end());
    // 		    }
    // 		} else {
    // 		    if (indep) {
    // 			scores[Triple(a, b, c)] = score;
    // 			colliders[Triple(a, b, c)] = (std::find(s.begin(), s.end(), b) == s.end());
    // 		    }
    // 		    else {
    // 			scores[Triple(a, b, c)] = 0;
    // 			colliders[Triple(a, b, c)] = 0;
    // 		    }
    // 		}
    // 	    }

    // 	    DepthChoiceGenerator cg2(adjc.size(), depth);
    // 	    std::vector<int> *comb3;
    // 	    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
    // 		std::vector<Variable*> s = GraphUtils::asList(*comb3, adjc);

    // 		double score;
    // 		bool indep = test->isIndependent(a, c, s, &score);
        
    // 		// if (verbose) {
    // 		// 	Rcpp::Rcout << a->getName() << " _||_ " << c->getName() << " | [";
    // 		// 	for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
    // 		// 	Rcpp::Rcout << "]\n";
    // 		// }
    // 		if (scores.count(Triple(a, b, c))) {
    // 		    if (indep && (score > scores[Triple(a, b, c)])) {
    // 			scores.at(Triple(a, b, c)) = score;
    // 			colliders.at(Triple(a, b, c)) = (std::find(s.begin(), s.end(), b) == s.end());
    // 		    }
    // 		} else {
    // 		    if (indep) {
    // 			scores[Triple(a, b, c)] = score;
    // 			colliders[Triple(a, b, c)] = (std::find(s.begin(), s.end(), b) == s.end());
    // 		    }
    // 		    else {
    // 			scores[Triple(a, b, c)] = 0;
    // 			colliders[Triple(a, b, c)] = 0;
    // 		    }
    // 		}
    // 	    }
    //     }
    // }

    std::vector<std::thread> threads;

    threads.push_back(std::thread( [this] { producer(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( [this] { consumer(); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

}

std::vector<Variable*>  SepsetProducerMaxP::getSepset(Variable* a, Variable* b) {
    std::vector<Variable*> sepset;
    double pval = 0;

    // if (extraSepsets != null) {
    std::vector<Variable*> possibleDsep = *extraSepsets.get(a, b);
    if (!possibleDsep.empty()) {
	double score;
	bool indep = test->isIndependent(a, b, possibleDsep, &score);
	pval = score;
	sepset = possibleDsep;
    }
	// }

    std::vector<Variable*> adja = graph.getAdjacentNodes(a);
    std::vector<Variable*> adjb = graph.getAdjacentNodes(b);

    DepthChoiceGenerator cg1(adja.size(), depth);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
	std::vector<Variable*> s = GraphUtils::asList(*comb2, adja);
	// if (verbose) {
	//     Rcpp::Rcout << a->getName() << " ? " << b->getName() << " [";
	//     for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
	//     Rcpp::Rcout << "]" << std::endl;
	// }
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep && (score > pval)) {
	    // if (verbose) {
	    // 	Rcpp::Rcout << a->getName() << " _||_ " << b->getName() << " [";
	    // 	for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
	    // 	Rcpp::Rcout << "]\t" << score << "\t" << indep << std::endl;
	    // }
	    pval = score;
	    sepset = s;
	}
    }

    DepthChoiceGenerator cg2(adjb.size(), depth);
    std::vector<int> *comb3;
    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
	std::vector<Variable*> s = GraphUtils::asList(*comb3, adjb);
	// if (verbose) {
	//     Rcpp::Rcout << a->getName() << " ? " << b->getName() << " [";
	//     for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
	//     Rcpp::Rcout << "]" << std::endl;
	// }
	double score;
	bool indep = test->isIndependent(a, b, s, &score);
	if (indep && (score > pval)) {
	    // if (verbose) {
	    // 	Rcpp::Rcout << a->getName() << " _||_ " << b->getName() << " [";
	    // 	for (Variable* n : s) Rcpp::Rcout << n->getName() << ",";
	    // 	Rcpp::Rcout << "] " << score << std::endl;
	    // }
	    pval = score;
	    sepset = s;
	}
    }
    
    return sepset;
}

bool SepsetProducerMaxP::isCollider(Variable* i, Variable* j, Variable* k) {
    // std::vector<Variable*> sepset = getSepset(i,k);
    // return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
    // if (verbose) Rcpp::Rcout << "Checking " << Triple(i,j,k) << "\n";
    if (colliders.count(Triple(i,j,k)) == 0) {
	// if (verbose) Rcpp::Rcout << "getSepset called\n";
	std::vector<Variable*> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) == sepset.end());
    }
    // if (verbose) Rcpp::Rcout << Triple(i,j,k) << " mapped: " << colliders[Triple(i,j,k)] << "\n";
    return colliders[Triple(i,j,k)];
}

bool SepsetProducerMaxP::isNoncollider(Variable* i, Variable* j, Variable* k) {
    // std::vector<Variable*> sepset = getSepset(i,k);
    // return (std::find(sepset.begin(), sepset.end(), j) != sepset.end());
    if (colliders.count(Triple(i,j,k)) == 0) {
	// if (verbose) Rcpp::Rcout << "getSepset called\n";
	std::vector<Variable*> sepset = getSepset(i,k);
	return (std::find(sepset.begin(), sepset.end(), j) != sepset.end());
    }
    return !colliders[Triple(i,j,k)];
}
