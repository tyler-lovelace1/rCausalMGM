#include "OrientCollidersMaxP.hpp"

#include <fstream>

OrientCollidersMaxP::OrientCollidersMaxP(IndependenceTest *test, EdgeListGraph *graph) : taskQueue(MAX_QUEUE_SIZE) {
    if (test == NULL) 
        throw std::invalid_argument("independenceTest may not be NULL.");

    this->independenceTest = test;

    if (graph == NULL) 
        throw std::invalid_argument("graph may not be NULL.");

    this->graph = graph;
}

void OrientCollidersMaxP::setDepth(int depth) {
    if (depth < -1)
        throw std::invalid_argument("Depth must be -1 (unlimited) or >= 0.");

    if (depth > 1000)
        throw std::invalid_argument("Depth must be <= 1000.");

    this->depth = depth;
}

void OrientCollidersMaxP::producer() {
    for (Variable* b : graph->getNodes()) {
        
        std::vector<Variable*> adjacentNodes = graph->getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            continue;
        }

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
            Variable* a = adjacentNodes[(*choice)[0]];
            Variable* c = adjacentNodes[(*choice)[1]];

            // Skip triples that are shielded.
            if (graph->isAdjacentTo(a, c)) {
                continue;
            }

            taskQueue.push(Triple(a, b, c));
        }

    }

    // Poison pill
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(Triple(NULL, NULL, NULL));
    }
}

void OrientCollidersMaxP::consumer() {
    while(true) {
        Triple triple = taskQueue.pop();
        Variable* a = triple.getX();
        Variable* b = triple.getY();
        Variable* c = triple.getZ();

        // Poison Pill
        if (a == NULL && b == NULL && c == NULL) break;

        if (useHeuristic) {
            if (existsShortPath(a, c, maxPathLength)) {
                testColliderMaxP(a, b, c);
            } else {
                testColliderHeuristic(a, b, c);
            }
        } else {
            testColliderMaxP(a, b, c);
        }
    }
}

void OrientCollidersMaxP::testColliderMaxP(Variable* a, Variable* b, Variable* c) {
    std::ofstream logfile;
    logfile.open("../test_results/orientCollidersMaxP.log", std::ios_base::app);

    std::vector<Variable*> adja = graph->getAdjacentNodes(a);
    std::vector<Variable*> adjc = graph->getAdjacentNodes(c);

    double score = std::numeric_limits<double>::infinity();
    std::vector<Variable*> S;

    DepthChoiceGenerator cg1(adja.size(), -1);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
        std::vector<Variable*> s = GraphUtils::asList(*comb2, adja);

        double _score;
        independenceTest->isIndependent(c, a, s, &_score);

        logfile << "Testing " << c->getName() << " || " << a->getName() << " _||_ {";
        for (Variable* n : s) logfile << " " << n->getName() << " ";
        logfile << "} with p = " << _score << std::endl;

        if (_score < score) {
            score = _score;
            S = s;
        }
    }

    DepthChoiceGenerator cg2(adjc.size(), -1);
    std::vector<int> *comb3;
    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
        std::vector<Variable*> s = GraphUtils::asList(*comb3, adjc);

        double _score;
        independenceTest->isIndependent(c, a, s, &_score);

        logfile << "Testing " << c->getName() << " || " << a->getName() << " _||_ {";
        for (Variable* n : s) logfile << " " << n->getName() << " ";
        logfile << "} with p = " << _score << std::endl;

        if (_score < score) {
            score = _score;
            S = s;
        }
    }

    // If !S.contains(b)
    if (std::find(S.begin(), S.end(), b) == S.end()) {
        std::unique_lock<std::mutex> mapLock(mapMutex);
        scores[Triple(a, b, c)] = score;
    }

    logfile.close();

}

void OrientCollidersMaxP::testColliderHeuristic(Variable* a, Variable* b, Variable* c) {
    std::ofstream logfile;
    logfile.open("../test_results/orientCollidersMaxP.log", std::ios_base::app);
    logfile << "DUDEK USING HEURISTIC " << Triple(a, b, c) << std::endl;

    std::vector<Variable*> empty = {};
    double s1;
    independenceTest->isIndependent(a, c, empty, &s1);

    std::vector<Variable*> bList = {b};
    double s2;
    independenceTest->isIndependent(a, c, bList, &s2);

    bool mycollider2 = s2 > s1;

    // Skip triples that are shielded.
    if (graph->isAdjacentTo(a, c)) return;

    if (graph->getEdges(a, b).size() > 1 || graph->getEdges(b, c).size() > 1) {
        return;
    }

    if (mycollider2) {
        std::unique_lock<std::mutex> mapLock(mapMutex);
        scores[Triple(a, b, c)] = std::abs(s2);
    }

    logfile.close();
}

void OrientCollidersMaxP::orientCollider(Variable* a, Variable* b, Variable* c) {
    if (wouldCreateBadCollider(a, b)) return;
    if (wouldCreateBadCollider(c, b)) return;
    if (graph->getEdges(a, b).size() > 1) return;
    if (graph->getEdges(b, c).size() > 1) return;
    graph->removeEdge(a, b);
    graph->removeEdge(c, b);
    graph->addDirectedEdge(a, b);
    graph->addDirectedEdge(c, b);
    Rcpp::Rcout << "ORIENTED SUCCESSFULLY" << std::endl;
}

bool OrientCollidersMaxP::wouldCreateBadCollider(Variable* x, Variable* y) {
    std::unordered_set<Variable*> empty = {};
    std::unordered_set<Variable*> ySet = {y};

    for (Variable* z : graph->getAdjacentNodes(y)) {
        if (x == z) continue;

        if (!graph->isAdjacentTo(x, z) &&
                graph->getEndpoint(z, y) == ENDPOINT_ARROW &&
                !sepset(x, z, empty, ySet)) {
                    return true;
                }
    }

    return false;
}

// Finds a sepset containing the nodes in 'containing' but not the nodes in 'notContaining'
// Returns true of that set exists, false if not
// If the set exists, it will be placed into output
// If the set is not needed, set output = NULL
bool OrientCollidersMaxP::sepset(Variable* a, Variable* c, std::unordered_set<Variable*>& containing, std::unordered_set<Variable*>& notContaining, std::vector<Variable*>* output /* = NULL */) {
    std::vector<Variable*> adj = graph->getAdjacentNodes(a);
    for (Variable* node : graph->getAdjacentNodes(c)) adj.push_back(node);
    adj.erase(std::remove(adj.begin(), adj.end(), a), adj.end());  // adj.remove(a)
    adj.erase(std::remove(adj.begin(), adj.end(), c), adj.end());  // adj.remove(c)

    int endD = std::min((depth == -1 ? 1000 : depth), (int) adj.size());
    for (int d = 0; d <= endD; d++) {
        if (d <= adj.size()) {
            ChoiceGenerator cg(adj.size(), d);
            std::vector<int> *choice;

            for (choice = cg.next(); choice != NULL; choice = cg.next()) {
                std::unordered_set<Variable*> v2 = GraphUtils::asSet(*choice, adj);
                for (Variable* n : containing) v2.insert(n);
                for (Variable* n : notContaining) v2.erase(n);
                v2.erase(a);
                v2.erase(c);

                std::vector<Variable*> v2List(v2.begin(), v2.end());
                double p2;
                independenceTest->isIndependent(a, c, v2List, &p2);

                if (p2 < 0) {
                    if (output != NULL) {
                        *output = v2List;
                    }
                    return true;
                }

            }
        }
    }

    return false;
}

// Returns true if there is an undirected path from x to either y or z within the given number of steps.
bool OrientCollidersMaxP::existsShortPath(Variable* x, Variable* z, int bound) {
    std::ofstream logfile;
    logfile.open("../test_results/orientCollidersMaxP.log", std::ios_base::app);
    logfile << "Testing if a path exists between " << x->getName() << " and " << z->getName() << std::endl;
    
    std::queue<Variable*> q;
    std::unordered_set<Variable*> v;
    q.push(x);
    v.insert(x);
    Variable* e = NULL;
    int distance = 0;

    while(!q.empty()) {
        Variable* t = q.front();
        q.pop();
        logfile << "t = " << t->getName() << std::endl;

        if (e == NULL) logfile << "e = NULL" << std::endl;
        else logfile << "e = " << e->getName() << std::endl;

        if (e == t) {
            e = NULL;
            distance++;
            logfile << "distance = " << distance << " bound = " << bound << std::endl;
            if (distance > (bound == -1 ? 1000 : bound)) {
                logfile.close();
                return false;
            }
        }

        for (Variable* u : graph->getAdjacentNodes(t)) {
            logfile << "u = " << u->getName() << std::endl;
            Edge edge = graph->getEdge(t, u);
            Variable* c = Edge::traverse(t, edge);
            if (c == NULL) continue;
            logfile << "c = " << c->getName() << std::endl;

            if (c == z && distance > 2) {
                logfile.close();
                return true;
            } 

            if (v.count(c) == 0) {
                v.insert(c);
                q.push(c);

                if (e == NULL) {
                    e = u;
                }
            }
        }
    }

    logfile.close();

    return false;
}

void OrientCollidersMaxP::addColliders() {
    std::ofstream logfile;

    logfile.open("../test_results/orientCollidersMaxP.log");

    logfile << "Starting colliders" << std::endl;

    logfile.close();

    scores.clear();

    std::vector<std::thread> threads;

    threads.push_back(std::thread( [this] { producer(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(std::thread( [this] { consumer(); } ));
    }

    for (int i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    // Get a list of keys in scores
    std::vector<Triple> tripleList;
    tripleList.reserve(scores.size());
    for (auto pair : scores) {
        tripleList.push_back(pair.first);

        // Print scores
        Rcpp::Rcout << pair.first << " score = " << pair.second << std::endl;
    }

    // Most independent ones first.
    std::sort(tripleList.begin(), tripleList.end(),
        [&](const Triple& t1, const Triple& t2) {
            return scores[t1] > scores[t2];
        });

    for (Triple triple : tripleList) {
        Variable* a = triple.getX();
        Variable* b = triple.getY();
        Variable* c = triple.getZ();

        if (!(graph->getEndpoint(b, a) == ENDPOINT_ARROW || graph->getEndpoint(b, c) == ENDPOINT_ARROW)) {
            Rcpp::Rcout << "orienting collider " << Triple(a, b, c) << std::endl;
            orientCollider(a, b, c);
        }
    }

    
}
