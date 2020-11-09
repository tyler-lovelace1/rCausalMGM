#include "OrientCollidersMaxP.hpp"

#include <fstream>
#include "GraphUtils.hpp"

OrientCollidersMaxP::OrientCollidersMaxP(IndependenceTest *test, EdgeListGraph *graph) : taskQueue(MAX_QUEUE_SIZE) {
    if (test == NULL) 
        throw std::invalid_argument("independenceTest may not be NULL.");

    this->independenceTest = test;

    if (graph == NULL) 
        throw std::invalid_argument("graph may not be NULL.");

    this->graph = graph;

    Rcpp::Rcout << "OCMP paralellism = " << parallelism << std::endl;
    if (parallelism == 0) {
        parallelism = 4;
        Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
    }
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

    // Poison pill
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(IndependenceTask(NULL, NULL, NULL, {}));
    }
}

void OrientCollidersMaxP::consumer() {
    std::ofstream logfile;
    logfile.open("../test_results/orientCollidersMaxP.log", std::ios_base::app);

    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a == NULL && it.b == NULL && it.c == NULL) break;

        double score;
        independenceTest->isIndependent(it.a, it.c, it.s, &score);

        logfile << "Testing " << it.a->getName() << " || " << it.c->getName() << " _||_ {";
        for (Variable* n : it.s) logfile << " " << n->getName() << " ";
        logfile << "} with p = " << score << std::endl;

        std::unique_lock<std::mutex> mapLock(mapMutex);
        if (score > scores[Triple(it.a, it.b, it.c)]) {
            scores[Triple(it.a, it.b, it.c)] = score;
            colliders[Triple(it.a, it.b, it.c)] = (std::find(it.s.begin(), it.s.end(), it.b) == it.s.end());
        }
        
    }

    logfile.close();
}

void OrientCollidersMaxP::testColliderMaxP(Variable* a, Variable* b, Variable* c) {
    
    std::vector<Variable*> adja = graph->getAdjacentNodes(a);
    std::vector<Variable*> adjc = graph->getAdjacentNodes(c);

    std::unique_lock<std::mutex> mapLock(mapMutex);
    scores[Triple(a, b, c)] = 0;
    mapMutex.unlock();

    DepthChoiceGenerator cg1(adja.size(), -1);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
        std::vector<Variable*> s = GraphUtils::asList(*comb2, adja);

        // If !s.contains(b) TODO - should this condition be here?
        taskQueue.push(IndependenceTask(a, b, c, s));
    }

    DepthChoiceGenerator cg2(adjc.size(), -1);
    std::vector<int> *comb3;
    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
        std::vector<Variable*> s = GraphUtils::asList(*comb3, adjc);
        // If !s.contains(b) TODO - should this condition be here?
        taskQueue.push(IndependenceTask(c, b, a, s));
    }

}

void OrientCollidersMaxP::testColliderHeuristic(Variable* a, Variable* b, Variable* c) {
    // TODO - should this be here?
    if (graph->getEdges(a, b).size() > 1 || graph->getEdges(b, c).size() > 1) {
        return;
    }

    std::unique_lock<std::mutex> mapLock(mapMutex);
    scores[Triple(a, b, c)] = 0;
    mapMutex.unlock();

    taskQueue.push(IndependenceTask(a, b, c, {}));
    taskQueue.push(IndependenceTask(a, b, c, {b}));
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
    // Rcpp::Rcout << "ORIENTED SUCCESSFULLY" << std::endl;
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

//TODO - what should this be testing?
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
        // logfile << "t = " << t->getName() << std::endl;

        // if (e == NULL) logfile << "e = NULL" << std::endl;
        // else logfile << "e = " << e->getName() << std::endl;

        if (e == t) {
            e = NULL;
            distance++;
            // logfile << "distance = " << distance << " bound = " << bound << std::endl;
            if (distance > (bound == -1 ? 1000 : bound)) {
                logfile << "No" << std::endl;
                logfile.close();
                return false;
            }
        }

        for (Variable* u : graph->getAdjacentNodes(t)) {
            // logfile << "u = " << u->getName() << std::endl;
            if (u == z && distance > 2) {
                logfile << "Yes" << std::endl;
                logfile.close();
                return true;
            } 

            if (v.count(u) == 0) {
                v.insert(u);
                q.push(u);

                if (e == NULL) {
                    e = u;
                }
            }
        }
    }

    logfile << "No" << std::endl;
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
    std::vector<Triple> colliderList;
    colliderList.reserve(scores.size());
    for (auto pair : scores) {
        // We only care about triples that are colliders
        if (colliders[pair.first]) {
            colliderList.push_back(pair.first);

            // Print scores
            Rcpp::Rcout << pair.first << " score = " << pair.second << std::endl;
        }
    }

    // Most independent ones first.
    std::sort(colliderList.begin(), colliderList.end(),
        [&](const Triple& t1, const Triple& t2) {
            return scores[t1] > scores[t2];
        });

    for (Triple triple : colliderList) {
        Variable* a = triple.getX();
        Variable* b = triple.getY();
        Variable* c = triple.getZ();

        if (!(graph->getEndpoint(b, a) == ENDPOINT_ARROW || graph->getEndpoint(b, c) == ENDPOINT_ARROW)) {
            // Rcpp::Rcout << "orienting collider " << Triple(a, b, c) << std::endl;
            orientCollider(a, b, c);
        }
    }

    
}
