#include "OrientCollidersMaxP.hpp"

#include <fstream>
#include "GraphUtils.hpp"

OrientCollidersMaxP::OrientCollidersMaxP(IndependenceTest *test, EdgeListGraph *graph, int threads) : taskQueue(MAX_QUEUE_SIZE) {
    if (test == NULL) 
        throw std::invalid_argument("independenceTest may not be NULL.");

    this->independenceTest = test;

    if (graph == NULL) 
        throw std::invalid_argument("graph may not be NULL.");

    this->graph = graph;

    if (threads > 0) parallelism = threads;
    else {
        parallelism = std::thread::hardware_concurrency();
        if (parallelism == 0) {
            parallelism = 4;
            Rcpp::Rcout << "Couldn't detect number of processors. Defaulting to 4" << std::endl;
        }
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
    for (const Node& b : graph->getNodes()) {
        
        std::vector<Node> adjacentNodes = graph->getAdjacentNodes(b);

        if (adjacentNodes.size() < 2) {
            continue;
        }

	std::sort(adjacentNodes.begin(), adjacentNodes.end());

        ChoiceGenerator cg(adjacentNodes.size(), 2);
        std::vector<int> *choice;

        for (choice = cg.next(); choice != NULL; choice = cg.next()) {
	    const Node& a = adjacentNodes[(*choice)[0]];
            const Node& c = adjacentNodes[(*choice)[1]];

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

	if (RcppThread::isInterrupted()) {
	    break;
	}
    }

    // Poison pill
    for (int i = 0; i < parallelism; i++) {
        taskQueue.push(IndependenceTask());
    }
}

void OrientCollidersMaxP::consumer() {
    while(true) {
        IndependenceTask it = taskQueue.pop();

        // Poison Pill
        if (it.a.isNull() && it.b.isNull() && it.c.isNull()) break;

        double score = independenceTest->getAlpha() + 1e-5;
        bool indep = independenceTest->isIndependent(it.a, it.c, it.s, &score);


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

void OrientCollidersMaxP::testColliderMaxP(const Node& a, const Node& b, const Node& c) {
    
    std::vector<Node> adja = graph->getAdjacentNodes(a);
    std::vector<Node> adjc = graph->getAdjacentNodes(c);

    DepthChoiceGenerator cg1(adja.size(), -1);
    std::vector<int> *comb2;
    for (comb2 = cg1.next(); comb2 != NULL; comb2 = cg1.next()) {
        std::vector<Node> s = GraphUtils::asList(*comb2, adja);

        // If !s.contains(b) TODO - should this condition be here?
        taskQueue.push(IndependenceTask(a, b, c, s));
    }

    DepthChoiceGenerator cg2(adjc.size(), -1);
    std::vector<int> *comb3;
    for (comb3 = cg2.next(); comb3 != NULL; comb3 = cg2.next()) {
        std::vector<Node> s = GraphUtils::asList(*comb3, adjc);
        // If !s.contains(b) TODO - should this condition be here?
	if (s.empty()) continue;
        taskQueue.push(IndependenceTask(c, b, a, s));
    }

}

void OrientCollidersMaxP::testColliderHeuristic(const Node& a, const Node& b, const Node& c) {
    // TODO - should this be here?
    if (graph->getEdges(a, b).size() > 1 || graph->getEdges(b, c).size() > 1) {
        return;
    }

    taskQueue.push(IndependenceTask(a, b, c, {}));
    taskQueue.push(IndependenceTask(a, b, c, {b}));
}

void OrientCollidersMaxP::orientCollider(const Node& a, const Node& b, const Node& c) {
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

bool OrientCollidersMaxP::wouldCreateBadCollider(const Node& x, const Node& y) {
    std::unordered_set<Node> empty = {};
    std::unordered_set<Node> ySet = {y};

    for (const Node& z : graph->getAdjacentNodes(y)) {
        if (x == z) continue;

        // if (!graph->isAdjacentTo(x, z) &&
	//     graph->getEndpoint(z, y) == ENDPOINT_ARROW &&
	//     !sepset(x, z, empty, ySet)) {
	//     return true;
	// }
	if (!graph->isAdjacentTo(x, z) &&
	    graph->getEndpoint(z, y) == ENDPOINT_ARROW &&
	    !colliders.at(Triple(x, y, z))) {
	    return true;
	}
    }

    return false;
}

// Finds a sepset containing the nodes in 'containing' but not the nodes in 'notContaining'
// Returns true of that set exists, false if not
// If the set exists, it will be placed into output
// If the set is not needed, set output = NULL
bool OrientCollidersMaxP::sepset(const Node& a, const Node& c, std::unordered_set<Node>& containing, std::unordered_set<Node>& notContaining, std::vector<Node>* output /* = NULL */) {
    std::vector<Node> adj = graph->getAdjacentNodes(a);
    for (const Node& node : graph->getAdjacentNodes(c)) adj.push_back(node);
    adj.erase(std::remove(adj.begin(), adj.end(), a), adj.end());  // adj.remove(a)
    adj.erase(std::remove(adj.begin(), adj.end(), c), adj.end());  // adj.remove(c)

    int endD = std::min((depth == -1 ? 1000 : depth), (int) adj.size());
    for (int d = 0; d <= endD; d++) {
        if (d <= adj.size()) {
            ChoiceGenerator cg(adj.size(), d);
            std::vector<int> *choice;

            for (choice = cg.next(); choice != NULL; choice = cg.next()) {
                std::unordered_set<Node> v2 = GraphUtils::asSet(*choice, adj);
                for (const Node& n : containing) v2.insert(n);
                for (const Node& n : notContaining) v2.erase(n);
                v2.erase(a);
                v2.erase(c);

                std::vector<Node> v2List(v2.begin(), v2.end());
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
bool OrientCollidersMaxP::existsShortPath(const Node& x, const Node& z, int bound) {
    
    std::queue<Node> q;
    std::unordered_set<Node> v;
    q.push(x);
    v.insert(x);
    Node e = Node();
    int distance = 0;

    while(!q.empty()) {
        Node t = q.front();
        q.pop();
        
        if (e == t) {
	    e = Node();
            distance++;
            if (distance > (bound == -1 ? 1000 : bound)) {
                return false;
            }
        }

        for (const Node& u : graph->getAdjacentNodes(t)) {
            if (u == z && distance > 2) {
                return true;
            } 

            if (v.count(u) == 0) {
                v.insert(u);
                q.push(u);

                if (e.isNull()) {
                    e = u;
                }
            }
        }
    }

    return false;
}

void OrientCollidersMaxP::addColliders() {

    scores.clear();

    std::vector<RcppThread::Thread> threads;

    threads.push_back(RcppThread::Thread( [this] { producer(); } ));

    for (int i = 0; i < parallelism; i++) {
        threads.push_back(RcppThread::Thread( [this] { consumer(); } ));
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
        }
    }

    // Most independent ones first.
    std::sort(colliderList.begin(), colliderList.end(),
	      [&](const Triple& t1, const Triple& t2) {
		  return scores[t1] > scores[t2];
	      });

    for (Triple triple : colliderList) {
        const Node& a = triple.getX();
        const Node& b = triple.getY();
        const Node& c = triple.getZ();

        if (!(graph->getEndpoint(b, a) == ENDPOINT_ARROW || graph->getEndpoint(b, c) == ENDPOINT_ARROW)) {
            orientCollider(a, b, c);
        }
    }

    
}
