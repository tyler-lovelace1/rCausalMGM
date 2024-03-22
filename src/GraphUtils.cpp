// [[Rcpp::depends(BH)]]
#include "GraphUtils.hpp"

std::vector<std::string> GraphUtils::splitString(std::string s, const std::string& delim) {
    std::vector<std::string> tokens;

    std::size_t find;

    while ((find = s.find(delim)) != std::string::npos) {
        std::string token = s.substr(0, find);
        tokens.push_back(token);
        s.erase(0, find + delim.length());
    }

    tokens.push_back(s);

    return tokens;
}

/**
 * Constructs a list of nodes from the given <code>nodes</code> list at the
 * given indices in that list.
 *
 * @param indices The indices of the desired nodes in <code>nodes</code>.
 * @param nodes The list of nodes from which we select a sublist.
 * @return the The sublist selected.
 */
std::vector<Node> GraphUtils::asList(std::vector<int>& indices, std::vector<Node>& nodes) {
    std::vector<Node> list;

    for (int i : indices) {
        list.push_back(nodes[i]);
    }

    return list;
}

std::unordered_set<Node> GraphUtils::asSet(std::vector<int>& indices, std::vector<Node>& nodes) {
    std::unordered_set<Node> set;

    for (int i : indices) {
        set.insert(nodes[i]);
    }

    return set;
}

EdgeListGraph GraphUtils::completeGraph(EdgeListGraph& graph) {
    EdgeListGraph graph2(graph.getNodes());

    graph2.removeEdges();

    std::vector<Node> nodes = graph2.getNodes();

    for (int i = 0; i < nodes.size(); i++) {
        for (int j = i+1; j < nodes.size(); j++) {
            Node node1 = nodes[i];
            Node node2 = nodes[j];
            graph2.addUndirectedEdge(node1, node2);
        }
    }

    return graph2;
}

EdgeListGraph GraphUtils::undirectedGraph(EdgeListGraph& graph) {
    EdgeListGraph graph2(graph.getNodes());

    for (Edge edge : graph.getEdges()) {
        graph2.addUndirectedEdge(edge.getNode1(), edge.getNode2());
    }

    return graph2;
}

std::set<Node> GraphUtils::possibleDsep2(const Node& x, const Node& y, EdgeListGraph& graph, int maxPathLength) {
    std::set<Node> dsep;

    // RcppThread::Rcout << "  Checking for Possible-D-Sep(" << x << "," << y << "):\n";

    std::queue<NodePair> Q;
    std::set<NodePair> V;

    std::map<Node, std::set<Node>> previous;
    previous[x] = {};

    int distance = 0;

    NodePair firstPairAtDepth = std::pair<Node,Node>(Node(),Node());

    for (Node b : graph.getAdjacentNodes(x)) {
	if (b == y) continue;

	NodePair edge = std::pair<Node, Node>(x, b);

	if (firstPairAtDepth.first.isNull() && firstPairAtDepth.second.isNull()) {
	    firstPairAtDepth = edge;
	}

	Q.push(edge);
	V.insert(edge);
	
	dsep.insert(b);

	previous[x].insert(b);
    }

    while (!Q.empty()) {
	NodePair edge = Q.front();
	Q.pop();

	if (firstPairAtDepth == edge) {
	    firstPairAtDepth = std::pair<Node,Node>(Node(),Node());
	    distance++;
	    if (distance > 0 && distance > (maxPathLength == -1 ? 1000 : maxPathLength)) {
		break;
	    }
	}

	Node a = edge.first;
	Node b = edge.second;

	// RcppThread::Rcout << "    For b = " << b << " previous = { ";
	// if (previous.count(b)>0) {
	//     for (Node n : previous[b]) {
	// 	RcppThread::Rcout << n << " ";
	//     }
	// }
	// RcppThread::Rcout << "}\n";

	// RcppThread::Rcout << "    Searching for semidirected paths from the above nodes to " << b << " or " << x << std::endl;

	if (existOnePathWithPossibleParents(previous, graph, b, x, b)) {
	    // RcppThread::Rcout << "      Paths found\n";
	    dsep.insert(b);
	}
	
	// dsep.insert(b);

	for (Node c : graph.getAdjacentNodes(b)) {
	    if (c == a) continue;
	    if (c == x) continue;
	    if (c == y) continue;

	    if (previous.count(c)==0) {
		previous[c] = { b };
	    } else {
		previous[c].insert(b);
	    }

	    if (graph.isDefCollider(a,b,c) || graph.isAdjacentTo(a,c)) {
		NodePair newEdge = std::pair<Node, Node>(a, c);
		
		if (V.count(newEdge)!=0) continue;

		V.insert(newEdge);
		Q.push(newEdge);

		if (firstPairAtDepth.first.isNull() && firstPairAtDepth.second.isNull()) {
		    firstPairAtDepth = newEdge;
		}
	    }
	}
    }

    // RcppThread::Rcout << "Poss-D-Sep(" << x << "," << y << ")\n";

    // for (auto it = previous.begin(); it != previous.end(); it++) {
    // 	RcppThread::Rcout << "    " << it->first << "  :  { ";
    // 	for (Node n : it->second) {
    // 	    RcppThread::Rcout << n << " ";
    // 	}
    // 	RcppThread::Rcout << "}\n";
    // }

    dsep.erase(x);
    dsep.erase(y);

    return dsep;
}

bool GraphUtils::existOnePathWithPossibleParents(std::map<Node,std::set<Node>>& previous,
						 EdgeListGraph& graph, const Node& w,
						 const Node& x, const Node& b) {
    if (w == x) return true;

    if (previous.count(w)==0) return false;
    
    for (Node r : previous[w]) {
	if (r == b || r == x) continue;

	if (graph.existsSemiDirectedPathFromTo(r, x) ||
	    graph.existsSemiDirectedPathFromTo(r, b)) {
	    
	    return true;
	}
    }
    
    return false;
}

std::unordered_set<Node> GraphUtils::possibleDsep(Node x, Node y, EdgeListGraph& graph, int maxPathLength) {
    std::unordered_set<Node> dsep;

    std::queue<boost::optional<std::pair<Node,Node>>> Q;
    std::unordered_set<std::pair<Node,Node>, boost::hash<std::pair<Node, Node> > > V;

    std::unordered_map<Node, std::vector<Node>> previous;
    std::vector<Node> null_vector = {};
    previous.insert(std::pair<Node, std::vector<Node>>(x, null_vector));

    boost::optional<std::pair<Node,Node>> e = {};
    int distance = 0;

    for (Node b : graph.getAdjacentNodes(x)) {
        if (b == y) {
            continue;
        }
        boost::optional<std::pair<Node,Node> > edge = std::pair<Node,Node>(x, b);
        if (!e) {
            e = edge;
        }
        Q.push(edge);
        V.insert(*edge);
        addToList(previous, b, x);
        dsep.insert(b);
    }
    while (!Q.empty()) {
        boost::optional<std::pair<Node, Node> > t = Q.front();
        Q.pop();

        if (e == t) {
            e = {};
            distance++;
            if (distance > 0 && distance > (maxPathLength == -1 ? 1000 : maxPathLength)) {
                break;
            }
        }

        Node a = t->first;
        Node b = t->second;
        if (existOnePathWithPossibleParents(previous, b, x, b, graph)) {
            dsep.insert(b);
        }

        for (Node c : graph.getAdjacentNodes(b)) {
            if (c == a) {
                continue;
            }
            if (c == x) {
                continue;
            }
            if (c == y) {
                continue;
            }

            addToList(previous, b, c);

            if (graph.isDefCollider(a, b, c) || graph.isAdjacentTo(a, c)) {
                boost::optional<std::pair<Node, Node> > u = std::pair<Node, Node>(a,c);
                if (std::count(V.begin(), V.end(), *u) != 0) {
                    continue;
                }

                V.insert(*u);
                Q.push(u);

                if (!e) {
                    e = u;
                }
            }
        }
    }
    dsep.erase(x);
    dsep.erase(y);
    return dsep;
}

void GraphUtils::addToList(std::unordered_map<Node, std::vector<Node>> previous, Node b, Node c) {
    std::vector<Node> list = previous[c];

    if (list.empty()) {
        std::vector<Node> null_list = {};
        list = null_list;
    }

    list.push_back(b);
}

bool GraphUtils::existOnePathWithPossibleParents(std::unordered_map<Node, std::vector<Node>> previous, Node w, Node x, Node b, EdgeListGraph& graph) {
    if (w == x) {
        return true;
    }
    const std::vector<Node> p = previous[w];
    if (p.empty()) {
        return false;
    }

    for (Node r : p) {
        if (r == b || r == x) {
            continue;
        }

        if ((existsSemidirectedPath(r, x, graph))
                || existsSemidirectedPath(r, b, graph)) {
            if (existOnePathWithPossibleParents(previous, r, x, b, graph)) {
                return true;
            }
        }
    }
    return false;
}

bool GraphUtils::existsSemidirectedPath(Node from, Node to, EdgeListGraph& G) {
    std::queue<Node> Q;

    std::unordered_set<Node> V;
    Q.push(from);
    V.insert(from);

    while (!Q.empty()) {
        Node t = Q.front();
        Q.pop();
        if (t == to) {
            return true;
        }
        for (Node u : G.getAdjacentNodes(t)) {
            Edge edge = G.getEdge(t, u);
            Node c = Edge::traverseSemiDirected(t, edge);

            if (c.isNull()) {
                continue;
            }
            if (std::count(V.begin(), V.end(), c) != 0) {
                continue;
            }

            V.insert(c);
            Q.push(c);
        }
    }
    return false;
}

bool GraphUtils::existsPossibleColliderPath(Node from, Node to, std::unordered_set<Node> Z, EdgeListGraph& G) {
    std::queue<Node> Q;

    if (Z.count(to) > 0) Z.erase(to);

    std::unordered_set<Node> V;
    Q.push(from);
    V.insert(from);

    while (!Q.empty()) {
        Node a = Q.front();
        Q.pop();
        if (a == to) {
            return true;
        }
        for (Node b : G.getAdjacentNodes(a)) {
	    
	    if (V.count(b) > 0) {
		continue;
	    }
	    
	    if (Z.count(b) > 0) {
		for (Node c : G.getAdjacentNodes(b)) {
		    if (G.isDefCollider(a,b,c)) {
			if (c == to) {
			    return true;
			}
			Edge edge = G.getEdge(a, b);
			if (edge.isBidirected()) {
			    V.insert(b);
			    Q.push(b);
			}
		    }
		}
	    } else {
		Edge edge = G.getEdge(a, b);
		if (a == from || edge.isBidirected() || Edge::traverseSemiDirected(b, edge) == a) {
		    if (b == to) {
			return true;
		    }
		}
	    }
        }
    }
    return false;
}

std::vector<Node> GraphUtils::getSepset(const Node& x, const Node& y, EdgeListGraph& graph) {
    if (!x.isObserved())
	throw std::runtime_error(x.getName() + " is an unobserved variable");
    if (!y.isObserved())
	throw std::runtime_error(y.getName() + " is an unobserved variable");
    
    std::vector<Node> sepset = { Node() };
    if (x == y) return sepset;

    std::vector<Node> z, oldZ;
    std::set<Node> path;
    std::set<Triple> colliders;
    bool first = true;

    while (first || z != oldZ) {
	first = false;
	oldZ = z;

	path = { x };
	colliders = {};

	for (const Node& b : graph.getAdjacentNodes(x)) {
	    if (sepsetPathFound(x, b, y, path, z, colliders, graph)) {
		return sepset;
	    }
	}

	// Rcpp::Rcout << "z for " << x << " and " << y << ":\n  ";
	// for (const Node& node : z) {
	//     Rcpp::Rcout << node.getName() << " ";
	// }
	// Rcpp::Rcout << std::endl;
    }

    // // std::unordered_set<Triple> colliders;
    // std::list<Node> path;
    // std::vector<Node> z, oldZ;
    // std::set<std::pair<std::list<Node>,std::vector<Node>>> visited;
    // std::stack<Node> nodeStack;
    // std::stack<std::list<Node>> pathStack;
    // std::stack<std::vector<Node>> sepsetStack;
    // // bool first = true;

    // // while (z != oldZ || first) {
    // // 	oldZ = z;
    // // 	first = false;
    // // 	visited.clear();

    // for (Node b : graph.getAdjacentNodes(x)) {
    // 	path = {x};
    // 	nodeStack.push(b);
    // 	pathStack.push(path);
    // 	sepsetStack.push(z);
    // 	visited.insert(std::pair<std::list<Node>,std::vector<Node>>(path, z));
    // }
    
    // while (!nodeStack.empty()) {
    // 	Node b = nodeStack.top();
    // 	nodeStack.pop();
    // 	path = pathStack.top();
    // 	pathStack.pop();
    // 	z = sepsetStack.top();
    // 	sepsetStack.pop();
    // 	Node a = path.back();

    // 	// Rcpp::Rcout << "nodeStack size = " << nodeStack.size()
    // 	// 	    << "\nNode " << b << std::endl;

    // 	if (std::find(path.begin(), path.end(), b) != path.end()) continue;

    // 	path.push_back(b);

    // 	if (visited.count(std::pair<std::list<Node>,std::vector<Node>>(path, z)) > 0)
    // 	    continue;

    // 	visited.insert(std::pair<std::list<Node>,std::vector<Node>>(path, z));

    // 	if (b == y) {
    // 	    Rcpp::Rcout << "Node y reached" << std::endl;
    // 	    break;
    // 	}

    // 	std::vector<Node> passNodes = getPassNodes(a, b, z, graph);

    // 	// Rcpp::Rcout << "Pass Nodes from " << x << " to " << y << " with z = [ ";
    // 	// for (const Node& node : z) {
    // 	//     Rcpp::Rcout << node.getName() << " ";
    // 	// }
    // 	// Rcpp::Rcout << "] :\n  ";
    // 	// for (const Node& node : passNodes) {
    // 	//     Rcpp::Rcout << node.getName() << " ";
    // 	// }

    // 	if (!b.isObserved() || std::find(z.begin(), z.end(), b) != z.end()) {
    // 	    for (const Node& c : passNodes) {
    // 		nodeStack.push(c);
    // 		pathStack.push(path);
    // 		sepsetStack.push(z);
    // 	    }
	    
    // 	} else {
    // 	    for (const Node& c : passNodes) {
    // 		nodeStack.push(c);
    // 		pathStack.push(path);
    // 		sepsetStack.push(z);
    // 	    }

    // 	    path.pop_back();
    // 	    z.push_back(b);

    // 	    passNodes = getPassNodes(a, b, z, graph);

    // 	    // Rcpp::Rcout << "Pass Nodes from " << x << " to " << y << " with z = [ ";
    // 	    // for (const Node& node : z) {
    // 	    // 	Rcpp::Rcout << node.getName() << " ";
    // 	    // }
    // 	    // Rcpp::Rcout << "] :\n  ";
    // 	    // for (const Node& node : passNodes) {
    // 	    // 	Rcpp::Rcout << node.getName() << " ";
    // 	    // }

    // 	    for (const Node& c : passNodes) {
    // 		nodeStack.push(c);
    // 		pathStack.push(path);
    // 		sepsetStack.push(z);
    // 	    }
    // 	}
    // }
	
    // 	if (path.back() == y)
    // 	    sepset = z;
    // }

    // Rcpp::Rcout << "Final z for " << x << " and " << y << ":\n  ";
    // for (const Node& node : z) {
    // 	Rcpp::Rcout << node.getName() << " ";
    // }
    // Rcpp::Rcout << std::endl;

    // if (path.count(y))
    // 	sepset = z;

    return z;
}

bool GraphUtils::sepsetPathFound(const Node& a, const Node& b, const Node& y,
				 std::set<Node>& path, std::vector<Node>& z,
				 std::set<Triple>& colliders, EdgeListGraph& graph) {

    if (b==y) return true;

    if (path.count(b)) return false;

    if (path.size() > 10) return false;

    path.insert(b);

    std::vector<Node> passNodes;

    if (!b.isObserved() || std::find(z.begin(), z.end(), b) != z.end()) {
	passNodes = getPassNodes(a, b, z, graph);
	for (const Node& c : passNodes) {
	    if (sepsetPathFound(b, c, y, path, z, colliders, graph)) {
		path.erase(b);
		return true;
	    }
	}

	path.erase(b);
	return false;
    } else {
	bool found1 = false;
	std::set<Triple> colliders1;

	passNodes = getPassNodes(a, b, z, graph);

	for (const Node& c : passNodes) {
	    if (sepsetPathFound(b, c, y, path, z, colliders1, graph)) {
	        found1 = true;
		break;
	    }
	}

	if (!found1) {
	    path.erase(b);
	    colliders.insert(colliders1.begin(), colliders1.end());
	    return false;
	}

	z.push_back(b);
	bool found2 = false;
	std::set<Triple> colliders2;

	passNodes = getPassNodes(a, b, z, graph);

	for (const Node& c : passNodes) {
	    if (sepsetPathFound(b, c, y, path, z, colliders2, graph)) {
	        found2 = true;
		break;
	    }
	}

	if (!found2) {
	    path.erase(b);
	    colliders.insert(colliders2.begin(), colliders2.end());
	    return false;
	}

	path.erase(b);
	z.pop_back();
	return true;
    }
}

std::vector<Node> GraphUtils::getPassNodes(const Node& a, const Node& b,
					   std::vector<Node>& z,
					   EdgeListGraph& graph) {
    std::vector<Node> passNodes;

    for (Node c : graph.getAdjacentNodes(b)) {
	if (a == c) continue;

	if (reachable(a, b, c, z, graph)) {
	    passNodes.push_back(c);
	}
    }

    return passNodes;
}

bool GraphUtils::reachable(const Node& a, const Node& b, const Node& c,
			   std::vector<Node>& z,
			   EdgeListGraph& graph) {
    bool collider = graph.isDefCollider(a,b,c);

    if (!collider && std::find(z.begin(), z.end(), b) == z.end()) {
	return true;
    }

    bool ancestor = isAncestor(b, z, graph);

    // Rcpp::Rcout << "Is " << b << " an ancestor of [ ";
    // for (const Node& node : z) {
    // 	Rcpp::Rcout << node.getName() << " ";
    // }
    // Rcpp::Rcout << "] : " << ancestor << std::endl;

    return collider && ancestor;
}


bool GraphUtils::isAncestor(const Node& b, std::vector<Node>& z, EdgeListGraph& graph) {
    if (std::find(z.begin(), z.end(), b) != z.end()) return true;

    std::queue<Node> nodeQueue;
    std::unordered_set<Node> visited;

    for (const Node& node : z) {
	nodeQueue.push(node);
	visited.insert(node);
    }

    while (!nodeQueue.empty()) {
	Node t = nodeQueue.front();
	nodeQueue.pop();

	if (t==b) return true;

	for (Node c : graph.getParents(t)) {
	    if (visited.count(c) == 0) {
		nodeQueue.push(c);
		visited.insert(c);
	    }
	}
    }

    return false;
}


std::vector<Node> GraphUtils::getInducingPath(const Node& x, const Node& y,
					      EdgeListGraph& graph) {
    if (!x.isObserved())
	throw std::runtime_error(x.getName() + " is an unobserved variable");
    if (!y.isObserved())
	throw std::runtime_error(y.getName() + " is an unobserved variable");

    std::vector<Node> inducingPath = {};
    if (x == y) return inducingPath;
    if (graph.getAdjacentNodes(x).size() == 0) return inducingPath;
    if (graph.getAdjacentNodes(y).size() == 0) return inducingPath;

    std::list<Node> path;
    std::set<std::list<Node>> visited;
    std::stack<Node> nodeStack;
    std::stack<std::list<Node>> pathStack;
    
    for (Node b : graph.getAdjacentNodes(x)) {
	path = {x};
	nodeStack.push(b);
	pathStack.push(path);
	visited.insert(path);
    }
    
    while (!nodeStack.empty()) {
	Node b = nodeStack.top();
	nodeStack.pop();
	path = pathStack.top();
	pathStack.pop();
	Node a = path.back();
	
	path.push_back(b);

	if (visited.count(path) > 0) continue;
	
	visited.insert(path);

	if (b == y) break;

	for (Node c : graph.getAdjacentNodes(b)) {
	    if (c == a) continue;

	    if (b.isObserved()) {
		if (!graph.isDefCollider(a, b, c)) {
		    continue;
		}
	    }

	    if (graph.isDefCollider(a, b, c)) {
		if (!(graph.isAncestorOf(b, x) || graph.isAncestorOf(b, y))) {
		    continue;
		}
	    }

	    if (std::find(path.begin(), path.end(), c) != path.end()) continue;

	    nodeStack.push(c);
	    pathStack.push(path);
	}

	// path.pop_back();
    }

    if (path.back() == y)
	inducingPath = std::vector<Node>(path.begin(), path.end());

    return inducingPath;
}
