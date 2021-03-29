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
std::vector<Variable*> GraphUtils::asList(std::vector<int>& indices, std::vector<Variable*>& nodes) {
    std::vector<Variable*> list;

    for (int i : indices) {
        list.push_back(nodes[i]);
    }

    return list;
}

std::unordered_set<Variable*> GraphUtils::asSet(std::vector<int>& indices, std::vector<Variable*>& nodes) {
    std::unordered_set<Variable*> set;

    for (int i : indices) {
        set.insert(nodes[i]);
    }

    return set;
}

EdgeListGraph GraphUtils::completeGraph(EdgeListGraph& graph) {
    EdgeListGraph graph2(graph.getNodes());

    graph2.removeEdges();

    std::vector<Variable*> nodes = graph2.getNodes();

    for (int i = 0; i < nodes.size(); i++) {
        for (int j = i+1; j < nodes.size(); j++) {
            Variable* node1 = nodes[i];
            Variable* node2 = nodes[j];
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

std::unordered_set<Variable*> GraphUtils::possibleDsep(Variable* x, Variable* y, EdgeListGraph& graph, int maxPathLength) {
    std::unordered_set<Variable*> dsep;

    std::queue<boost::optional<std::pair<Variable*,Variable*>>> Q;
    std::unordered_set<std::pair<Variable*,Variable*>, boost::hash<std::pair<Variable*, Variable*> > > V;

    std::unordered_map<Variable*, std::vector<Variable*>> previous;
    std::vector<Variable*> null_vector = {};
    previous.insert(std::pair(x, null_vector));

    boost::optional<std::pair<Variable*,Variable*>> e = {};
    int distance = 0;

    for (Variable* b : graph.getAdjacentNodes(x)) {
        if (b == y) {
            continue;
        }
        boost::optional<std::pair<Variable*,Variable*> > edge = std::pair<Variable*,Variable*>(x, b);
        if (!e) {
            e = edge;
        }
        Q.push(edge);
        V.insert(*edge);
        addToList(previous, b, x);
        dsep.insert(b);
    }
    while (!Q.empty()) {
        boost::optional<std::pair<Variable*, Variable*> > t = Q.front();
        Q.pop();

        if (e == t) {
            e = {};
            distance++;
            if (distance > 0 && distance > (maxPathLength == -1 ? 1000 : maxPathLength)) {
                break;
            }
        }

        Variable* a = t->first;
        Variable* b = t->second;
        if (existOnePathWithPossibleParents(previous, b, x, b, graph)) {
            dsep.insert(b);
        }

        for (Variable* c : graph.getAdjacentNodes(b)) {
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
                boost::optional<std::pair<Variable*, Variable*> > u = std::pair<Variable*, Variable*>(a,c);
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

void GraphUtils::addToList(std::unordered_map<Variable*, std::vector<Variable*>> previous, Variable* b, Variable* c) {
    std::vector<Variable*> list = previous[c];

    if (list.empty()) {
        std::vector<Variable*> null_list = {};
        list = null_list;
    }

    list.push_back(b);
}

bool GraphUtils::existOnePathWithPossibleParents(std::unordered_map<Variable*, std::vector<Variable*>> previous, Variable* w, Variable* x, Variable* b, EdgeListGraph& graph) {
    if (w == x) {
        return true;
    }
    const std::vector<Variable*> p = previous[w];
    if (p.empty()) {
        return false;
    }

    for (Variable* r : p) {
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

bool GraphUtils::existsSemidirectedPath(Variable* from, Variable* to, EdgeListGraph& G) {
    std::queue<Variable*> Q;

    std::unordered_set<Variable*> V;
    Q.push(from);
    V.insert(from);

    while (!Q.empty()) {
        Variable* t = Q.front();
        Q.pop();
        if (t == to) {
            return true;
        }
        for (Variable* u : G.getAdjacentNodes(t)) {
            Edge edge = G.getEdge(t, u);
            Variable* c = Edge::traverseSemiDirected(t, edge);

            if (c == NULL) {
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
