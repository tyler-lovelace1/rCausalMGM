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
