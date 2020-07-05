#include "GraphUtils.hpp"

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