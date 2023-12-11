#include "Knowledge.hpp"

Knowledge::Knowledge(std::vector<Node> nodes, Rcpp::List& knowledge) : nodes(nodes) {
    
    if (knowledge.size()!=0) {
	std::vector<std::string> names = knowledge.names();

	if (std::find(names.begin(), names.end(), "tiers") != names.end() &&
	    !Rf_isNull(knowledge["tiers"])) {
	    std::set<Node> usedNodes;
	    std::set<Node> nodeSet(nodes.begin(), nodes.end());
	    Rcpp::List tierList = knowledge["tiers"];

	    for (int i = 0; i < tierList.size(); i++) {
		tiers.push_back(std::set<Node>());
		std::vector<std::string> nodeNames = tierList[i];
		for (std::string nodeName : nodeNames) {
		    Node n(new ContinuousVariable(nodeName));
		    if (usedNodes.count(n)) {
			throw std::invalid_argument("Invalid Knowledge: Node " + nodeName + " is in more than one tier.");
		    } else {
			auto it = nodeSet.find(n);
			if (it == nodeSet.end()) {
			    throw std::invalid_argument("Invalid Knowledge: Node " + nodeName + " is not in the dataset.");
			}
			tiers[i].insert(*it);
			usedNodes.insert(*it);
		    }
		}
	    }

	    if (usedNodes != nodeSet) {
		throw std::invalid_argument("Invalid Knowledge: All variables in the dataset must be in the tier list if it is not NULL.");
	    }
	
	    if (std::find(names.begin(), names.end(), "forbiddenWithinTier") != names.end() &&
		!Rf_isNull(knowledge["forbiddenWithinTier"])) {
		Rcpp::LogicalVector _forbiddenWithinTier = knowledge["forbiddenWithinTier"];
		forbiddenWithinTier = Rcpp::as<std::vector<bool>>(_forbiddenWithinTier);
		if (forbiddenWithinTier.size() != tiers.size()) {
		    throw std::invalid_argument("Invalid Knowledge: The list of tiers and the logical vector forbiddenWithinTier must be the same length.");
		}
	    } else {
		for (int i = 0; i < tiers.size(); i++) {
		    forbiddenWithinTier.push_back(false);
		}
	    }
	}
    
	if (std::find(names.begin(), names.end(), "forbidden") != names.end() &&
	    !Rf_isNull(knowledge["forbidden"])) {
	    Rcpp::List forbiddenList = knowledge["forbidden"];
	    for (int i = 0; i < forbiddenList.size(); i++) {
		std::vector<std::string> nodeNames = forbiddenList[i];
		if (nodeNames.size() != 2) {
		    throw std::invalid_argument("Invalid Knowledge: Forbidden edges must contain the names of two nodes.");
		}
		forbidden.insert(
		    std::pair<Node,Node>(Node(new ContinuousVariable(nodeNames[0])),
					 Node(new ContinuousVariable(nodeNames[1])))
		    );
	    }
	}
    
	if (std::find(names.begin(), names.end(), "required") != names.end() &&
	    !Rf_isNull(knowledge["required"])) {
	    Rcpp::List requiredList = knowledge["required"];
	    for (int i = 0; i < requiredList.size(); i++) {
		std::vector<std::string> nodeNames = requiredList[i];
		if (nodeNames.size() != 2) {
		    throw std::invalid_argument("Invalid Knowledge: Required edges must contain the names of two nodes.");
		}
		required.insert(
		    std::pair<Node,Node>(Node(new ContinuousVariable(nodeNames[0])),
					 Node(new ContinuousVariable(nodeNames[1])))
		    );
	    }
	}

	for (int i = 0; i < nodes.size(); i++) {
	    for (int j = 0; j < i; j++) {
		if (i==j) continue;
		if (isForbidden(nodes[i], nodes[j]) && isRequired(nodes[i], nodes[j])) {
		    throw std::invalid_argument("Invalid Knowledge: The edge " + nodes[i].getName() + " *-> " + nodes[j].getName() + " is both forbidden and required.");
		}
		if (isForbidden(nodes[j], nodes[i]) && isRequired(nodes[j], nodes[i])) {
		    throw std::invalid_argument("Invalid Knowledge: The edge " + nodes[j].getName() + " *-> " + nodes[i].getName() + " is both forbidden and required.");
		}
	    }
	}
    }
}

bool Knowledge::isForbidden(const Node& node1, const Node& node2) {
    return isForbiddenByTiers(node1, node2) ||
	forbidden.count(std::pair<Node,Node>(node1, node2));
}

bool Knowledge::isForbiddenByTiers(const Node& node1, const Node& node2) {
    int tier1=-1, tier2=-1;

    if (tiers.empty()) return false;

    for (int i = 0; i < tiers.size(); i++) {
	if (tiers[i].count(node1)) tier1 = i;
	if (tiers[i].count(node2)) tier2 = i;
	if (tier1 >= 0 && tier2 >= 0) break;
    }

    if (tier1==tier2) return forbiddenWithinTier[tier1];

    return tier1 > tier2;
}

bool Knowledge::isRequired(const Node& node1, const Node& node2) {
    return required.count(std::pair<Node,Node>(node1, node2));
}

bool Knowledge::noEdgeRequired(const Node& node1, const Node& node2) {
    return !(isRequired(node1, node2) || isRequired(node2, node1));
}
