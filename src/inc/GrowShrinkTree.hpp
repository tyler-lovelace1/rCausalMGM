#ifndef GROWSHRINKTREE_HPP_
#define GROWSHRINKTREE_HPP_

#include <memory>
#include "Score.hpp"

struct GrowShrinkNode {

    struct cmp_ptr {
        bool operator()(const GrowShrinkNode* lhs,
			const GrowShrinkNode* rhs) const {
            return *lhs < *rhs;
        }
    };
    
    Node node;
    
    double score = 0.0;
    std::atomic<double> shrinkScore;
    std::atomic_bool shrinkFlag;
    std::vector<std::atomic_bool> visited;
    std::vector<Node> condSet;
    std::vector<Node> shrinkCondSet;

    std::set<GrowShrinkNode*,cmp_ptr> children;
    
    GrowShrinkNode() {
	shrinkFlag.store(true);
	shrinkScore.store(0.0);
	// score.store(0.0);
    }

    GrowShrinkNode(int numNodes) : visited(std::vector<std::atomic_bool>(numNodes)) {
	shrinkFlag.store(true);
	shrinkScore.store(0.0);
	// score.store(0.0);
	for (int i = 0; i < numNodes; i++) {
	    visited[i].store(false);
	    // visited[i] = false;
	}
    }

    GrowShrinkNode(Node node,
		   std::vector<Node> condSet,
		   int numNodes) : node(node),
				   condSet(condSet),
				   visited(std::vector<std::atomic_bool>(numNodes)) {
	// score.store(0.0);
	for (int i = 0; i < numNodes; i++) {
	    visited[i].store(false);
	    // visited[i] = false;
	}
	this->condSet.push_back(node);
	shrinkFlag.store(false);
	// shrinkFlag = false;
    }

    // GrowShrinkNode(GrowShrinkNode& rhs) : node(rhs.node),
    // 					  score(rhs.score),
    // 					  visited(std::vector<std::atomic_bool>(rhs.visited.size())),
    // 					  condSet(rhs.condSet),
    // 					  shrinkCondSet(rhs.shrinkCondSet) {

    // 	shrinkScore.store(rhs.shrinkScore.load());
    // 	shrinkFlag.store(rhs.shrinkFlag.load());

    // 	for (int i = 0; i < visited.size(); i++) {
    // 	    visited[i].store(visited[i].load());
    // 	}

    // 	for (auto it = rhs.children.begin(); it != rhs.children.end(); it++) {
    // 	    children.insert(std::make_shared<GrowShrinkNode>(*(*it)));
    // 	}
    // }

    // GrowShrinkNode(GrowShrinkNode&& rhs) : node(std::move(rhs.node)),
    // 					   score(std::move(rhs.score)),
    // 					   shrinkScore(std::move(rhs.shrinkScore)),
    // 					   shrinkFlag(std::move(rhs.shrinkFlag)),
    // 					   visited(std::move(rhs.visited)),
    // 					   condSet(std::move(rhs.condSet)),
    // 					   shrinkCondSet(std::move(rhs.shrinkCondSet)) {

    // 	for (auto it = rhs.children.begin(); it != rhs.children.end(); it++) {
    // 	    children.insert(std::exchange((*it), nullptr));
    // 	}
    // }

    ~GrowShrinkNode() {
	// Rcpp::Rcout << "GrowShrinkNode " << node << " : " << score << " deleted.\n";
	for (auto it = children.begin(); it != children.end(); it++) {
	    delete (*it);
	}
	children.clear();
    }

    // ~GrowShrinkNode() {
    // 	for (auto it = children.begin(); it != children.end(); it++) {
    // 	    delete *it;
    // 	}
    // }

    // std::set<std::shared_ptr<GrowShrinkNode>,cmp_ptr> getChildren() { return children; };

    void addChild(GrowShrinkNode* newChild) {
	auto insertFlag = children.insert(newChild);
	if (!insertFlag.second) {
	    delete newChild;
	}
    }

    // void addChildren(std::set<std::shared_ptr<GrowShrinkNode>,cmp_ptr> newChildren) {
    // 	for (auto it = newChildren.begin(); it != newChildren.end(); it++) {
    // 	    addChild(*it);
    // 	    // auto insertFlag = children.insert(*it);
    // 	    // if (!insertFlag.second) {
    // 	    // 	delete *it;
    // 	    // }
    // 	}
    // }

    void shrink(Score* scorer, Node target, bool verbose) {
	if (!shrinkFlag.load()) {    
	    double oldScore = 1e20, curScore = score, tempScore;
	    bool changeFlag = false;

	    if (verbose) {
		Rcpp::Rcout << "    Score = " << score << " : { ";
		for (Node n : condSet) {
		    Rcpp::Rcout << n << " ";
		}
		Rcpp::Rcout << "}\n";
	    }
    
	    std::list<Node> activeList(condSet.begin(), condSet.end());

	    std::list<Node>::iterator removeIt;

	    while (curScore < oldScore) {
		changeFlag = false;
		oldScore = curScore;
		for (auto it = activeList.begin(); it != activeList.end(); it++) {
		    std::vector<Node> tempActive;
		    for (auto jt = activeList.begin(); jt != activeList.end(); jt++) {
			if (jt != it) {
			    tempActive.push_back(*jt);
			}
		    }
	    
		    tempScore = scorer->localScore(target, tempActive);
	    
		    if (tempScore < curScore) {
			curScore = tempScore;
			removeIt = it;
			changeFlag = true;
		    }	
		}
		if (changeFlag) {
		    activeList.erase(removeIt);
		    if (verbose) {
			Rcpp::Rcout << "    Score = " << curScore << " : { ";
			for (Node n : activeList) {
			    Rcpp::Rcout << n << " ";
			}
			Rcpp::Rcout << "}\n";
		    }
		}
	    }

	    shrinkCondSet = std::vector<Node>(activeList.begin(), activeList.end());
	    shrinkScore.store(curScore);
	    shrinkFlag.store(true);
	    // shrinkScore = curScore;
	    // shrinkFlag = shrinkFlag
	}
    }

    // void clear() {
    // 	for (auto it = children.begin(); it != children.end(); it++) {
    // 	    (*it)->clear();
    // 	}
    // 	children.clear();
    // 	// delete this;
    // }

    bool operator==(const GrowShrinkNode& rhs) const {
	return score == rhs.score && node == rhs.node;
    }
    
    bool operator<(const GrowShrinkNode& rhs) const {
	if (condSet.size() == rhs.condSet.size()) {
	    if (score == rhs.score) {
		return node < rhs.node;
	    }
	    return score < rhs.score;
	}
	return condSet.size() < rhs.condSet.size();
    }
};

class GrowShrinkTree {
private:
    Score* scorer;
    bool verbose = false;
    int numNodes;
    Node target;
    std::unique_ptr<GrowShrinkNode> root;
    std::unordered_map<Node,int> node2idx;
    std::mutex treeMutex;
    
    // std::map<GrowShrinkNode, std::set<Node>> childSet;
    // std::map<GrowShrinkNode, std::set<GrowShrinkNode>> children;
    // std::mutex treeMutex;
    
public:

    GrowShrinkTree() {}

    GrowShrinkTree(Score* scorer, Node target);

    // GrowShrinkTree(GrowShrinkTree& rhs) {
    // 	scorer = rhs.scorer;
    // 	verbose = rhs.verbose;
    // 	numNodes = rhs.numNodes;
    // 	target = rhs.target;
    // 	node2idx = rhs.node2idx;

    // 	root = new GrowShrinkNode(*rhs.root);
    // }

    // GrowShrinkTree(GrowShrinkTree&& rhs) : root(std::exchange(rhs.root, nullptr)) {
    // 	scorer = rhs.scorer;
    // 	verbose = std::move(rhs.verbose);
    // 	numNodes = std::move(rhs.numNodes);
    // 	target = std::move(rhs.target);
    // 	node2idx = std::move(rhs.node2idx);

    // 	// root = new GrowShrinkNode(std::move(*rhs.root));
    // }

    // GrowShrinkTree& operator=(GrowShrinkTree& rhs) {
    // 	scorer = rhs.scorer;
    // 	verbose = rhs.verbose;
    // 	numNodes = rhs.numNodes;
    // 	target = rhs.target;
    // 	node2idx = rhs.node2idx;

    // 	delete root;
    // 	root = new GrowShrinkNode(*rhs.root);

    //     return *this;
    // }

    // GrowShrinkTree& operator=(GrowShrinkTree&& rhs) {
    //     if (this != &rhs) {
    // 	    scorer = rhs.scorer;
    // 	    verbose = std::move(rhs.verbose);
    // 	    numNodes = std::move(rhs.numNodes);
    // 	    target = std::move(rhs.target);
    // 	    node2idx = std::move(rhs.node2idx);

    //         delete root;
    //         root = rhs.root;
    //         rhs.root = nullptr;
    //     }
    //     return *this;
    // }

    // ~GrowShrinkTree() {
    // 	Rcpp::Rcout << "GrowShrinkTree " << target << " deleted.\n";
    // 	root.reset();
    // 	// delete root;
    // }

    void setVerbose(bool verbose) { this->verbose = verbose; }

    void setPenalty(bool penalty) { this->scorer->setPenalty(penalty); }

    // std::unique_ptr<GrowShrinkNode> grow(std::vector<Node>& candidates);

    std::vector<Node> growShrink(std::vector<Node>& candidates, double* scoreReturn);

    // std::vector<Node> shrink(std::vector<Node> active, double score, double* scoreReturn = NULL);

    std::vector<Node> search(std::vector<Node> candidates, double* scoreReturn = NULL);

    std::vector<Node> search(double* scoreReturn = NULL);

    void reset() {
	// childSet.clear();
	// children.clear();
	// delete root;
	// Rcpp::Rcout << "GrowShrinkTree reset called\n";
	// if (root.unique())  Rcpp::Rcout << "GrowShrinkTree root is unique\n";
	// else                Rcpp::Rcout << "GrowShrinkTree root has " << root.use_count() << " owners\n";
	// root->clear();
	root.reset(new GrowShrinkNode(numNodes));
    }

    friend Rcpp::StringVector GrowShrinkTreeTest(const Rcpp::DataFrame &df, std::string target);
    friend Rcpp::StringVector GrowShrinkTreeSubSetTest(const Rcpp::DataFrame &df, std::string target, int numSub);

    friend std::vector<double> GrowShrinkTreeParallelSubSetTest(const Rcpp::DataFrame &df, std::string target, int numSub, int threads);

};

#endif /* GROWSHRINKTREE_HPP_ */
