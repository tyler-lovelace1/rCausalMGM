#ifndef MEEKRULES_HPP_
#define MEEKRULES_HPP_

#include "EdgeListGraph.hpp"
#include "ChoiceGenerator.hpp"
#include "Knowledge.hpp"
// #include "SepsetProducer.hpp"
#include "SearchGraphUtils.hpp"
#include <stack>

/**
 * Implements Meek's complete orientation rule set for PC (Chris Meek (1995), "Causal inference and causal explanation
 * with background knowledge"), modified for Conservative PC to check noncolliders against recorded noncolliders before
 * orienting.
 * <p>
 * For now, the fourth rule is always performed.
 *
 * @author Joseph Ramsey
 * @author converted to C++ by Max Dudek
 */
class MeekRules {

private:
    // True if cycles are to be aggressively prevented. May be
    // expensive for large graphs (but also useful for large graphs).
    bool aggressivelyPreventCycles = true;

    // If knowledge is available.
    bool useRule4 = false;

    /**
     * The SepsetMap being constructed.
     * Needs to be made
     * Probably need to pass in a pointer
     * Use virtual for classes to be inherited but not implemented
     * use of protected: 'anything it needs to store'
     */
    // SepsetProducer sepsets;

    Knowledge knowledge;

    // The stack of nodes to be visited.
    std::stack<Node> directStack;

    // The initial list of nodes to visit.
    std::vector<Node> nodes;

    // The list of nodes actually visited.
    std::unordered_set<Node> visited;

    // Edges already oriented by the algorithm to avoid repeats and prevent cycles.
    std::unordered_set<Edge> oriented;

    // True if unforced parents should be undirected before orienting.
    bool shouldUndirectUnforcedEdges = false;

    void orientUsingMeekRulesLocally(EdgeListGraph& graph);

    void runMeekRules(Node node, EdgeListGraph& graph);

    /**
     * Meek's rule R1: if a-->b, b---c, and a not adj to c, then a-->c
     */
    void meekR1(Node b, EdgeListGraph& graph);
    void r1Helper(Node a, Node b, Node c, EdgeListGraph& graph);

    /**
     * If a-->b-->c, a--c, then b-->c.
     */
    void meekR2(Node c, EdgeListGraph& graph);
    void r2Helper(Node a, Node b, Node c, EdgeListGraph& graph);

    /**
     * Meek's rule R3. If a--b, a--c, a--d, c-->b, d-->b, then orient a-->b.
     */
    void meekR3(Node a, EdgeListGraph& graph);
    bool isKite(Node a, Node d, Node b, Node c, EdgeListGraph& graph);

    void meekR4(Node a, EdgeListGraph& graph);

    void direct(Node a, Node c, EdgeListGraph& graph);

    bool isUnshieldedNoncollider(Node a, Node b, Node c, EdgeListGraph& graph);

    bool isArrowpointAllowed(Node from, Node to);

    bool doesNotCreateCycle(Node from, Node to, EdgeListGraph& graph);

    void undirectUnforcedEdges(Node y, EdgeListGraph& graph);

public:

    MeekRules() {}

    // MeekRules(SepsetProducer& sp) : sepsets(sp) {}

    // void meekR0(EdgeListGraph& graph);

    void orientImplied(EdgeListGraph& graph);

    void orientImplied(EdgeListGraph& graph, std::vector<Node>& nodes);

    bool isAggressivelyPreventCycles() { return aggressivelyPreventCycles; }
    void setAggressivelyPreventCycles(bool aggressivelyPreventCycles) { this->aggressivelyPreventCycles = aggressivelyPreventCycles; }

    bool isUndirectUnforcedEdges() { return shouldUndirectUnforcedEdges; }
    void setUndirectUnforcedEdges(bool shouldUndirectUnforcedEdges) { this->shouldUndirectUnforcedEdges = shouldUndirectUnforcedEdges; }

    std::unordered_set<Node> getVisited() { return visited; }

    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }


};

#endif /* MEEKRULES_HPP_ */
