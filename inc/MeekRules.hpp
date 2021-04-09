#ifndef MEEKRULES_HPP_
#define MEEKRULES_HPP_

#include "EdgeListGraph.hpp"
#include "ChoiceGenerator.hpp"
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
    //True if cycles are to be aggressively prevented. May be expensive for large graphs (but also useful for large
    //graphs).
    bool aggressivelyPreventCycles = false;

    // If knowledge is available.
    bool useRule4 = false;

    // The stack of nodes to be visited.
    std::stack<Variable*> directStack;

    // The initial list of nodes to visit.
    std::vector<Variable*> nodes;

    // The lsit of nodes actually visited.
    std::unordered_set<Variable*> visited;

    // Edges already oriented by the algorithm to avoid repeats and prevent cycles.
    std::unordered_set<Edge> oriented;

    // True if unforced parents should be undirected before orienting.
    bool shouldUndirectUnforcedEdges = false;

    void orientUsingMeekRulesLocally(EdgeListGraph& graph);

    void runMeekRules(Variable* node, EdgeListGraph& graph);

    /**
     * Meek's rule R1: if a-->b, b---c, and a not adj to c, then a-->c
     */
    void meekR1(Variable* b, EdgeListGraph& graph);
    void r1Helper(Variable* a, Variable* b, Variable* c, EdgeListGraph& graph);

    /**
     * If a-->b-->c, a--c, then b-->c.
     */
    void meekR2(Variable* c, EdgeListGraph& graph);
    void r2Helper(Variable* a, Variable* b, Variable* c, EdgeListGraph& graph);

    /**
     * Meek's rule R3. If a--b, a--c, a--d, c-->b, d-->b, then orient a-->b.
     */
    void meekR3(Variable* a, EdgeListGraph& graph);
    bool isKite(Variable* a, Variable* d, Variable* b, Variable* c, EdgeListGraph& graph);

    void meekR4(Variable* a, EdgeListGraph& graph);

    void direct(Variable* a, Variable* c, EdgeListGraph& graph);

    static bool isUnshieldedNoncollider(Variable* a, Variable* b, Variable* c, EdgeListGraph& graph);

    static bool isArrowpointAllowed(Variable* from, Variable* to);

    void undirectUnforcedEdges(Variable* y, EdgeListGraph& graph);

public:

    MeekRules() {}

    void orientImplied(EdgeListGraph& graph);

    void orientImplied(EdgeListGraph& graph, std::vector<Variable*>& nodes);

    bool isAggressivelyPreventCycles() { return aggressivelyPreventCycles; }
    void setAggressivelyPreventCycles(bool aggressivelyPreventCycles) { this->aggressivelyPreventCycles = aggressivelyPreventCycles; }

    bool isUndirectUnforcedEdges() { return shouldUndirectUnforcedEdges; }
    void setUndirectUnforcedEdges(bool shouldUndirectUnforcedEdges) { this->shouldUndirectUnforcedEdges = shouldUndirectUnforcedEdges; }

    std::unordered_set<Variable*> getVisited() { return visited; }

};

#endif /* MEEKRULES_HPP_ */