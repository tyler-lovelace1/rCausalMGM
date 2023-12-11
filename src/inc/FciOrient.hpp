#ifndef FCIORIENT_HPP_
#define FCIORIENT_HPP_

/**
 * Extends Erin Korber's implementation of the Fast Causal Inference algorithm
 * (found in FCI.java) with Jiji Zhang's Augmented FCI rules (found in sec. 4.1
 * of Zhang's 2006 PhD dissertation, "Causal Inference and Reasoning in Causally
 * Insufficient Systems").
 * <p>
 * This class is based off a copy of FCI.java taken from the repository on
 * 2008/12/16, revision 7306. The extension is done by extending
 * doFinalOrientation() with methods for Zhang's rules R5-R10 which implements
 * the augmented search. (By a remark of Zhang's, the rule applications can be
 * staged in this way.)
 *
 * @author Erin Korber, June 2004
 * @author Alex Smith, December 2008
 * @author Joseph Ramsey
 * @author Choh-Man Teng
 * @author Jack Fiore, September 2020
 */

#include "EdgeListGraph.hpp"
#include "Knowledge.hpp"
#include "IndependenceTest.hpp"
#include "SepsetMap.hpp"
#include "ChoiceGenerator.hpp"
#include "SepsetProducer.hpp"

class FciOrient {

private:

    /**
     * The SepsetMap being constructed.
     * Needs to be made
     * Probably need to pass in a pointer
     * Use virtual for classes to be inherited but not implemented
     * use of protected: 'anything it needs to store'
     */
    SepsetProducer sepsets;

    Knowledge knowledge;

    bool changeFlag = true;

    /**
     * flag for complete rule set, true if should use complete rule set, false
     * otherwise.
     */
    bool completeRuleSetUsed = false;

    /**
     * True iff the possible dsep search is done.
     */
    bool possibleDsepSearchDone = true;

    /**
     * True iff the formation of cycles and almost cycles is strictly enforced.
     */
    bool aggressivelyPreventCycles = true;

    /**
     * The maximum length for any discriminating path. -1 if unlimited;
     * otherwise, a positive integer.
     */
    int maxPathLength = -1;

    /**
     *  For profiling, to find relative effectiveness of each rule
     */
    std::unordered_map<std::string,std::string> whyOrient;

    EdgeListGraph truePag;
    EdgeListGraph dag;
    bool skipDiscriminatingPathRule = false;

    //===========================PRIVATE METHODS=========================//
    std::vector<Node> getSepset(const Node& i, const Node& k) { return sepsets.getSepset(i, k); }

    void printWrongColliderMessage(const Node& a, const Node& b, const Node& c, std::string location, EdgeListGraph& graph);

    void spirtesFinalOrientation(EdgeListGraph& graph);

    void zhangFinalOrientation(EdgeListGraph& graph);

     /// R1, away from collider
     // If a*->bo-*c and a, c not adjacent then a*->b->c
    void ruleR1(const Node& a, const Node& b, const Node& c, EdgeListGraph& graph);

    bool isNoncollider(const Node& a,
		       const Node& b,
		       const Node& c) { return sepsets.isNoncollider(a, b, c); }

    bool isCollider(const Node& a,
		    const Node& b,
		    const Node& c) { return sepsets.isCollider(a, b, c); }

     //if a*-oc and either a-->b*->c or a*->b-->c, then a*->c
     // This is Zhang's rule R2.
    void ruleR2(const Node& a, const Node& b, const Node& c, EdgeListGraph& graph);

    /**
     * The triangles that must be oriented this way (won't be done by another
     * rule) all look like the ones below, where the dots are a collider path
     * from L to A with each node on the path (except L) a parent of C.
     * <pre>
     *          B
     *         xo           x is either an arrowhead or a circle
     *        /  \
     *       v    v
     * L....A --> C
     * </pre>
     * <p>
     * This is Zhang's rule R4, discriminating undirectedPaths.
     */
    void ruleR4A(EdgeListGraph& graph);

    /**
     * a method to search "back from a" to find a DDP. It is called with a
     * reachability list (first consisting only of a). This is breadth-first,
     * utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. </p>
     * The body of a DDP consists of colliders that are parents of c.
     */
    void reachablePathFind(const Node& a, const Node& b, const Node& c,
                           std::list<Node> reachable, EdgeListGraph& graph);

    /**
     * Orients the edges inside the definte discriminating path triangle. Takes
     * the left endpoint, and a,b,c as arguments.
     */
    void doDdpOrientation(const Node& d, const Node& a, const Node& b, const Node& c, EdgeListGraph& graph);

    /**
     * Orients the edges inside the definte discriminating path triangle. Takes
     * the left endpoint, and a,b,c as arguments.
     */
    bool doDdpOrientation(const Node& d, const Node& a, const Node& b, const Node& c,
                             std::unordered_map<Node, Node> previous, EdgeListGraph& graph);


    void printDdp(const Node& d, std::vector<Node> path,
		  const Node& a, const Node& b, const Node& c,
		  EdgeListGraph& graph);


    std::vector<Node> getPath(const Node& c, std::unordered_map<Node, Node> previous);

    /**
     * Orients every edge on a path as undirected (i.e. A---B).
     * <p>
     * DOES NOT CHECK IF SUCH EDGES ACTUALLY EXIST: MAY DO WEIRD THINGS IF
     * PASSED AN ARBITRARY LIST OF NODES THAT IS NOT A PATH.
     *
     * @param path The path to orient as all tails.
     */
    void orientTailPath(std::vector<Node> path, EdgeListGraph& graph);

    /**
     * Gets a list of every uncovered partially directed path between two nodes
     * in the graph.
     * <p>
     * Probably extremely slow.
     *
     * @param n1 The beginning node of the undirectedPaths.
     * @param n2 The ending node of the undirectedPaths.
     * @return A list of uncovered partially directed undirectedPaths from n1 to
     * n2.
     */
    std::vector<std::vector<Node>> getUcPdPaths(const Node& n1,
						const Node& n2,
						EdgeListGraph& graph);

    /**
     * Used in getUcPdPaths(n1,n2) to perform a breadth-first search on the
     * graph.
     * <p>
     * ASSUMES soFar CONTAINS AT LEAST ONE NODE!
     * <p>
     * Probably extremely slow.
     *
     * @param curr The getModel node to test for addition.
     * @param soFar The getModel partially built-up path.
     * @param end The node to finish the undirectedPaths at.
     * @param ucPdPaths The getModel list of uncovered p.d. undirectedPaths.
     */
    void getUcPdPsHelper(const Node& curr, std::vector<Node> soFar, const Node& end,
                         std::vector<std::vector<Node>> ucPdPaths, EdgeListGraph& graph);

    /**
     * Gets a list of every uncovered circle path between two nodes in the graph
     * by iterating through the uncovered partially directed undirectedPaths and
     * only keeping the circle undirectedPaths.
     * <p>
     * Probably extremely slow.
     *
     * @param n1 The beginning node of the undirectedPaths.
     * @param n2 The ending node of the undirectedPaths.
     * @return A list of uncovered circle undirectedPaths between n1 and n2.
     */
    std::vector<std::vector<Node>> getUcCirclePaths(const Node& n1,
						    const Node& n2,
						    EdgeListGraph& graph);

    /**
     * Tries to apply Zhang's rule R8 to a pair of nodes A and C which are
     * assumed to be such that Ao->C.
     * <p>
     * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
     * <p>
     * R8: If Ao->C and A-->B-->C or A--oB-->C, then A-->C.
     *
     * @param a The node A.
     * @param c The node C.
     * @return Whether or not R8 was successfully applied.
     */
    bool ruleR8(const Node& a, const Node& c, EdgeListGraph& graph);

    /**
     * Tries to apply Zhang's rule R9 to a pair of nodes A and C which are
     * assumed to be such that Ao->C.
     * <p>
     * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
     * <p>
     * R9: If Ao->C and there is an uncovered p.d. path u=<A,B,..,C> such that
     * C,B nonadjacent, then A-->C.
     *
     * @param a The node A.
     * @param c The node C.
     * @return Whether or not R9 was succesfully applied.
     */
    bool ruleR9(const Node& a, const Node& c, EdgeListGraph& graph);

    /**
     * Tries to apply Zhang's rule R10 to a pair of nodes A and C which are
     * assumed to be such that Ao->C.
     * <p>
     * MAY HAVE WEIRD EFFECTS ON ARBITRARY NODE PAIRS.
     * <p>
     * R10: If Ao->C, B-->C<--D, there is an uncovered p.d. path u1=<A,M,...,B>
     * and an uncovered p.d. path u2=
     * <A,N,...,D> with M != N and M,N nonadjacent then A-->C.
     *
     * @param a The node A.
     * @param c The node C.
     * @return Whether or not R10 was successfully applied.
     */
    bool ruleR10(const Node& a, const Node& c, EdgeListGraph& graph);

    /**
     * Helper method. Appears to check if an arrowpoint is permitted by
     * background knowledge.
     *
     * @param x The possible other node.
     * @param y The possible point node.
     * @return Whether the arrowpoint is allowed.
     */
    bool isArrowpointAllowed(const Node& x, const Node& y, EdgeListGraph& graph);

    bool doesNotCreateCycle(Node from, Node to, EdgeListGraph& graph);

    bool doesNotCreateAlmostCycle(Node from, Node to, EdgeListGraph& graph);

    //============================CONSTRUCTORS============================//

public:

    FciOrient() {}
    /**
     * Constructs a new FCI search for the given independence test and
     * background knowledge.
     */
    FciOrient(SepsetProducer& sepsets);

    FciOrient(SepsetProducer& sepsets, std::unordered_map<std::string,std::string> whyOrient);


    //========================PUBLIC METHODS==========================//
    EdgeListGraph orient(EdgeListGraph& graph);

    SepsetProducer getSepsets() { return sepsets; }

    /**
     * @return true if Zhang's complete rule set should be used, false if only
     * R1-R4 (the rule set of the original FCI) should be used. False by
     * default.
     */
    bool isCompleteRuleSetUsed() { return completeRuleSetUsed; }

    /**
     * @param completeRuleSetUsed set to true if Zhang's complete rule set
     * should be used, false if only R1-R4 (the rule set of the original FCI)
     * should be used. False by default.
     */
    void setCompleteRuleSetUsed(bool completeRuleSetUsed) { this->completeRuleSetUsed = completeRuleSetUsed; }


    /**
     * Orients colliders in the graph. (FCI Step C)
     * <p>
     * Zhang's step F3, rule R0.
     */
    void ruleR0(EdgeListGraph& graph);

    /**
     * Orients the graph according to rules in the graph (FCI step D).
     * <p>
     * Zhang's step F4, rules R1-R10.
     */
    void doFinalOrientation(EdgeListGraph& graph);


    //Does all 3 of these rules at once instead of going through all
    // triples multiple times per iteration of doFinalOrientation.
    void rulesR1R2cycle(EdgeListGraph& graph);


    /**
     * Implements the double-triangle orientation rule, which states that if
     * D*-oB, A*->B<-*C and A*-oDo-*C, then
     * D*->B.
     * <p>
     * This is Zhang's rule R3.
     */
    void ruleR3(EdgeListGraph& graph);


    /**
     * The triangles that must be oriented this way (won't be done by another
     * rule) all look like the ones below, where the dots are a collider path
     * from L to A with each node on the path (except L) a parent of C.
     * <pre>
     *          B
     *         xo           x is either an arrowhead or a circle
     *        /  \
     *       v    v
     * L....A --> C
     * </pre>
     * <p>
     * This is Zhang's rule R4, discriminating undirectedPaths.
     */
    void ruleR4B(EdgeListGraph& graph);

    /**
     * a method to search "back from a" to find a DDP. It is called with a
     * reachability list (first consisting only of a). This is breadth-first,
     * utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. </p>
     * The body of a DDP consists of colliders that are parents of c.
     */
    void ddpOrient(const Node& a, const Node& b, const Node& c, EdgeListGraph& graph);


    /**
     * Implements Zhang's rule R5, orient circle undirectedPaths: for any Ao-oB,
     * if there is an uncovered circle path u =
     * <A,C,...,D,B> such that A,D nonadjacent and B,C nonadjacent, then A---B
     * and orient every edge on u undirected.
     */
    void ruleR5(EdgeListGraph& graph);

    /**
     * Implements Zhang's rules R6 and R7, applies them over the graph once.
     * Orient single tails. R6: If A---Bo-*C then A---B--*C. R7: If A--oBo-*C
     * and A,C nonadjacent, then A--oB--*C
     */
    void ruleR6R7(EdgeListGraph& graph);

    /**
     * Implements Zhang's rules R8, R9, R10, applies them over the graph once.
     * Orient arrow tails. I.e., tries R8, R9, and R10 in that sequence on each
     * Ao->C in the graph.
     */
    void rulesR8R9R10(EdgeListGraph& graph);

    void fciOrientbk(Knowledge& knowledge, EdgeListGraph& graph);

    void setKnowledge(Knowledge& knowledge) { this->knowledge = knowledge; }

    bool isPossibleDsepSearchDone() { return possibleDsepSearchDone; }

    void setPossibleDsepSearchDone(bool possibleDsepSearchDone) { this->possibleDsepSearchDone = possibleDsepSearchDone; }

    /**
     * @return the maximum length of any discriminating path, or -1 of
     * unlimited.
     */
    int getMaxPathLength() { return maxPathLength; }

    /**
     * @param maxPathLength the maximum length of any discriminating path, or -1
     * if unlimited.
     */
    void setMaxPathLength(int maxPathLength);


    void setTruePag(EdgeListGraph& truePag) { this->truePag = truePag; }

    /**
     * The true PAG if available. Can be null.
     */
    EdgeListGraph getTruePag() { return truePag; }

    void setChangeFlag(bool changeFlag) { this->changeFlag = changeFlag; }

    /**
     * change flag for repeat rules
     */
    bool isChangeFlag() { return changeFlag; }

    void setskipDiscriminatingPathRule(bool skip) { this->skipDiscriminatingPathRule = skip; }

};

#endif /* FCIORIENT_HPP_ */
