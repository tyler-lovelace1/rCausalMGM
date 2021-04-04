#include "FciOrient.hpp"
#include <algorithm>
#include <queue>


//============================CONSTRUCTORS============================//

/**
 * Constructs a new FCI search for the given independence test and
 * background knowledge.
 */
FciOrient::FciOrient(SepsetProducer* sepsets) {
  this->sepsets = sepsets;
  this->whyOrient = std::unordered_map<std::string,std::string>();

  //Different sepsetProducers
  // if (sepsets instanceof SepsetsGreedy) {
  //     SepsetsGreedy _sepsets = (SepsetsGreedy) sepsets;
  //     this->dag = _sepsets.getDag();
  // } else if (sepsets instanceof DagSepsets) {
  //     DagSepsets _sepsets = (DagSepsets) sepsets;
  //     this->edag = _sepsets.getDag();
  // }
}

FciOrient::FciOrient(SepsetProducer* sepsets, std::unordered_map<std::string,std::string> whyOrient) {
  this->sepsets = sepsets;
  this->whyOrient = whyOrient;
//   if (sepsets instanceof SepsetsGreedy) {
//       SepsetsGreedy _sepsets = (SepsetsGreedy) sepsets;
//       this.dag = _sepsets.getDag();
//   } else if (sepsets instanceof  DagSepsets) {
//       DagSepsets _sepsets = (DagSepsets) sepsets;
//       this.dag = _sepsets.getDag();
//   }
}

//============================PUBLIC METHODS============================//

EdgeListGraph FciOrient::orient(EdgeListGraph& graph) {
  ruleR0(graph);

  // Step CI D. (Zhang's step F4.)
  doFinalOrientation(graph);

  // graph.closeInducingPaths();   //to make sure it's a legal PAG

  return graph;
}

/**
 * Orients colliders in the graph. (FCI Step C)
 * <p>
 * Zhang's step F3, rule R0.
 */
void FciOrient::ruleR0(EdgeListGraph& graph) {
  graph.reorientAllWith(ENDPOINT_CIRCLE);
  // fciOrientbk(knowledge, graph, graph.getNodes());

  std::vector<Variable*> nodes = graph.getNodes();

  for (Variable* b : nodes) {
      std::vector<Variable*> adjacentNodes = graph.getAdjacentNodes(b);

      if (adjacentNodes.size() < 2) {
          continue;
      }

      ChoiceGenerator cg(adjacentNodes.size(), 2);
      std::vector<int>  *combination;

      for (combination = cg.next(); combination != NULL; combination = cg.next()) {
          Variable* a = adjacentNodes.at(combination->at(0));
          Variable* c = adjacentNodes.at(combination->at(1));

          // Skip triples that are shielded.
          if (graph.isAdjacentTo(a, c)) {
              continue;
          }

          if (graph.isDefCollider(a, b, c)) {
              continue;
          }

          if (sepsets->isCollider(a, b, c)) {
              if (!isArrowpointAllowed(a, b, graph)) {
                  continue;
              }

              if (!isArrowpointAllowed(c, b, graph)) {
                  continue;
              }

              graph.setEndpoint(a, b, ENDPOINT_ARROW);
              graph.setEndpoint(c, b, ENDPOINT_ARROW);

              /**FOR ANALYZING CONSISTENCY**/
              std::vector<Variable*> temp = sepsets->getSepset(a,c);
              if(temp.empty()) {
                  whyOrient.insert(std::pair<std::string, std::string>(a->getName() + "," + b->getName(), "0,null"));
                  whyOrient.insert(std::pair<std::string, std::string>(c->getName() + "," + b->getName(), "0,null"));
              }
              else
              {
                  std::string x = "0";
                  for(Variable* p:temp)
                  {
                      x+=(","+p->getName());
                  }
                  whyOrient.insert(std::pair<std::string, std::string>(a->getName()+","+b->getName(),x));
                  whyOrient.insert(std::pair<std::string, std::string>(c->getName()+","+b->getName(),x));
              }

      }
    }
  }
}

/**
 * Orients the graph according to rules in the graph (FCI step D).
 * <p>
 * Zhang's step F4, rules R1-R10.
 */
void FciOrient::doFinalOrientation(EdgeListGraph& graph) {
  if (completeRuleSetUsed) {
      zhangFinalOrientation(graph);
  } else {
      spirtesFinalOrientation(graph);
  }
}


//Does all 3 of these rules at once instead of going through all
// triples multiple times per iteration of doFinalOrientation.
void FciOrient::rulesR1R2cycle(EdgeListGraph& graph) {
  std::vector<Variable*> nodes = graph.getNodes();

  for (Variable* B : nodes) {
      std::vector<Variable*> adj = graph.getAdjacentNodes(B);

      if (adj.size() < 2) {
          continue;
      }

      ChoiceGenerator cg(adj.size(), 2);
      std::vector<int> *combination;

      for (combination = cg.next(); combination != NULL; combination = cg.next()) {
          Variable* A = adj.at(combination->at(0));
          Variable* C = adj.at(combination->at(1));

          //choice gen doesnt do diff orders, so must switch A & C around.
          ruleR1(A, B, C, graph);
          ruleR1(C, B, A, graph);
          ruleR2(A, B, C, graph);
          ruleR2(C, B, A, graph);
      }
  }
}


/**
 * Implements the double-triangle orientation rule, which states that if
 * D*-oB, A*->B<-*C and A*-oDo-*C, then
 * D*->B.
 * <p>
 * This is Zhang's rule R3.
 */
void FciOrient::ruleR3(EdgeListGraph& graph) {
  std::vector<Variable*> nodes = graph.getNodes();


  for (Variable* B : nodes) {

      std::vector<Variable*> intoBArrows = graph.getNodesInTo(B, ENDPOINT_ARROW);
      std::vector<Variable*> intoBCircles = graph.getNodesInTo(B, ENDPOINT_CIRCLE);

      for (Variable* D : intoBCircles) {
          if (intoBArrows.size() < 2) {
              continue;
          }

          ChoiceGenerator gen(intoBArrows.size(), 2);
          std::vector<int> *combination;

          for (combination = gen.next(); combination != NULL; combination = gen.next()) {
              Variable* A = intoBArrows.at(combination->at(0));
              Variable* C = intoBArrows.at(combination->at(1));

              if (graph.isAdjacentTo(A, C)) {
                  continue;
              }

              if (!graph.isAdjacentTo(A, D)
                      || !graph.isAdjacentTo(C, D)) {
                  continue;
              }

              if (!sepsets->isNoncollider(A, D, C)) {
                  continue;
              }

              if (graph.getEndpoint(A, D) != ENDPOINT_CIRCLE) {
                  continue;
              }

              if (graph.getEndpoint(C, D) != ENDPOINT_CIRCLE) {
                  continue;
              }

              if (!isArrowpointAllowed(D, B, graph)) {
                  continue;
              }

              graph.setEndpoint(D, B, ENDPOINT_ARROW);

              whyOrient.insert(std::pair<std::string, std::string>(D->getName() + "," + B->getName(),"3,"+A->getName()+","+C->getName()));

              changeFlag = true;
          }
      }
  }
}


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
void FciOrient::ruleR4B(EdgeListGraph& graph) {
  if (skipDiscriminatingPathRule) {
      return;
  }

  std::vector<Variable*> nodes = graph.getNodes();

  for (Variable* b : nodes) {

      //potential A and C candidate pairs are only those
      // that look like this:   A<-*Bo-*C
      std::vector<Variable*> possA = graph.getNodesOutTo(b, ENDPOINT_ARROW);
      std::vector<Variable*> possC = graph.getNodesInTo(b, ENDPOINT_CIRCLE);

      for (Variable* a : possA) {
          for (Variable* c : possC) {
              if (!graph.isParentOf(a, c)) {
                  continue;
              }

              if (graph.getEndpoint(b, c) != ENDPOINT_ARROW) {
                  continue;
              }

              ddpOrient(a, b, c, graph);
          }
      }
  }
}

/**
 * a method to search "back from a" to find a DDP. It is called with a
 * reachability list (first consisting only of a). This is breadth-first,
 * utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. </p>
 * The body of a DDP consists of colliders that are parents of c.
 */
void FciOrient::ddpOrient(Variable* a, Variable* b, Variable* c, EdgeListGraph& graph) {
  std::queue<Variable*> Q;

  std::unordered_set<Variable*> V;

  Variable* e = NULL;
  int distance = 0;

  std::unordered_map<Variable*, Variable*> previous;

  std::vector<Variable*> cParents = graph.getParents(c);

  Q.push(a);
  V.insert(a);
  V.insert(b);
  previous.insert(std::pair<Variable*,Variable*>(a, b));

  while (!Q.empty()) {
      Variable* t = Q.front();
      Q.pop();

      if (e == NULL || e == t) {
          e = t;
          distance++;
          if (distance > 0 && (distance > (maxPathLength == -1) ? 1000 : maxPathLength)) {
              return;
          }
      }

      const std::vector<Variable*> nodesInTo = graph.getNodesInTo(t, ENDPOINT_ARROW);

      for (Variable* d : nodesInTo) {
          if (std::count(V.begin(), V.end(), d) != 0) {
              continue;
          }

          previous.insert(std::pair<Variable*, Variable*>(d, t));
          Variable* p = previous.at(t);

          if (!graph.isDefCollider(d, t, p)) {
              continue;
          }

          previous.insert(std::pair<Variable*, Variable*>(d, t));

          if (!graph.isAdjacentTo(d, c)) {
              if (doDdpOrientation(d, a, b, c, previous, graph)) {
                  return;
              }
          }
          if (std::count(cParents.begin(), cParents.end(), d) != 0) {
              Q.push(d);
              V.insert(d);
          }
      }
  }
}

/**
 * Implements Zhang's rule R5, orient circle undirectedPaths: for any Ao-oB,
 * if there is an uncovered circle path u =
 * <A,C,...,D,B> such that A,D nonadjacent and B,C nonadjacent, then A---B
 * and orient every edge on u undirected.
 */
void FciOrient::ruleR5(EdgeListGraph& graph) {
  std::vector<Variable*> nodes = graph.getNodes();

  for (Variable* a : nodes) {
      std::vector<Variable*> adjacents = graph.getNodesInTo(a, ENDPOINT_CIRCLE);

      for (Variable* b : adjacents) {
          if (!(graph.getEndpoint(a, b) == ENDPOINT_CIRCLE)) {
              continue;
          }
          // We know Ao-oB.

          std::vector<std::vector<Variable*>> ucCirclePaths = getUcCirclePaths(a, b, graph);

          for (std::vector<Variable*> u : ucCirclePaths) {
              if (u.size() < 3) {
                  continue;
              }

              Variable* c = u.at(1);
              Variable* d = u.at(u.size() - 2);

              if (graph.isAdjacentTo(a, d)) {
                  continue;
              }
              if (graph.isAdjacentTo(b, c)) {
                  continue;
              }
              // We know u is as required: R5 applies!


              graph.setEndpoint(a, b, ENDPOINT_TAIL);
              graph.setEndpoint(b, a, ENDPOINT_TAIL);
              orientTailPath(u, graph);
              changeFlag = true;
          }
      }
  }
}

/**
 * Implements Zhang's rules R6 and R7, applies them over the graph once.
 * Orient single tails. R6: If A---Bo-*C then A---B--*C. R7: If A--oBo-*C
 * and A,C nonadjacent, then A--oB--*C
 */
void FciOrient::ruleR6R7(EdgeListGraph& graph) {
  std::vector<Variable*> nodes = graph.getNodes();

  for (Variable* b : nodes) {
      std::vector<Variable*> adjacents = graph.getAdjacentNodes(b);

      if (adjacents.size() < 2) {
          continue;
      }

      ChoiceGenerator cg(adjacents.size(), 2);
      std::vector<int>  *combination;

      for (combination = cg.next(); combination != NULL; combination = cg.next()) {
          Variable* a = adjacents.at(combination->at(0));
          Variable* c = adjacents.at(combination->at(1));

          if (graph.isAdjacentTo(a, c)) {
              continue;
          }

          if (!(graph.getEndpoint(b, a) == ENDPOINT_TAIL)) {
              continue;
          }
          if (!(graph.getEndpoint(c, b) == ENDPOINT_CIRCLE)) {
              continue;
          }
          // We know A--*Bo-*C.

          if (graph.getEndpoint(a, b) == ENDPOINT_TAIL) {

              // We know A---Bo-*C: R6 applies!
              graph.setEndpoint(c, b, ENDPOINT_TAIL);

              changeFlag = true;
          }

          if (graph.getEndpoint(a, b) == ENDPOINT_CIRCLE) {
//                    if (graph.isAdjacentTo(a, c)) continue;


              // We know A--oBo-*C and A,C nonadjacent: R7 applies!
              graph.setEndpoint(c, b, ENDPOINT_TAIL);
              changeFlag = true;
          }

      }
  }
}

/**
 * Implements Zhang's rules R8, R9, R10, applies them over the graph once.
 * Orient arrow tails. I.e., tries R8, R9, and R10 in that sequence on each
 * Ao->C in the graph.
 */
void FciOrient::rulesR8R9R10(EdgeListGraph& graph) {
  std::vector<Variable*> nodes = graph.getNodes();

  for (Variable* c : nodes) {
      std::vector<Variable*> intoCArrows = graph.getNodesInTo(c, ENDPOINT_ARROW);

      for (Variable* a : intoCArrows) {
          if (!(graph.getEndpoint(c, a) == ENDPOINT_CIRCLE)) {
              continue;
          }
          // We know Ao->C.

          // Try each of R8, R9, R10 in that order, stopping ASAP.
          if (!ruleR8(a, c, graph)) {
              bool b = ruleR9(a, c, graph);

              if (!b) {
                  ruleR10(a, c, graph);
              }
          }
      }
  }
}

/**
 * @param maxPathLength the maximum length of any discriminating path, or -1
 * if unlimited.
 */
void FciOrient::setMaxPathLength(int maxPathLength) {
  if (maxPathLength < -1) {
      throw std::invalid_argument("Max path length must be -1 (unlimited) or >= 0: " + std::to_string(maxPathLength));
  }

  this->maxPathLength = maxPathLength;
}

//============================PRIVATE METHODS============================//

void FciOrient::printWrongColliderMessage(Variable* a, Variable* b, Variable* c, std::string location, EdgeListGraph& graph) {
  // if (&truePag != NULL && graph.isDefCollider(a, b, c) && !truePag.isDefCollider(a, b, c)) {
  if (truePag.getNumNodes() >= 1 && graph.isDefCollider(a, b, c) && !truePag.isDefCollider(a, b, c)) {
      Rcpp::Rcout << location << ": Orienting collider by mistake: " << a << "*->" << b << "<-*" << c;
  }
}

void FciOrient::spirtesFinalOrientation(EdgeListGraph& graph) {
  changeFlag = true;
  bool firstTime = true;

  while (changeFlag) {
      changeFlag = false;
      rulesR1R2cycle(graph);
      ruleR3(graph);

      // R4 requires an arrow orientation.
      if (changeFlag || firstTime /* && !knowledge.isEmpty() */) {
          ruleR4B(graph);
          firstTime = false;
      }
  }
}

void FciOrient::zhangFinalOrientation(EdgeListGraph& graph) {
  changeFlag = true;
  bool firstTime = true;

  while (changeFlag) {
      changeFlag = false;
      rulesR1R2cycle(graph);
      ruleR3(graph);

      // R4 requires an arrow orientation.
      if (changeFlag || firstTime /* && !knowledge.isEmpty() */ ) {
          ruleR4B(graph);
          firstTime = false;
      }
  }

  if (isCompleteRuleSetUsed()) {
      // Now, by a remark on page 100 of Zhang's dissertation, we apply rule
      // R5 once.
      ruleR5(graph);

      // Now, by a further remark on page 102, we apply R6,R7 as many times
      // as possible.
      changeFlag = true;

      while (changeFlag) {
          changeFlag = false;
          ruleR6R7(graph);
      }

      // Finally, we apply R8-R10 as many times as possible.
      changeFlag = true;

      while (changeFlag) {
          changeFlag = false;
          rulesR8R9R10(graph);
      }

  }
}

 /// R1, away from collider
 // If a*->bo-*c and a, c not adjacent then a*->b->c
void FciOrient::ruleR1(Variable* a, Variable* b, Variable* c, EdgeListGraph& graph) {
  if (graph.isAdjacentTo(a, c)) {
      return;
  }

  if (graph.getEndpoint(a, b) == ENDPOINT_ARROW && graph.getEndpoint(c, b) == ENDPOINT_CIRCLE) {
      if (!isArrowpointAllowed(b, c, graph)) {
          return;
      }

      graph.setEndpoint(c, b, ENDPOINT_TAIL);
      graph.setEndpoint(b, c, ENDPOINT_ARROW);
      whyOrient.insert(std::pair<std::string, std::string>(b->getName() + "," + c->getName(),"1"));
      changeFlag = true;
  }
}


 //if a*-oc and either a-->b*->c or a*->b-->c, then a*->c
 // This is Zhang's rule R2.
void FciOrient::ruleR2(Variable* a, Variable* b, Variable* c, EdgeListGraph& graph) {
  if ((graph.isAdjacentTo(a, c))
          && (graph.getEndpoint(a, c) == ENDPOINT_CIRCLE)) {

      if ((graph.getEndpoint(a, b) == ENDPOINT_ARROW)
              && (graph.getEndpoint(b, c) == ENDPOINT_ARROW) && ((graph.getEndpoint(b, a) == ENDPOINT_TAIL)
              || (graph.getEndpoint(c, b) == ENDPOINT_TAIL))) {

          if (!isArrowpointAllowed(a, c, graph)) {
              return;
          }

          graph.setEndpoint(a, c, ENDPOINT_ARROW);
          whyOrient.insert(std::pair<std::string, std::string>(a->getName() + "," + c->getName(), "2," + b->getName()));

          changeFlag = true;
      }
  }
}

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
void FciOrient::ruleR4A(EdgeListGraph& graph) {
  std::vector<Variable*> nodes = graph.getNodes();

  for (Variable* b : nodes) {

      //potential A and C candidate pairs are only those
      // that look like this:   A<-*Bo-*C
      std::vector<Variable*> possA = graph.getNodesOutTo(b, ENDPOINT_ARROW);
      std::vector<Variable*> possC = graph.getNodesInTo(b, ENDPOINT_CIRCLE);

      for (Variable* a : possA) {
          for (Variable* c : possC) {
              if (!graph.isParentOf(a, c)) {
                  continue;
              }

              if (graph.getEndpoint(b, c) != ENDPOINT_ARROW) {
                  continue;
              }

              std::list<Variable*> reachable;
              reachable.push_back(a);
          }
      }
  }
}

/**
 * a method to search "back from a" to find a DDP. It is called with a
 * reachability list (first consisting only of a). This is breadth-first,
 * utilizing "reachability" concept from Geiger, Verma, and Pearl 1990. </p>
 * The body of a DDP consists of colliders that are parents of c.
 */
void FciOrient::reachablePathFind(Variable* a, Variable* b, Variable* c,
                                  std::list<Variable*> reachable, EdgeListGraph& graph) {
//        Map<Node, Node> next = new HashMap<Node, Node>();   // RFCI: stores the next node in the disciminating path
//        // path containing the nodes in the traiangle
//        next.put(a, b);
//        next.put(b, c);



  std::vector<Variable*> v = graph.getParents(c);
  std::unordered_set<Variable*> cParents(v.begin(), v.end());


  // Needed to avoid cycles in failure case.
  std::unordered_set<Variable*> visited;
  visited.insert(b);
  visited.insert(c);

  Variable* e = reachable.front();
  int distance = 0;

  // We don't want to include a,b,or c on the path, so they are added to
  // the "visited" set.  b and c are added explicitly here; a will be
  // added in the first while iteration.
  while (reachable.size() > 0) {
      Variable* x = reachable.front();
      reachable.pop_front();
      visited.insert(x);

      if (e == x) {
          e = x;
          distance++;
          const int maxPathLength_ = (maxPathLength == -1) ? 1000 : maxPathLength;

          if (distance > 0 && distance > maxPathLength_) {
              continue;
          }
      }

      // Possible DDP path endpoints.
      std::vector<Variable*> pathExtensions = graph.getNodesInTo(x, ENDPOINT_ARROW);
      for (Variable* var: visited) {
        std::remove(pathExtensions.begin(), pathExtensions.end(), var);
      }

      for (Variable* d : pathExtensions) {
          // If d is reachable and not adjacent to c, its a DDP
          // endpoint, so do DDP orientation. Otherwise, if d <-> c,
          // add d to the list of reachable nodes.
          if (!graph.isAdjacentTo(d, c)) {
              // Check whether <a, b, c> should be reoriented given
              // that d is not adjacent to c; if so, orient and stop.
              doDdpOrientation(d, a, b, c, graph);
              return;
          } else if (cParents.count(d)) {
              if (graph.getEndpoint(x, d) == ENDPOINT_ARROW) {
                  reachable.push_back(d);
              }
          }
      }
  }
}

/**
 * Orients the edges inside the definte discriminating path triangle. Takes
 * the left endpoint, and a,b,c as arguments.
 */
void FciOrient::doDdpOrientation(Variable* d, Variable* a, Variable* b, Variable* c, EdgeListGraph& graph) {

  std::vector<Variable*> sepset = getSepset(d, c);

  if (sepset.empty()) {
      return;
  }

  if (std::count(sepset.begin(), sepset.end(), b) != 0) {
      graph.setEndpoint(c, b, ENDPOINT_TAIL);
      changeFlag = true;

  } else {
      if (!isArrowpointAllowed(a, b, graph)) {
          return;
      }

      if (!isArrowpointAllowed(c, b, graph)) {
          return;
      }

      graph.setEndpoint(a, b, ENDPOINT_ARROW);
      graph.setEndpoint(c, b, ENDPOINT_ARROW);
      changeFlag = true;
  }
}

/**
 * Orients the edges inside the definte discriminating path triangle. Takes
 * the left endpoint, and a,b,c as arguments.
 */
bool FciOrient::doDdpOrientation(Variable* d, Variable* a, Variable* b, Variable* c,
                                 std::unordered_map<Variable*, Variable*> previous, EdgeListGraph& graph) {

  if (dag.getNumNodes() != 0) {
      if (dag.isAncestorOf(b, c)) {
          graph.setEndpoint(c, b, ENDPOINT_TAIL);
          whyOrient.insert(std::pair<std::string, std::string>(c->getName() + "," + b->getName(),"4"));
          changeFlag = true;
      } else {
          if (!isArrowpointAllowed(a, b, graph)) {
              return false;
          }

          if (!isArrowpointAllowed(c, b, graph)) {
              return false;
          }

          graph.setEndpoint(a, b, ENDPOINT_ARROW);
          graph.setEndpoint(c, b, ENDPOINT_ARROW);
          whyOrient.insert(std::pair<std::string, std::string>(a->getName() + "," + b->getName(),"4"));
          whyOrient.insert(std::pair<std::string, std::string>(c->getName() + "," + b->getName(),"4"));

          changeFlag = true;
      }

      return true;
  }

  if (graph.isAdjacentTo(d, c)) {
      throw std::invalid_argument("Illegal Argument");
  }

  std::vector<Variable*> path = getPath(d, previous);

  bool ind = getSepsets()->isIndependent(d, c, path);

  // std::vector<Variable*> path2 (path);
  // path2.remove(b);
  std::vector<Variable*> path2(path);
  if (std::find(path2.begin(), path2.end(), b) != path2.end()) {
    path2.erase(std::find(path2.begin(), path2.end(), b));
  }
  for(Variable* var: path) {
    if (var != b) {
      path2.push_back(var);
    }
  }


  bool ind2 = getSepsets()->isIndependent(d, c, path2);

  if (!ind && !ind2) {
      std::vector<Variable*> sepset = getSepsets()->getSepset(d, c);

      if (sepset.empty()) {
          return false;
      }

      ind = (std::count(sepset.begin(), sepset.end(), b));
  }

  //        printDdp(d, path, a, b, c, graph);
  if (ind) {
  //            if (sepset.contains(b)) {
      graph.setEndpoint(c, b, ENDPOINT_TAIL);

      changeFlag = true;
      return true;
  } else {
      if (!isArrowpointAllowed(a, b, graph)) {
          return false;
      }

      if (!isArrowpointAllowed(c, b, graph)) {
          return false;
      }

      graph.setEndpoint(a, b, ENDPOINT_ARROW);
      graph.setEndpoint(c, b, ENDPOINT_ARROW);

      changeFlag = true;
      return true;
  }


}

// graph is only used in the deleted verbose section, making it unnecessary here.
void FciOrient::printDdp(Variable* d, std::vector<Variable*> path, Variable* a, Variable* b, Variable* c, EdgeListGraph& graph) {
  std::vector<Variable*> nodes;
  nodes.push_back(d);
  for(Variable* var: path) {
    nodes.push_back(var);
  }
  nodes.push_back(a);
  nodes.push_back(b);
  nodes.push_back(c);
}


std::vector<Variable*> FciOrient::getPath(Variable* c, std::unordered_map<Variable*, Variable*> previous) {
  std::vector<Variable*> l;

  Variable* p = c;

  do {
      p = previous.at(p);

      if (p != NULL) {
          l.push_back(p);
      }
  } while (p != NULL);

  return l;
}

/**
 * Orients every edge on a path as undirected (i.e. A---B).
 * <p>
 * DOES NOT CHECK IF SUCH EDGES ACTUALLY EXIST: MAY DO WEIRD THINGS IF
 * PASSED AN ARBITRARY LIST OF NODES THAT IS NOT A PATH.
 *
 * @param path The path to orient as all tails.
 */
void FciOrient::orientTailPath(std::vector<Variable*> path, EdgeListGraph& graph) {
  for (int i = 0; i < path.size() - 1; i++) {
      Variable* n1 = path.at(i);
      Variable* n2 = path.at(i + 1);

      graph.setEndpoint(n1, n2, ENDPOINT_TAIL);
      graph.setEndpoint(n2, n1, ENDPOINT_TAIL);
      changeFlag = true;
  }
}

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
std::vector<std::vector<Variable*>> FciOrient::getUcPdPaths(Variable* n1, Variable* n2, EdgeListGraph& graph) {
  std::vector<std::vector<Variable*>> ucPdPaths;

  std::vector<Variable*> soFar;

  soFar.push_back(n1);

  std::vector<Variable*> adjacencies = graph.getAdjacentNodes(n1);
  for (Variable* curr : adjacencies) {
      getUcPdPsHelper(curr, soFar, n2, ucPdPaths, graph);
  }

  return ucPdPaths;
}

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
void FciOrient::getUcPdPsHelper(Variable* curr, std::vector<Variable*> soFar, Variable* end,
                                std::vector<std::vector<Variable*>> ucPdPaths, EdgeListGraph& graph) {

  if (std::count(soFar.begin(), soFar.end(), curr) != 0) {
      return;
  }

  Variable* prev = soFar.at(soFar.size() - 1);
  if (graph.getEndpoint(prev, curr) == ENDPOINT_TAIL
          || graph.getEndpoint(curr, prev) == ENDPOINT_ARROW) {
      return; // Adding curr would make soFar not p.d.
  } else if (soFar.size() >= 2) {
      Variable* prev2 = soFar.at(soFar.size() - 2);
      if (graph.isAdjacentTo(prev2, curr)) {
          return; // Adding curr would make soFar not uncovered.
      }
  }

  soFar.push_back(curr); // Adding curr is OK, so let's do it.

  if (curr == end) {
      // We've reached the goal! Save soFar as a path.
      std::vector<Variable*> soFar_copy(soFar);
      ucPdPaths.push_back(soFar_copy);
  } else {
      // Otherwise, try each node adjacent to the getModel one.
      std::vector<Variable*> adjacents = graph.getAdjacentNodes(curr);
      for (Variable* next : adjacents) {
          getUcPdPsHelper(next, soFar, end, ucPdPaths, graph);
      }
  }
  soFar.pop_back();
}


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
std::vector<std::vector<Variable*>> FciOrient::getUcCirclePaths(Variable* n1, Variable* n2, EdgeListGraph& graph) {
  std::vector<std::vector<Variable*>> ucCirclePaths;
  std::vector<std::vector<Variable*>> ucPdPaths = getUcPdPaths(n1, n2, graph);

  for (std::vector<Variable*> path : ucPdPaths) {
      for (int i = 0; i < path.size() - 1; i++) {
          Variable* j = path.at(i);
          Variable* sj = path.at(i + 1);

          if (!(graph.getEndpoint(j, sj) == ENDPOINT_CIRCLE)) {
              break;
          }
          if (!(graph.getEndpoint(sj, j) == ENDPOINT_CIRCLE)) {
              break;
          }
          // This edge is OK, it's all circles.

          if (i == path.size() - 2) {
              // We're at the last edge, so this is a circle path.
              ucCirclePaths.push_back(path);
          }
      }
  }

  return ucCirclePaths;
}

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
bool FciOrient::ruleR8(Variable* a, Variable* c, EdgeListGraph& graph) {

  std::vector<Variable*> intoCArrows = graph.getNodesInTo(c, ENDPOINT_ARROW);

  for (Variable* b : intoCArrows) {
      // We have B*->C.
      if (!graph.isAdjacentTo(a, b)) {
          continue;
      }
      if (!graph.isAdjacentTo(b, c)) {
          continue;
      }

      // We have A*-*B*->C.
      if (!(graph.getEndpoint(b, a) == ENDPOINT_TAIL)) {
          continue;
      }
      if (!(graph.getEndpoint(c, b) == ENDPOINT_TAIL)) {
          continue;
      }
      // We have A--*B-->C.

      if (graph.getEndpoint(a, b) == ENDPOINT_TAIL) {
          continue;
      }
      // We have A-->B-->C or A--oB-->C: R8 applies!


      graph.setEndpoint(c, a, ENDPOINT_TAIL);
      changeFlag = true;
      return true;
  }

  return false;
}

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
bool FciOrient::ruleR9(Variable* a, Variable* c, EdgeListGraph& graph) {
  std::vector<std::vector<Variable*>> ucPdPsToC = getUcPdPaths(a, c, graph);

  for (std::vector<Variable*> u : ucPdPsToC) {
      Variable* b = u.at(1);
      if (graph.isAdjacentTo(b, c)) {
          continue;
      }
      if (b == c) {
          continue;
      }
      // We know u is as required: R9 applies!


      graph.setEndpoint(c, a, ENDPOINT_TAIL);
      changeFlag = true;
      return true;
  }

  return false;
}

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
bool FciOrient::ruleR10(Variable* a, Variable* c, EdgeListGraph& graph) {
  std::vector<Variable*> intoCArrows = graph.getNodesInTo(c, ENDPOINT_ARROW);

  for (Variable* b : intoCArrows) {
      if (b == a) {
          continue;
      }

      if (!(graph.getEndpoint(c, b) == ENDPOINT_TAIL)) {
          continue;
      }
      // We know Ao->C and B-->C.

      for (Variable* d : intoCArrows) {
          if (d == a || d == b) {
              continue;
          }

          if (!(graph.getEndpoint(d, c) == ENDPOINT_TAIL)) {
              continue;
          }
          // We know Ao->C and B-->C<--D.

          std::vector<std::vector<Variable*>> ucPdPsToB = getUcPdPaths(a, b, graph);
          std::vector<std::vector<Variable*>> ucPdPsToD = getUcPdPaths(a, d, graph);
          for (std::vector<Variable*> u1 : ucPdPsToB) {
              Variable* m = u1.at(1);
              for (std::vector<Variable*> u2 : ucPdPsToD) {
                  Variable* n = u2.at(1);

                  if (m == n) {
                      continue;
                  }
                  if (graph.isAdjacentTo(m, n)) {
                      continue;
                  }
                  // We know B,D,u1,u2 as required: R10 applies!


                  graph.setEndpoint(c, a, ENDPOINT_TAIL);
                  changeFlag = true;
                  return true;
              }
          }
      }
  }

  return false;
}



/**
 * Helper method. Appears to check if an arrowpoint is permitted by
 * background knowledge.
 *
 * @param x The possible other node.
 * @param y The possible point node.
 * @return Whether the arrowpoint is allowed.
 */
bool FciOrient::isArrowpointAllowed(Variable* x, Variable* y, EdgeListGraph& graph) {
  if (graph.getEndpoint(x, y) == ENDPOINT_ARROW) {
      return true;
  }

  if (graph.getEndpoint(x, y) == ENDPOINT_TAIL) {
      return false;
  }

  if (graph.getEndpoint(y, x) == ENDPOINT_ARROW) {
      // switched the commented sections
      return true;
      // if (!knowledge.isForbidden(x.getName(), y.getName())) {
      //     return true;
      // }
  }

  if (graph.getEndpoint(y, x) == ENDPOINT_TAIL) {
      // added return true as did above in place of knowledge
      return true;
      // if (!knowledge.isForbidden(x.getName(), y.getName())) {
      //     return true;
      // }
  }

  return graph.getEndpoint(y, x) == ENDPOINT_CIRCLE;
}
