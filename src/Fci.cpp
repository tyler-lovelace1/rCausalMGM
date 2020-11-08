#include "Fci.hpp"

//============================CONSTRUCTORS============================//

/**
 * Constructs a new FCI search for the given independence test and background knowledge.
 */
Fci::Fci(IndependenceTest independenceTest) {

}

/**
 * Constructs a new FCI search for the given independence test and background knowledge and a list of variables to
 * search over.
 */
Fci::Fci(IndependenceTest independenceTest, std::vector<Variable*> searchVars) {

}

//========================PUBLIC METHODS==========================//

EdgeListGraph Fci::search(IFas fas) {
  EdgeListGraph temp;
  return temp;
}

/**
 * @param maxPathLength the maximum length of any discriminating path, or -1 if unlimited.
 */
void Fci::setMaxPathLength(int maxPathLength) {
  return;
}

//========================PRIVATE METHODS==========================//

void Fci::buildIndexing(std::vector<Variable*> nodes) {
  return;
}

/**
 * Orients according to background knowledge
 */
void fciOrientbk(/*IKnowledge bk,*/ EdgeListGraph graph, std::vector<Variable*> variables) {
  return;
}
