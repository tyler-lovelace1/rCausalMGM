#include "IndependenceTestRandom.hpp"

/**
 * @return true about 25% of the time
 */
bool IndependenceTestRandom::isIndependent(Variable* x, Variable* y, std::vector<Variable*>& z) {
    uint seed = 13 * std::hash<std::string>()(x->getName()) + 17 * std::hash<std::string>()(y->getName());

    for (Variable* node : z) {
        seed += 19 * std::hash<std::string>()(node->getName());
    }

    std::srand(seed);

    bool independent = (std::rand() % 4 == 0);

    if (independent) {
        Rcpp::Rcout << "Variable " << x->getName() << " is independent of " << y->getName() << " conditioned on [";
        for (Variable* node : z) {
            Rcpp::Rcout << node->getName() << " - ";
        } 
        Rcpp::Rcout << "]" << std::endl;
    }

    return independent;
}
