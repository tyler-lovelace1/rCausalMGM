#include "MGM.hpp"

MGM::MGM(arma::mat& x, arma::mat& y, std::vector<Variable*>& variables, std::vector<int>& l, std::vector<double>& lambda) {
    
    if (l.size() != y.n_cols)
        throw std::invalid_argument("length of l doesn't match number of variables in Y");

    if (y.n_rows != x.n_rows)
        throw std::invalid_argument("different number of samples for x and y");

    //lambda should have 3 values corresponding to cc, cd, and dd
    if (lambda.size() != 3)
        throw std::invalid_argument("Lambda should have three values for cc, cd, and dd edges respectively");
    
    this->xDat = x;
    this->yDat = y;
    this->l = l;
    this->p = x.n_cols;
    this->q = y.n_cols;
    this->n = x.n_rows;
    this->variables = variables;

    this->lambda = arma::vec(lambda);
    fixData();
    initParameters();
    calcWeights();
    makeDummy();

}
