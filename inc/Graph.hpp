#ifndef GRAPH_HPP_
#define GRAPH_HPP_

// [[Rcpp::depends(BH)]]

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "Variable.hpp"

typedef boost::adjacency_list<boost::setS, 
                              boost::vecS, 
                              boost::bidirectionalS,
                              Variable*,
                              boost::no_property> Graph;


#endif /* GRAPH_HPP_ */

