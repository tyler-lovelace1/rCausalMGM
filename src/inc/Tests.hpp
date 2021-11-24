#ifndef TESTS_HPP
#define TESTS_HPP

#include <RcppArmadillo.h>

class Tests {

public:
    static void testConcurrentQueue();

    static void testMGMFunctions(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testMGMTiming(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testPcStable(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testCpcStable(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testPcMax(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testSTEPS(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testGraphFromFile(const Rcpp::DataFrame &df, const std::string& filename, const int maxDiscrete = 5);

};

#endif /* TESTS_HPP */