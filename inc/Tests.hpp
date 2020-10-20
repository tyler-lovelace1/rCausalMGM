#ifndef TESTS_HPP
#define TESTS_HPP

#include "MGM.hpp"
#include "IndTestMulti.hpp"
#include "SepsetMap.hpp"
#include "ChoiceGenerator.hpp"
#include "IndependenceTestRandom.hpp"
#include "PcStable.hpp"
#include "CpcStable.hpp"
#include "PcMax.hpp"
#include "BlockingQueue.hpp"
#include "STEPS.hpp"
#include <thread>
#include <atomic>
#include <cstdlib>


class Tests {

public:
    static void testConcurrentQueue();

    static void testMGMFunctions(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testMGMTiming(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testPcStable(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testCpcStable(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testPcMax(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

    static void testSTEPS(const Rcpp::DataFrame &df, const int maxDiscrete = 5);

};

#endif /* TESTS_HPP */