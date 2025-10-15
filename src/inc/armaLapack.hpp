/*
 * Simple header file including RcppArmadillo headers and restoring
 * the macro NDEBUG if it is unset by RcppArmadillo (see comment below)
 *
 * @author Alain Hauser
 * $Id: armaLapack.hpp 13 2011-10-31 16:49:25Z alhauser $
 */

#ifndef ARMALAPACK_HPP_
#define ARMALAPACK_HPP_

#ifdef NDEBUG
// Make sure Armadillo does use "assert"...
#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// ... but restore the NDEBUG macro (unset by RcppArmadillo.h)
#define NDEBUG
#else
// Simply include RcppArmadillo
#include <RcppArmadillo.h>
#endif /* NDEBUG */

// For compiler that do not know uint...
typedef unsigned int uint;

// static int randWrapper(const int n) { return std::floor(R::unif_rand()*n); }

struct R_RNG_Engine {
    typedef uint result_type;
    
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

    result_type operator()() {
        return static_cast<result_type>(R::unif_rand() * max());
    }
};

#endif /* ARMALAPACK_HPP_ */
