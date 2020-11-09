#include "DepthChoiceGenerator.hpp"

/**
 * Constructs a new choice generator for a choose b. Once this
 * initialization has been performed, successive calls to next() will
 * produce the series of combinations.  To begin a new series at any time,
 * call this init method again with new values for a and b.
 *
 * @param a     the number of objects being selected from.
 * @param depth the maximum number of objects selected.
 */
DepthChoiceGenerator::DepthChoiceGenerator(int a, int depth) {
    if ((a < 0) || (depth < -1))
        throw std::invalid_argument("For 'a choose depth', a must be nonnegative and depth >= -1: a = " + std::to_string(a) + " depth = " + std::to_string(depth));

    this->a = a;
    this->b = 0;

    this->depth = depth;
    effectiveDepth = depth;
    if ((depth == -1) || (depth > a)) effectiveDepth = a;

    initialize();

    
}

void DepthChoiceGenerator::initialize() {
    choiceLocal = std::vector<int>(b);
    choiceReturned = std::vector<int>(b);
    diff = a-b;

    // Initialize the choice array with successive integers [0 1 2 ...].
    // Set the value at the last index one less than it would be in such
    // a series, ([0 1 2 ... b - 2]) so that on the first call to next()
    // the first combination ([0 1 2 ... b - 1]) is returned correctly.

    if (b > 0) {
        for (int i = 0; i < b-1; i++) {
            choiceLocal[i] = i;
        }

        choiceLocal[b-1] = b-2;
    }
        
    begun = false;
}

/**
 * @return the next combination in the series, or NULL if the series is
 * finished.
 */
std::vector<int>* DepthChoiceGenerator::next() {
    int i = b;

    // Scan from the right for the first index whose value is less than
    // its expected maximum (i + diff) and perform the fill() operation
    // at that index.
    while(--i > -1) {
        if (choiceLocal[i] < i + diff) {
            fill(i);
            begun = true;
            choiceReturned = choiceLocal;
            return &choiceReturned;
        }
    }

    if (begun) {
        b++;

        if (b > effectiveDepth) return NULL;

        initialize();
        return next();
    } else {
        begun = true;
        choiceReturned = choiceLocal;
        return &choiceReturned;
    }
}

/**
 * Fills the 'choice' array, from index 'index' to the end of the array,
 * with successive integers starting with choice[index] + 1.
 *
 * @param index the index to begin this incrementing operation.
 */
void DepthChoiceGenerator::fill(int index) {
    choiceLocal[index]++;

    for (int i = index + 1; i < b; i++) {
        choiceLocal[i] = choiceLocal[i-1] + 1;
    }
}

/**
 * This static method will print the series of combinations for a choose depth
 * to System.out.
 *
 * @param a     the number of objects being selected from.
 * @param depth the maximum number of objects selected.
 */
void DepthChoiceGenerator::testPrint(int a, int depth) {
    DepthChoiceGenerator cg(a, depth);

    std::vector<int>* choice;

    Rcpp::Rcout << "\nPrinting combinations for " << a << " choose (depth) " << depth << ":" << std::endl;

    while ((choice = cg.next()) != NULL) {
        for (int aChoice : *choice) {
            Rcpp::Rcout << aChoice << "\t";
        }

        Rcpp::Rcout << std::endl;
    }

    Rcpp::Rcout << std::endl;
}