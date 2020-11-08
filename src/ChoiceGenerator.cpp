#include "ChoiceGenerator.hpp"

/**
 * Constructs a new choice generator for a choose b. Once this
 * initialization has been performed, successive calls to next() will
 * produce the series of combinations.  To begin a new series at any time,
 * call this init method again with new values for a and b.
 *
 * @param a the number of objects being selected from.
 * @param b the number of objects in the desired selection.
 */
ChoiceGenerator::ChoiceGenerator(int a, int b) {
    if ((a < 0) || (b < 0) || (a < b))
        throw std::invalid_argument("For 'a choose b', a and b must be nonnegative with a >= b: a = " + std::to_string(a) + " b = " + std::to_string(b));

    this->a = a;
    this->b = b;
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
std::vector<int>* ChoiceGenerator::next() {
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
        return NULL;
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
void ChoiceGenerator::fill(int index) {
    choiceLocal[index]++;

    for (int i = index + 1; i < b; i++) {
        choiceLocal[i] = choiceLocal[i-1] + 1;
    }
}

/**
 * This static method will print the series of combinations for a choose b
 * to System.out.
 *
 * @param a the number of objects being selected from.
 * @param b the number of objects in the desired selection.
 */
void ChoiceGenerator::testPrint(int a, int b) {
    ChoiceGenerator cg(a, b);

    std::vector<int>* choice;

    Rcpp::Rcout << "\nPrinting combinations for " << a << " choose " << b << ":" << std::endl;

    while ((choice = cg.next()) != NULL) {
        for (int aChoice : *choice) {
            Rcpp::Rcout << aChoice << "\t";
        }

        Rcpp::Rcout << std::endl;
    }

    Rcpp::Rcout << std::endl;
}