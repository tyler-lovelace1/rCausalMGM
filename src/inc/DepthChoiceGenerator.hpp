#ifndef DEPTHCHOICEGENERATOR_HPP_
#define DEPTHCHOICEGENERATOR_HPP_

#include "armaLapack.hpp"

/**
 * Generates (nonrecursively) all of the combinations of a choose b, where a, b
 * are positive integers and a >= b.  The values of a and b are given in the
 * constructor, and the sequence of choices is obtained by repeatedly calling
 * the next() method.  When the sequence is finished, NULL is returned.</p> </p>
 * <p>A valid combination for the sequence of combinations for a choose b
 * generated by this class is an array x[] of b integers i, 0 <= i < a, such
 * that x[j] < x[j + 1] for each j from 0 to b - 1.
 * <p>
 * Works by calling ChoiceGenerator with increasingly larger values of a.
 * <p>
 * To see what this class does, try calling DepthChoiceGenerator.testPrint(5, 3), for
 * instance.
 *
 * @author Joseph Ramsey
 */
class DepthChoiceGenerator {

private:
    /**
     * The number of objects being selected from.
     */
    int a;

    /**
     * The number of objects in the desired selection.
     */
    int b;

    /**
     * The difference between a and b (should be nonnegative).
     */
    int diff;

    /**
     * Maximum a.
     */
    int depth;

    /**
     * Effective maximum a.
     */
    int effectiveDepth;

    /**
     * The internally stored choice.
     */
    std::vector<int> choiceLocal;

    
    /**
     * The choice that is returned. Used, since the returned array can be
     * modified by the user.
     */
    std::vector<int> choiceReturned;

    /**
     * Indicates whether the next() method has been called since the last
     * initialization.
     */
    bool begun;

    /**
     * Fills the 'choice' array, from index 'index' to the end of the array,
     * with successive integers starting with choice[index] + 1.
     *
     * @param index the index to begin this incrementing operation.
     */
    void fill(int index);

    void initialize();

public:

    /**
     * Constructs a new choice generator for a choose b. Once this
     * initialization has been performed, successive calls to next() will
     * produce the series of combinations.  To begin a new series at any time,
     * call this init method again with new values for a and b.
     *
     * @param a     the number of objects being selected from.
     * @param depth the maximum number of objects selected.
     */
    DepthChoiceGenerator(int a, int depth);

    /**
     * @return the next combination in the series, or NULL if the series is
     * finished.
     */
    std::vector<int>* next();

    /**
     * This static method will print the series of combinations for a choose depth
     * to System.out.
     *
     * @param a     the number of objects being selected from.
     * @param depth the maximum number of objects selected.
     */
    static void testPrint(int a, int depth);

    int getA() { return a; }
    int getB() { return b; }
    int getDepth() { return depth; }
};

#endif /* DEPTHCHOICEGENERATOR_HPP_ */
