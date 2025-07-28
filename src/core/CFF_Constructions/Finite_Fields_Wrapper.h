#ifndef FFCFF_WRAPPER_HEADER
#define FFCFF_WRAPPER_HEADER

#include <math.h>

// to remove the FLINT library dependancy from this codebase, reimplement these functions

void populateFiniteField(int p, int k, int *addition_field, int *multiplication_field);
int compute_field_size(int p, int k);
int hornerPolynomialEvalOverFq(
    int polyLength,
    int polynomialCoefficients[polyLength],
    int x,
    int q,
    int addition_field[q][q],
    int multiplication_field[q][q]
);
void populateMultiplicativeInverses(int p, int k, int *multiplication_field, int *multiplicative_inverses);
void populateAdditiveInverses(int p, int k, int *addition_field, int *additive_inverses);
#endif