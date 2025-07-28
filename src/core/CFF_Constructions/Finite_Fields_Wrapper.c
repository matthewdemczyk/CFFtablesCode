#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include "flint/flint.h"
#include "flint/fq.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_mod_poly.h"

#include "../CombinatorialAlgorithms.h"


// largest power allowed (a field is size prime ^ k)
#define MAX_K 32

typedef struct
{
    fmpz_poly_t poly;
    fq_t ff;
} field_element_t;

int compute_field_size(int p, int k)
{
    int result = 1;
    for(int i = 0; i < k; i++)
    {
        if (result > INT_MAX / p)
        {
            return -1;
        }
        result *= p;
    }
    return result;
}

static inline int convertPolynomialToInteger(int base, int k, const int* number) {
    int ret = number[k - 1];
    for (int i = k - 2; i >= 0; i--)
    {
        ret = ret * base + number[i];
    }
    return ret;
}

static inline int convertPolynomialToInteger_opt(int base, int k, const int* number) {
    // handle small k values directly to avoid loop overhead
    switch (k) {
        case 1:
            return number[0];
        case 2:
            return number[1] * base + number[0];
        case 3:
            return (number[2] * base + number[1]) * base + number[0];
        case 4:
            return ((number[3] * base + number[2]) * base + number[1]) * base + number[0];
    }

    // for larger k, use horner's method
    const int *ptr = number + k - 1;
    int ret = *ptr--;

    // unroll the loop by 4
    while (ptr > number + 3) {
        ret = ret * base + *ptr--;
        ret = ret * base + *ptr--;
        ret = ret * base + *ptr--;
        ret = ret * base + *ptr--;
    }

    // handle remaining elements
    while (ptr >= number) {
        ret = ret * base + *ptr--;
    }

    return ret;
}

// helper function to used in populateFiniteField
static inline int* get_element(int *array, int row, int col, int field_size) {
    return &array[row * field_size + col];
}

void populateFiniteField(int p, int k, int *addition_field, int *multiplication_field) {
    if (k > MAX_K || k <= 0 || p <= 1) {
        //return -1;
        exit(1); // invalid parameters
    }

    int field_size = compute_field_size(p, k);
    if (field_size <= 0) {
        //return -1;
        exit(1); // field size computation failed
    }

    // initialize context
    fq_ctx_t ctx;
    fmpz_t p_fmpz;
    fmpz_init_set_ui(p_fmpz, p);
    fq_ctx_init(ctx, p_fmpz, k, "c");

    // pre-allocate working elements
    field_element_t *elements = malloc(field_size * sizeof(field_element_t));
    if (!elements) {
        fmpz_clear(p_fmpz);
        fq_ctx_clear(ctx);
        //return -1;
        exit(1);
    }

    int *temp_coeffs = malloc(k * sizeof(int));
    if (!temp_coeffs) {
        free(elements);
        fmpz_clear(p_fmpz);
        fq_ctx_clear(ctx);
        //return -1;
        exit(1);
    }

    // initialize all elements
    for (int i = 0; i < field_size; i++)
    {
        // initialize the ff element as a polynomial
        // need to create an fmpz poly to init a ff elem
        fmpz_poly_init(elements[i].poly);
        fmpz_poly_fit_length(elements[i].poly, k);
        fq_init(elements[i].ff, ctx);

        // convert index to polynomial coefficients
        int temp = i;
        for (int j = 0; j < k; j++) {
            temp_coeffs[j] = temp % p;
            temp /= p;
        }

        // set coefficients
        for (int j = 0; j < k; j++) {
            fmpz_poly_set_coeff_si(elements[i].poly, j, temp_coeffs[j]);
        }
        fq_set_fmpz_poly(elements[i].ff, elements[i].poly, ctx);
    }

    field_element_t result;
    fmpz_poly_init(result.poly);
    fq_init(result.ff, ctx);

    for (int i = 0; i < field_size; i++) {
        for (int j = i; j < field_size; j++) {
            // addition
            fq_add(result.ff, elements[i].ff, elements[j].ff, ctx);
            fq_get_fmpz_poly(result.poly, result.ff, ctx);

            for (int l = 0; l < k; l++) {
                temp_coeffs[l] = fmpz_poly_get_coeff_si(result.poly, l);
            }
            int add_result = convertPolynomialToInteger_opt(p, k, temp_coeffs);

            // use helper function for safe array access
            *get_element(addition_field, i, j, field_size) = add_result;
            *get_element(addition_field, j, i, field_size) = add_result;

            // multiplication
            fq_mul(result.ff, elements[i].ff, elements[j].ff, ctx);
            fq_get_fmpz_poly(result.poly, result.ff, ctx);

            for (int l = 0; l < k; l++) {
                temp_coeffs[l] = fmpz_poly_get_coeff_si(result.poly, l);
            }
            int mult_result = convertPolynomialToInteger_opt(p, k, temp_coeffs);

            // use helper function for safe array access
            *get_element(multiplication_field, i, j, field_size) = mult_result;
            *get_element(multiplication_field, j, i, field_size) = mult_result;
        }
    }
    // free memory
    fmpz_clear(p_fmpz);
    fq_clear(result.ff, ctx);
    fmpz_poly_clear(result.poly);

    for (int i = 0; i < field_size; i++) {
        fq_clear(elements[i].ff, ctx);
        fmpz_poly_clear(elements[i].poly);
    }

    free(elements);
    free(temp_coeffs);
    fq_ctx_clear(ctx);
    //return 0;
}
// horner's method
// https://en.wikipedia.org/wiki/Horner%27s_method#Polynomial_evaluation_and_long_division
int hornerPolynomialEvalOverFq(int polyLength, int polynomialCoefficients[polyLength], int x, int q, int addition_field[q][q], int multiplication_field[q][q])
{
    int polynomialSolution = polynomialCoefficients[0];
    for (int j = 1; j < polyLength; j++)
    {
        polynomialSolution = addition_field[polynomialCoefficients[j]][multiplication_field[polynomialSolution][x]];
    }
    return polynomialSolution;
}

// (could do this inside of populateFiniteField to avoid these loops)
void populateAdditiveInverses(int p, int k, int* addition_field, int* additive_inverses)
{
    int field_size = compute_field_size(p, k);
    for (int i = 0; i < field_size; i++)
    {
        for (int j = 0; j < field_size; j++)
        {
            if (*get_element(addition_field, i, j, field_size) == 0)
            {
                additive_inverses[i] = j;
            }
        }
    }
}

void populateMultiplicativeInverses(int p, int k, int* multiplication_field, int* multiplicative_inverses)
{
    int field_size = compute_field_size(p, k);
    multiplicative_inverses[0] = -1;
    for (int i = 1; i < field_size; i++)
    {
        for (int j = 1; j < field_size; j++)
        {
            if (*get_element(multiplication_field, i, j, field_size) == 1)
            {
                multiplicative_inverses[i] = j;
            }
        }
    }
}

// older, less optimized iteration of this code:

/*
// k tuple lex unrank (aka: convert base "base" number with k digits from "number" to base 10)
int convertPolynomialToInteger(int base, int k, int* number)
{
    int ret = 0;
    int q = 1;

    for (int i = k - 1; i >= 0; i--)
    {
        ret += number[i] * q;
        q *= base;
    }

    return ret;
}

void populateFiniteField(int p, int k, int addition_field[((int)pow(p,k))][((int)pow(p,k))], int multiplication_field[((int)pow(p,k))][((int)pow(p,k))])
{
    fmpz_poly_t add_result_poly, mult_result_poly, a_poly, b_poly;
    int a_coeffs[k];
    int b_coeffs[k];
    int add_result_coeffs[k];
    int mult_result_coeffs[k];
    fmpz_poly_init(a_poly);
    fmpz_poly_init(b_poly);
    fmpz_poly_init(add_result_poly);
    fmpz_poly_init(mult_result_poly);
    fmpz_poly_fit_length(a_poly, k);
    fmpz_poly_fit_length(b_poly, k);
    fmpz_poly_fit_length(add_result_poly, k);
    fmpz_poly_fit_length(mult_result_poly, k);

    // init the context for the finite field
    long k_slong = k;
    fq_ctx_t ctx;
    fmpz_t p_fmpz;
    fmpz_init(p_fmpz);
    fmpz_set_ui(p_fmpz, p);
    fq_ctx_init(ctx, p_fmpz, k_slong, "c");

    // init finite field elements
    fq_t a_ff, b_ff, add_result_ff, mult_result_ff;
    fq_init(a_ff, ctx);
    fq_init(b_ff, ctx);
    fq_init(add_result_ff, ctx);
    fq_init(mult_result_ff, ctx);


    setToAllZeros(k, a_coeffs);
    int x = 0,y;
    do
    {
        setToAllZeros(k, b_coeffs);
        for (int i = 0; i < k; i++)
        {
            fmpz_poly_set_coeff_si(a_poly, i, a_coeffs[i]);
        }
        fq_set_fmpz_poly(a_ff, a_poly, ctx);
        y = 0;
        do
        {
            for (int i = 0; i < k; i++)
            {
                fmpz_poly_set_coeff_si(b_poly, i, b_coeffs[i]);
            }
            fq_set_fmpz_poly(b_ff, b_poly, ctx);
            fq_add(add_result_ff, a_ff, b_ff, ctx);
            fq_mul(mult_result_ff, a_ff, b_ff, ctx);
            fq_get_fmpz_poly(add_result_poly, add_result_ff, ctx);
            fq_get_fmpz_poly(mult_result_poly, mult_result_ff, ctx);
            for (int asdf = 0; asdf < k; asdf++)
            {
                add_result_coeffs[asdf] = fmpz_poly_get_coeff_si(add_result_poly, asdf);
                mult_result_coeffs[asdf] = fmpz_poly_get_coeff_si(mult_result_poly, asdf);
            }
            int add_result = convertPolynomialToInteger(p, k, add_result_coeffs);
            int mult_result = convertPolynomialToInteger(p, k, mult_result_coeffs);
            multiplication_field[x][y] = mult_result;
            multiplication_field[y][x] = mult_result;

            addition_field[x][y] = add_result;
            addition_field[y][x] = add_result;
            y++;
        } while (nextLexicographicTuple(p,k,b_coeffs));
        x++;
    } while (nextLexicographicTuple(p,k,a_coeffs));

    // clean up flint stuff
    fmpz_clear(p_fmpz);
    fmpz_poly_clear(a_poly);
    fmpz_poly_clear(b_poly);
    fmpz_poly_clear(mult_result_poly);
    fmpz_poly_clear(add_result_poly);
    fq_clear(a_ff, ctx);
    fq_clear(b_ff, ctx);
    fq_clear(add_result_ff, ctx);
    fq_clear(mult_result_ff, ctx);
    fq_ctx_clear(ctx);
}
*/