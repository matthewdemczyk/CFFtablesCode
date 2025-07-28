#include "Reed_Solomon_Construction.h"

#include <string.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../CFF.h"
#include "../CFF_Tables.h"
#include "../CombinatorialAlgorithms.h"
#include "Finite_Fields_Wrapper.h"
#include "Short_Reed_Solomon_Construction.h"


void rsShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "Reed-Solomon");
}

void rsLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "RS($%hd^{%hd}$,$%hd$,$%hd$)", consParams[0], consParams[1], consParams[2], consParams[3]);
}

void rsConstructCFF(int d, int t)
{
    global_tables_array[d-1]->array[t].cff = reedSolomonConstruction(
        global_tables_array[d-1]->array[t].consParams[0],
        global_tables_array[d-1]->array[t].consParams[1],
        global_tables_array[d-1]->array[t].consParams[2],
        global_tables_array[d-1]->array[t].consParams[3]
    );
}

CFF_Construction_And_Name_Functions rsConstructionFunctions = {
    .shortSrcFormatter = rsShortSrcFormatter,
    .longSrcFormatter = rsLongSrcFormatter,
    .constructionFunction = rsConstructCFF
};

void addReedSolomonCodes(CFF_Table *table, int cff_d, bool prime_array[])
{
    CFF_n_type cff_n;
    long long q;
    short m, short_m, short_k, exp;
    int cff_t;
    for (short p = 2; (unsigned) p <= ceil(sqrt(t_max)); p++) //loop over q
    {
        if (prime_array[p])
        {
            q = p;
            exp = 1;
            while (q <= t_max)
            {
                for (int k = 2; k < q + 1; k++)
                {
                    m = (cff_d * (k - 1)) + 1;
                    if (m > q + 1)
                    {
                        break;
                    }
                    cff_t = m * q;
                    if (cff_t < 0) break; //this means overflow happened in above line
                    cff_n = pow(q, k);
                    updateTable(table, cff_t, cff_n, &rsConstructionFunctions, p, exp, k, m, 0);
                    for (int s = 1; s + 1 <= q && s < m && s < k; s++)
                    {
                        short_m = m - s;
                        short_k = k - s;
                        cff_t = short_m * q;
                        cff_n = pow(q, short_k);
                        updateTable(table, cff_t, cff_n, &shortRsConstructionFunctions, p, exp, k, m, s);
                    }
                }
                q *= p;
                exp++;
            }
        }
    }
}

CoverFreeFamily* reedSolomonConstruction(int p, int exp, int t, int m)
{
    // populate finite field add and mult tables
    int q = pow(p, exp);
    int multiplication_field[q][q];
    int addition_field[q][q];
    populateFiniteField(p, exp, (int*) addition_field, (int*) multiplication_field);

    // allocate cff memory and fill with zeros
    CoverFreeFamily *cff = allocateCFF(
        (m - 1) / (m - (m - t + 1)), // = d
        q * m,                       // = t
        ((int) pow(q, t))            // = n
    );

    // loop over all polynomials/codewords
    int polynomialCoefficients[t];
    setToAllZeros(t, polynomialCoefficients);
    int cn = 0; //codeword number
    int ln; //letter number
    do
    {
        ln = 0;
        //set * letter
        setCFFMatrixValue(cff, polynomialCoefficients[0], cn, 1);
        ln++;
        //evaluate polynomial from polynomialCoefficients with x=i over Fq using horner's method for poly eval
        for (int i = 0; i < m - 1; i++)
        {
            int polynomialSolution = hornerPolynomialEvalOverFq(t, polynomialCoefficients, i, q, addition_field, multiplication_field);
            setCFFMatrixValue(cff, (ln * q) + polynomialSolution, cn, 1);
            ln++;
        }
        cn++;
    } while (nextLexicographicTuple(q, t, polynomialCoefficients));

   return cff;
}