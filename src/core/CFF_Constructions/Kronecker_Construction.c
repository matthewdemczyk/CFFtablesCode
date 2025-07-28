#include "Kronecker_Construction.h"

#include <string.h>

#include "../CFF.h"
#include "../CFF_Tables.h"


void kroneckerShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "Kronecker");
}

void kroneckerLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "Kr($%hd$,$%hd$)", consParams[0], consParams[1]);
}

void kroneckerConstructCFF(int d, int t)
{
    CoverFreeFamily *right = getByT(d, global_tables_array[d-1]->array[t].consParams[0]);
    CoverFreeFamily *left = getByT(d, global_tables_array[d-1]->array[t].consParams[1]);
    global_tables_array[d-1]->array[t].cff = KroneckerProductConstruction(left, right);
}

CFF_Construction_And_Name_Functions kroneckerConstructionFunctions = {
    .shortSrcFormatter = kroneckerShortSrcFormatter,
    .longSrcFormatter = kroneckerLongSrcFormatter,
    .constructionFunction = kroneckerConstructCFF
};

CoverFreeFamily* KroneckerProductConstruction(CoverFreeFamily *cff_left, CoverFreeFamily *cff_right)
{
    // Verify if the cffs' parameters are valid
    if (cff_left->d != cff_right->d)
    {
        return NULL;
    }

    // Allocate memory for the product cff (and fills matrix with 0s)
    CoverFreeFamily *product_cff = allocateCFF(
        cff_left->d, // "d" of the product cff
        cff_left->t * cff_right->t, // "t" of the product cff
        cff_left->n * cff_right->n // "n" of the product cff
    );

    // loop over cff_right
    for (int n1 = 0; n1 < cff_right->n; n1++)
    {
        for (int t1 = 0; t1 < cff_right->t; t1++)
        {
            if (getCFFMatrixValue(cff_right, t1, n1) == 1)
            {
                // loop over cff_left
                for (int n2 = 0; n2 < cff_left->n; n2++)
                {
                    for (int s = 0; s < cff_left->t; s++)
                    {
                        if (getCFFMatrixValue(cff_left, s, n2) == 1)
                        {
                            setCFFMatrixValue(
                                product_cff,
                                (t1*cff_left->t) + s,
                                (n1*cff_left->n) + n2,
                                1
                            );
                        }
                    }
                }
            }
        }
    }

    return product_cff;
}

// the function to add kronecker product construction parameters to the tables is in the file
// src/core/CFF_Constructions/Optimized_Kronecker_Construction.c
// function: void applyPairConstructions(CFF_Table *table, CFF_Table *d_minus_one_table, int cff_d);