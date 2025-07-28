#include "Additive_Construction.h"

#include <string.h>

#include "../CFF.h"
#include "../CFF_Tables.h"


void addShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "Additive");
}

void addLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "Add($%hd$,$%hd$)", consParams[0], consParams[1]);
}

void addConstructCFF(int d, int t)
{
    CoverFreeFamily *right = getByT(d, global_tables_array[d-1]->array[t].consParams[0]);
    CoverFreeFamily *left = getByT(d, global_tables_array[d-1]->array[t].consParams[1]);
    global_tables_array[d-1]->array[t].cff = addConstruction(left, right);
}

CFF_Construction_And_Name_Functions addConstructionFunctions = {
    .shortSrcFormatter = addShortSrcFormatter,
    .longSrcFormatter = addLongSrcFormatter,
    .constructionFunction = addConstructCFF
};

CoverFreeFamily* addConstruction(CoverFreeFamily *leftCFF, CoverFreeFamily *rightCFF)
{
    CoverFreeFamily *result = allocateCFF(
        leftCFF->d,
        leftCFF->t + rightCFF->t,
        leftCFF->n + rightCFF->n
    );
    for (int r = 0; r < leftCFF->t; r++)
    {
        for (int c = 0; c < leftCFF->n; c++)
        {
            if (getCFFMatrixValue(leftCFF, r, c) == 1)
            {
                setCFFMatrixValue(result, r, c, 1);
            }
        }
    }
    for (int r = 0; r < rightCFF->t; r++)
    {
        for (int c = 0; c < rightCFF->n; c++)
        {
            if (getCFFMatrixValue(rightCFF, r, c) == 1)
            {
                setCFFMatrixValue(result, r + leftCFF->t, c + leftCFF->n, 1);
            }
        }
    }

    return result;
}

// the function to add additive construction parameters to the tables is in the file
// src/core/CFF_Constructions/Optimized_Kronecker_Construction.c
// function: void applyPairConstructions(CFF_Table *table, CFF_Table *d_minus_one_table, int cff_d);