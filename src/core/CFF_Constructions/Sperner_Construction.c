#include "Sperner_Construction.h"

#include <string.h>

#include "../CFF_Tables.h"
#include "../CFF.h"
#include "../CombinatorialAlgorithms.h"


void spernerShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "Sperner");
}

void spernerLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "Sp($%hd$)", consParams[0]);
}

void spernerConstructCFF(int d, int t)
{
    global_tables_array[d-1]->array[t].cff = SpernerConstruction((int) global_tables_array[d-1]->array[t].n);
}

CFF_Construction_And_Name_Functions spernerConstructionFunctions = {
    .shortSrcFormatter = spernerShortSrcFormatter,
    .longSrcFormatter = spernerLongSrcFormatter,
    .constructionFunction = spernerConstructCFF
};

CoverFreeFamily* SpernerConstruction(int n)
{
    // first, determine t = min{s : choose(s, s/2) >= n }
    int s = 0;
    while (((long long) choose(s, s / 2)) < n) // cast choose(s,s/2) to suppress compiler warnings
    {
        s++;
    }
    int t = s;

    // allocate memory for the CFF and set its matrix to all 0s
    CoverFreeFamily *cff = allocateCFF(1, t, n);

    // create the buffer for the current t/2-subset of {1 ... t}
    // and populate it with the first (lexicographic) set
    int subset[t / 2];
    for (int i = 0; i < t / 2; i++)
    {
        subset[i] = i;
    }

    // iterate over all of the t/2-subsets and assign them to a column in the cff
    int col = 0;
    do
    {
        // iterate over the subset and assign its points to the current CFF row
        for (int i = 0; i < t / 2; i++)
        {
            if (subset[i] < t && col < n)
            {
                setCFFMatrixValue(cff, subset[i], col, 1);
            }
        }
        col++;
    } while (nextCombination(t, t / 2, subset));
    return cff;
}