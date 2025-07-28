#include "ID_Construction.h"

#include <stdbool.h>
#include <string.h>

#include "../CFF_Tables.h"
#include "../CFF.h"


void idShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "ID");
}

void idLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "ID($%hd$)", consParams[0]);
}

void idConstructCFF(int d, int t)
{
    global_tables_array[d-1]->array[t].cff = idMatrix(t, d);
}

CFF_Construction_And_Name_Functions idConstructionFunctions = {
    .shortSrcFormatter = idShortSrcFormatter,
    .longSrcFormatter = idLongSrcFormatter,
    .constructionFunction = idConstructCFF
};

CoverFreeFamily* idMatrix(int n, int d)
{
    CoverFreeFamily *result = allocateCFF(d, n, n);
    for (int i = 0; i < n; i++)
    {
        setCFFMatrixValue(result, i, i, 1);
    }
    return result;
}