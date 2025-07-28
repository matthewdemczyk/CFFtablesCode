#include "Ext_By_One_Construction.h"

#include <string.h>
#include <stdbool.h>

#include "../CFF_Tables.h"
#include "../CFF.h"
#include "Additive_Construction.h"


void extShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "Extension by one");
}

void extLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "Extension by one of $%hd$", consParams[0]);
}

void extConstructCFF(int d, int t)
{
    CoverFreeFamily *toExtend = getByT(d, global_tables_array[d-1]->array[t].consParams[0]);
    global_tables_array[d-1]->array[t].cff = addConstruction(toExtend, getByT(d, 1));
}

CFF_Construction_And_Name_Functions extConstructionFunctions = {
    .shortSrcFormatter = extShortSrcFormatter,
    .longSrcFormatter = extLongSrcFormatter,
    .constructionFunction = extConstructCFF
};

void extendByOne(CFF_Table *table, int cff_d)
{
    for (int t = cff_d+1; t < table->numCFFs; t++)
    {
        updateTable(table, t + 1, table->array[t].n + 1, &extConstructionFunctions, t, 0, 0, 0, 0);
    }
}