#include "Fixed_CFFs.h"

#include <stdbool.h>

#include "../CFF_Tables.h"
#include "../CFF.h"
#include "ID_Construction.h"


void survShortSrcFormatter(char* strBuffer)
{
    strcpy(strBuffer, "Constant-weight binary code");
}

void survLongSrcFormatter(short* consParams, char* strBuffer)
{
    sprintf(strBuffer, "Survey CFF %hd", consParams[0]);
}

void survConstructCFF(int d, int t)
{
    //global_tables_array[d-1]->array[t].cff = readFromFile(d,t,global_tables_array[d-1]->array[t].n);
    global_tables_array[d-1]->array[t].cff = idMatrix(t, d);
}

CFF_Construction_And_Name_Functions survConstructionFunctions = {
    .shortSrcFormatter = survShortSrcFormatter,
    .longSrcFormatter = survLongSrcFormatter,
    .constructionFunction = survConstructCFF
};

void addSurveyCFFs(CFF_Table* table)
{
    table->array[10].n = 13;
    table->array[10].consParams[0] = 10;
    table->array[10].functions = &survConstructionFunctions;

    table->array[11].n = 17;
    table->array[11].consParams[0] = 11;
    table->array[11].functions = &survConstructionFunctions;


    table->array[12].n = 20;
    table->array[12].consParams[0] = 12;
    table->array[12].functions = &survConstructionFunctions;

    table->array[13].n = 26;
    table->array[13].consParams[0] = 13;
    table->array[13].functions = &survConstructionFunctions;

    table->array[14].n = 28;
    table->array[14].consParams[0] = 14;
    table->array[14].functions = &survConstructionFunctions;

    table->array[15].n = 42;
    table->array[15].consParams[0] = 15;
    table->array[15].functions = &survConstructionFunctions;

    table->array[16].n = 48;
    table->array[16].consParams[0] = 16;
    table->array[16].functions = &survConstructionFunctions;

    table->array[17].n = 68;
    table->array[17].consParams[0] = 17;
    table->array[17].functions = &survConstructionFunctions;

    table->array[18].n = 69;
    table->array[18].consParams[0] = 18;
    table->array[18].functions = &survConstructionFunctions;

    table->array[19].n = 76;
    table->array[19].consParams[0] = 19;
    table->array[19].functions = &survConstructionFunctions;

    table->array[20].n = 90;
    table->array[20].consParams[0] = 20;
    table->array[20].functions = &survConstructionFunctions;

    table->array[21].n = 120;
    table->array[21].consParams[0] = 21;
    table->array[21].functions = &survConstructionFunctions;

    table->array[22].n = 176;
    table->array[22].consParams[0] = 22;
    table->array[22].functions = &survConstructionFunctions;

    table->array[23].n = 253;
    table->array[23].consParams[0] = 23;
    table->array[23].functions = &survConstructionFunctions;
}
