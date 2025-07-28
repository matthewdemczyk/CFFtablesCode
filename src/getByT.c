#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/CFF_Tables.h"
#include "core/CFF.h"
#include "core/CFF_Constructions/Reed_Solomon_Construction.h"

int main(int argc, char *argv[])
{
    char longSrc[200];
    if (argc == 3)
    {
        int d = atoi(argv[1]);
        int t = atoi(argv[2]);
        makeTables(d, t, MAX_CFF_N, false, false);
        CoverFreeFamily* cff = getByT(d,t);
        global_tables_array[d-1]->array[t].functions->longSrcFormatter(global_tables_array[d-1]->array[t].consParams, longSrc);
        printf("%s\n",longSrc);
        printCoverFreeFamily(cff);
        printf("Starting CFF verification algorithm...\n");
        if (verifyCoverFreeFamily(cff))
        {
            printf("CFF passed verification\n");
        } else
        {
            printf("CFF **FAILED** verification\n");
        }
        freeGlobalTableArray();
    } else
    {
        printf("Give 2 parameters d,t\n");
        return 1;
    }
}