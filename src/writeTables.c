#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/CFF_Tables.h"


int main()
{
    makeTables(MAX_CFF_D, MAX_CFF_T, MAX_CFF_N, true, true);
    //makeTables(10, 200, 1000, true, true);
    writeTableCSVs();
    freeGlobalTableArray();
}