#include "../core/CFF.h"
#include "../core/CFF_Constructions/Sperner_Construction.h"

#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc != 2) {
        printf("give 1 parameter: t\n");
        return 1;
    }

    int t = atoi(argv[1]);

    CoverFreeFamily *cff = SpernerConstruction(t);
    printCoverFreeFamily(cff);
    puts("Starting verify...");
    if (verifyCoverFreeFamily(cff))
    {
        puts("cff valid");
    } else
    {
        puts("CFF INVALID");
    }
    freeCFF(cff);
    return 0;
}