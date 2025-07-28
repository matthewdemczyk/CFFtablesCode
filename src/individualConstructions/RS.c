#include "../core/CFF.h"
#include "../core/CFF_Constructions/Reed_Solomon_Construction.h"

#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        puts("give 4 parameters: prime, base, k, m");
        return 1;
    }

    int prime = atoi(argv[1]);
    int e = atoi(argv[2]);
    int k = atoi(argv[3]);
    int m = atoi(argv[4]);

    CoverFreeFamily *cff = reedSolomonConstruction(prime, e, k, m);
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