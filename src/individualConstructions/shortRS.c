#include "../core/CFF.h"
#include "../core/CFF_Constructions/Short_Reed_Solomon_Construction.h"

#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc != 6) {
        printf("give 5 parameters: prime, base, k, m, s\n");
        return 1;
    }

    int prime = atoi(argv[1]);
    int base = atoi(argv[2]);
    int k = atoi(argv[3]);
    int m = atoi(argv[4]);
    int s = atoi(argv[5]);

    CoverFreeFamily *cff = shortRSconstruction(prime, base, k, m, s);
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