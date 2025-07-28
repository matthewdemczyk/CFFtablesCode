#include "../core/CFF.h"
#include "../core/CFF_Constructions/ID_Construction.h"

#include <stdlib.h>

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("give 1 parameter: t\n");
        return 1;
    }

    int t = atoi(argv[1]);

    CoverFreeFamily* cff = idMatrix(t, t - 1);
    printCoverFreeFamily(cff);
    freeCFF(cff);
    return 0;
}