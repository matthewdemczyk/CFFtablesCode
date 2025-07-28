#include "STS_Construction.h"

#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

#include "../CFF_Tables.h"
#include "../CFF.h"


CFF_Construction_And_Name_Functions stsConstructionFunctions = {
    .shortSrcFormatter = stsShortSrcFormatter,
    .longSrcFormatter = stsLongSrcFormatter,
    .constructionFunction = stsConstructCFF
};

void stsShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "STS");
}

void stsLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "STS($%hd$)", consParams[0]);
}

void stsConstructCFF(int d, int t)
{
    (void) d; // suppress compiler warnings
    global_tables_array[1]->array[t].cff = fromSTS(t);
}

void addSTS(CFF_Table *table)
{
    CFF_n_type n;
    for (unsigned v = 2; v <= t_max; v++)
    {
        if (v % 6 == 3 || v % 6 == 1)
        {
            n = ((v * (v - 1)) / 6);
            updateTable(table, v, n, &stsConstructionFunctions, v, 0, 0, 0,0);
        }
    }
}

typedef struct SteinerTripleSystem
{
    int order;
    int numBlocks;
    int blocks[][3];
} SteinerTripleSystem;

int symmetricIdempotentQuasiGroupFunction
(
    int n, // order of the quasi group
    int x, // left operand
    int y  // right operand
)
{
    return (((n + 1) / 2) * (x + y)) % n;
}

int halfIdempotentQuasiGroupFunction
(
    int n, // order of the quasi group
    int x, // left operand
    int y  // right operand
)
{
    int t = (x + y) % n;
    if (t % 2 == 0)
    {
        return t / 2;
    } else
    {
        return (t + n - 1) / 2;
    }
}

SteinerTripleSystem* boseConstruction
(
    int v,                     // order of the STS (must be 3 mod 6)
    int (*ball)(int, int, int) // quasi group function
)
{
    int b = (v * (v - 1)) / 6; // b is the number of blocks of the STS
    SteinerTripleSystem *sts = malloc(sizeof(SteinerTripleSystem) + sizeof(int[b][3]));
    sts->numBlocks = b;
    sts->order = v;

    int n = (v - 3) / 6;
    int Q = (2 * n) + 1; // Q is the order of the quasi group

    // populate quasi group
    int quasiGroup[Q][Q];
    for (int x = 0; x < Q; x++)
    {
        for (int y = 0; y < Q; y++)
        {
            quasiGroup[x][y] = ball(Q, x, y);
        }
    }

    // create the blocks
    int blockNum = 0; //index of the current block in the array in the SteinerTripleSystem struct
    for (int x = 0; x <= 2 * n; x++)
    {
        for (int y = 0; y <= 2 * n; y++)
        {
            if (x == y)
            {   // Type 1 triple case
                sts->blocks[blockNum][0] = x + (Q * 0) + 1;
                sts->blocks[blockNum][1] = x + (Q * 1) + 1;
                sts->blocks[blockNum][2] = x + (Q * 2) + 1;
                blockNum++;
            } else if (x < y)
            {   // Type 2 triple case
                for (int i = 0; i < 3; i++)
                {
                    sts->blocks[blockNum][0] = x + (Q * i) + 1;
                    sts->blocks[blockNum][1] = y + (Q * i) + 1;
                    sts->blocks[blockNum][2] = quasiGroup[x][y] + (Q * ((i + 1) % 3)) + 1;
                    blockNum++;
                }
            }
        }
    }
    return sts;
}

SteinerTripleSystem* skolemConstruction
(
    int v,                     // order of the STS (must be 1 mod 6)
    int (*ball)(int, int, int) // quasi group function
)
{
    int b = (v * (v - 1)) / 6; // b is the number of blocks of the STS
    SteinerTripleSystem *sts = malloc(sizeof(SteinerTripleSystem) + sizeof(int[b][3]));
    sts->numBlocks = b;
    sts->order = v;

    int n = (v - 1) / 6;
    int Q = 2 * n; // Q is the order of the quasi group

    // populate quasi group
    int quasiGroup[Q][Q];
    for (int x = 0; x < Q; x++)
    {
        for (int y = 0; y < Q; y++)
        {
            quasiGroup[x][y] = ball(Q, x, y);
        }
    }

    // create the blocks
    int blockNum = 0; // index of the current block in the array in the SteinerTripleSystem struct
    int inf = v;
    for (int x = 0; x <= n-1; x++)
    {
        // Type 1 triple case
        sts->blocks[blockNum][0] = x + (Q * 0) + 1;
        sts->blocks[blockNum][1] = x + (Q * 1) + 1;
        sts->blocks[blockNum][2] = x + (Q * 2) + 1;
        blockNum++;

        // Type 2 triple case
        for (int i = 0; i < 3; i++)
        {
            sts->blocks[blockNum][0] = inf;
            sts->blocks[blockNum][1] = n + x + (Q * i) + 1;
            sts->blocks[blockNum][2] = x + (Q * ((i + 1) % 3)) + 1;
            blockNum++;
        }
    }
    for (int x = 0; x <= (2 * n) - 1; x++)
    {
        for (int y = 0; y <= (2 * n) - 1; y++)
        {
            if (x < y)
            {   // Type 3 triple case
                for (int i = 0; i < 3; i++)
                {
                    sts->blocks[blockNum][0] = x + (Q * i) + 1;
                    sts->blocks[blockNum][1] = y + (Q * i) + 1;
                    sts->blocks[blockNum][2] = quasiGroup[x][y] + (Q * ((i + 1) % 3)) + 1;
                    blockNum++;
                }
            }
        }
    }
    return sts;
}

CoverFreeFamily* fromSTS(int v)
{
    SteinerTripleSystem *sts;
    if (v % 6 == 1)
    {
        sts = skolemConstruction(v, halfIdempotentQuasiGroupFunction);
    } else if (v % 6 == 3)
    {
        sts = boseConstruction(v, symmetricIdempotentQuasiGroupFunction);
    } else
    {
        return NULL;
    }
    CoverFreeFamily *cff = allocateCFF(2, sts->order, sts->numBlocks);
    for (int i = 0; i < sts->numBlocks; i++)
    {
        for (int n = 0; n < 3; n++)
        {
            setCFFMatrixValue(cff, sts->blocks[i][n] - 1, i, 1);
        }
    }
    free(sts);
    return cff;
}
