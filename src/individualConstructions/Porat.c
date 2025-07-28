#include "../core/CFF.h"
#include "../core/CFF_Constructions/Porat_Construction.h"
#include "../core/CFF_Constructions/Finite_Fields_Wrapper.h"
#include "../core/CombinatorialAlgorithms.h"

#include <math.h>
#include <stdlib.h>

void printMatrix(int* m, int a, int b)
{
    printf("Matrix [rows=%d x cols=%d]:\n", a, b);
    int (*g)[b] = (int (*)[b])m;
    for (int i = 0; i < a; i++)
    {
        for (int j = 0; j < b; j++)
        {
            printf("%d ", g[i][j]);
        }
        printf("\n");
    }
}

// ensures that every codeword has at least gs.minDistance number of non-zero elements
bool verifyCodeWeight(GeneratorMatrixStruct gs)
{
    // 2d array pointer to gs.code
    int (*code)[gs.numCodewords] = (int (*)[gs.numCodewords]) gs.code;
    int codewordWeight;
    for (int i = 1; i < gs.numCodewords; i++)
    {
        codewordWeight = 0;
        for (int j = 0; j < gs.m; j++)
        {
            if (code[j][i] != 0)
            {
                codewordWeight++;
            }
        }
        if (codewordWeight < gs.minDistance)
        {
            return false;
        }

    }
    return true;
}

bool verifyMinDistance(GeneratorMatrixStruct gs)
{
    int (*code)[gs.numCodewords] = (int (*)[gs.numCodewords]) gs.code;
    int currentRowsToCompareDistance;
    int experimentalMinDistance = gs.m;
    int rowsToCompare[2] = {0,1};
    do
    {
        currentRowsToCompareDistance = 0;
        for (int i = 0; i < gs.m; i++)
        {
            if (code[i][rowsToCompare[0]] != code[i][rowsToCompare[1]])
            {
                currentRowsToCompareDistance++;
            }
        }
        if (currentRowsToCompareDistance < experimentalMinDistance)
        {
            experimentalMinDistance = currentRowsToCompareDistance;
        }
        if (currentRowsToCompareDistance < gs.minDistance)
        {
            printf("\nCodewords %d, %d have %d distance, but should have >=%d distance. Code is invalid.\n", rowsToCompare[0], rowsToCompare[1], currentRowsToCompareDistance, gs.minDistance);
            return false;
        }
    } while (kSubsetLexSuccessor(gs.numCodewords, 2, rowsToCompare));
    printf("\nCode is valid and has actual min distance = %d\n", experimentalMinDistance);
    return true;
}


int main(int argc, char *argv[])
{
    // paper example
    //int q = 3;
    //int k = 2;
    //int r = 3;
    //int m = 3;

    // 2-cff example from table. m=0 means m will be calculated
    //int q = 7;
    //int k = 2;
    //int r = 3;
    //int m = 0;

    int p;
    int a;
    int k;
    int r;
    int m;
    //int q;

    if (argc == 5)
    {
        p = atoi(argv[1]);
        a = atoi(argv[2]);
        k = atoi(argv[3]);
        r = atoi(argv[4]);
        m = 0;
        int q = (int)pow(p,a);
        if (q < 2 * r || q >= 4 * r)
        {
            printf("Parameters do not match thm requirments. q not in [2r, 4r)\n");
            printf("exiting early\n");
            //exit(1);
        }
    } else if (argc == 6)
    {
        p = atoi(argv[1]);
        a = atoi(argv[2]);
        k = atoi(argv[3]);
        r = atoi(argv[4]);
        m = atoi(argv[5]);
    } else
    {
        // use defaults if no sizes passed
        p = 3;
        a = 1;
        k = 2;
        r = 3;
        m = 3;
        printf("argc != 4,5. Using default test values for p=%d,a=%d,k=%d,r=%d,m=%d\n", p,a,k,r,m);
    }
    //q = compute_field_size(p,a);
    GeneratorMatrixStruct genMatrixStruct = constructGeneratorMatrixBestBrute(p, a, k, r, m);
    printf("Finished porat construction.\n");
    printf("\nGenerator Matrix:\n");
    printMatrix(genMatrixStruct.generatorMatrix, genMatrixStruct.m, genMatrixStruct.k);
    printf("\nCode:\n");
    printMatrix(genMatrixStruct.code, genMatrixStruct.m, genMatrixStruct.numCodewords);
    if (verifyMinDistance(genMatrixStruct))
    {
        //printf("Code has valid weight\n");
    } else
    {
        printf("Code is NOT valid\n");
    }
    CoverFreeFamily* cff = genmatrixToCFF(genMatrixStruct);
    printCoverFreeFamily(cff);
    printf("Starting CFF verify...\n");
    if (verifyCoverFreeFamily(cff))
    {
        printf("CFF is valid\n");
    } else
    {
        printf("CFF is NOT valid\n");
    }
    freeCFF(cff);
    freeGeneratorMatrixStruct(genMatrixStruct);

    genMatrixStruct = constructGeneratorMatrixBestOptimized(p, a, k, r, m);
    printf("Finished porat construction.\n");
    printf("\nGenerator Matrix:\n");
    printMatrix(genMatrixStruct.generatorMatrix, genMatrixStruct.m, genMatrixStruct.k);
    printf("\nCode:\n");
    printMatrix(genMatrixStruct.code, genMatrixStruct.m, genMatrixStruct.numCodewords);
    if (verifyMinDistance(genMatrixStruct))
    {
        //printf("Code has valid weight\n");
    } else
    {
        printf("Code is NOT valid\n");
    }
    cff = genmatrixToCFF(genMatrixStruct);
    printCoverFreeFamily(cff);
    printf("Starting CFF verify...\n");
    if (verifyCoverFreeFamily(cff))
    {
        printf("CFF is valid\n");
    } else
    {
        printf("CFF is NOT valid\n");
    }
    freeCFF(cff);
    freeGeneratorMatrixStruct(genMatrixStruct);
}