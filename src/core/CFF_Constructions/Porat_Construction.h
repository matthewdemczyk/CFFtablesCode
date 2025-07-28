#ifndef PORAT_CONSTRUCTION_HEADER
#define PORAT_CONSTRUCTION_HEADER

#include <stdbool.h>
#include <math.h>

#include "../CFF_Tables.h"
#include "../CFF.h"

void poratShortSrcFormatter(char *strBuffer);
void poratLongSrcFormatter(short *consParams, char *strBuffer);
void poratConstructCFF(int d, int t);

void addPoratCodes(CFF_Table *table, int cff_d, bool *prime_power_array);

// will be returned, completely filled, by the porat construction
typedef struct GeneratorMatrixStruct {
    int q;
    int k;
    int r;
    int m;
    double Hq;
    double delta;
    int* generatorMatrix;
    int* code;
    int numCodewords;
    int minDistance;
} GeneratorMatrixStruct;
CoverFreeFamily* genmatrixToCFF(GeneratorMatrixStruct gs);

double poratEntropyFunction(double q, double r);

// these are several iterations of the algorithm
GeneratorMatrixStruct constructGeneratorMatrixBruteMod(int q, int k, int r, int m);
GeneratorMatrixStruct constructGeneratorMatrixBruteFF(int p, int a, int k, int r, int m);
GeneratorMatrixStruct constructGeneratorMatrixBestBrute(int p, int a, int k, int r, int m);
GeneratorMatrixStruct constructGeneratorMatrixBestOptimized(int p, int a, int k, int r, int m);
void freeGeneratorMatrixStruct(GeneratorMatrixStruct genMatrixStruct);

#endif