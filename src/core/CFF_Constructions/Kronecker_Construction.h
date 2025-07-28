#ifndef KRONECKER_CONSTRUCTION_HEADER
#define KRONECKER_CONSTRUCTION_HEADER

#include "../CFF.h"
#include "../CFF_Tables.h"

void kroneckerShortSrcFormatter(char *strBuffer);
void kroneckerLongSrcFormatter(short *consParams, char *strBuffer);
void kroneckerConstructCFF(int d, int t);

CoverFreeFamily* KroneckerProductConstruction(CoverFreeFamily *cff_left, CoverFreeFamily *cff_right);

// must be declared extern because it's used in
// src/core/CFF_Constructions/Optimized_Kronecker_Construction.c
extern CFF_Construction_And_Name_Functions kroneckerConstructionFunctions;

#endif