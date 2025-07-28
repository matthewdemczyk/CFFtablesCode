#ifndef ADDITIVE_CFF_CONSTRUCTION_HEADER
#define ADDITIVE_CFF_CONSTRUCTION_HEADER

#include "../CFF.h"
#include "../CFF_Tables.h"

void addShortSrcFormatter(char *strBuffer);
void addLongSrcFormatter(short *consParams, char *strBuffer);
void addConstructCFF(int d, int t);

CoverFreeFamily* addConstruction(CoverFreeFamily *leftCFF, CoverFreeFamily *rightCFF);

// must be declared extern because it's used in
// src/core/CFF_Constructions/Optimized_Kronecker_Construction.c
extern CFF_Construction_And_Name_Functions addConstructionFunctions;

#endif