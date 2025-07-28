#ifndef SHORT_RS_CONSTRUCTION_HEADER
#define SHORT_RS_CONSTRUCTION_HEADER

#include "../CFF.h"
#include "../CFF_Tables.h"

void shortRsShortSrcFormatter(char *strBuffer);
void shortRsLongSrcFormatter(short *consParams, char *strBuffer);
void shortRsConstructCFF(int d, int t);

CoverFreeFamily* shortRSconstruction(int p, int exp, int k, int m, int s);

// must be declared extern because it's used in
// src/core/CFF_Constructions/Reed_Solomon_Construction.c
extern CFF_Construction_And_Name_Functions shortRsConstructionFunctions;

#endif