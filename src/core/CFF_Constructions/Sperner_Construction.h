#ifndef SPERNER_CONSTRUCTION_HEADER
#define SPERNER_CONSTRUCTION_HEADER

#include "../CFF.h"
#include "../CFF_Tables.h"

CoverFreeFamily* SpernerConstruction(int n);

void spernerShortSrcFormatter(char *strBuffer);
void spernerLongSrcFormatter(short *consParams, char *strBuffer);
void spernerConstructCFF(int d, int t);

extern CFF_Construction_And_Name_Functions spernerConstructionFunctions;

// sperner table is filled using a function in CFF_Tables.c
// (this is because d=1 is a special case)

#endif