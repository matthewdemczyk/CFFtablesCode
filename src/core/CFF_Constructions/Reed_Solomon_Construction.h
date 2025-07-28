#ifndef REED_SOLOMON_CONSTRUCTION_HEADER
#define REED_SOLOMON_CONSTRUCTION_HEADER
#include "../CFF_Tables.h"
#include "../CFF.h"

void rsShortSrcFormatter(char *strBuffer);
void rsLongSrcFormatter(short *consParams, char *strBuffer);
void rsConstructCFF(int d, int t);

void addReedSolomonCodes(CFF_Table *table, int cff_d, bool prime_array[]);

CoverFreeFamily* reedSolomonConstruction(int p, int exp, int t, int m);

#endif