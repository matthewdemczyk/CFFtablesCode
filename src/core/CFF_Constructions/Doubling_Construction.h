#ifndef DOUBLEING_CONSTRUCTION_HEADER
#define DOUBLEING_CONSTRUCTION_HEADER

#include <stdbool.h>

#include "../CFF.h"
#include "../CFF_Tables.h"

CoverFreeFamily* doublingConstruction(CoverFreeFamily *cff, int s);

void doublingShortSrcFormatter(char *strBuffer);
void doublingLongSrcFormatter(short *consParams, char *strBuffer);
void doublingConstructCFF(int d, int t);

void doublingConstructionFiller(CFF_Table *table_2, CFF_Table *table_1);

#endif