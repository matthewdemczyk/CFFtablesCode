#ifndef EXT_BY_ONE_CONSTRUCTION_HEADER
#define EXT_BY_ONE_CONSTRUCTION_HEADER

#include <stdbool.h>

#include "../CFF_Tables.h"
#include "../CFF.h"

void extShortSrcFormatter(char *strBuffer);
void extLongSrcFormatter(short *consParams, char *strBuffer);
void extConstructCFF(int d, int t);

void extendByOne(CFF_Table *table, int cff_d);

#endif