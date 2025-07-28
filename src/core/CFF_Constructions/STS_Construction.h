#ifndef STS_CONSTRUCTION_HEADER
#define STS_CONSTRUCTION_HEADER
#include "../CFF_Tables.h"
#include "../CFF.h"

void stsShortSrcFormatter(char *strBuffer);
void stsLongSrcFormatter(short *consParams, char *strBuffer);
void stsConstructCFF(int d, int t);

void addSTS(CFF_Table *table);

CoverFreeFamily* fromSTS(int v);

#endif