#ifndef FIXED_CFFS_HEADER
#define FIXED_CFFS_HEADER

#include <stdbool.h>
#include <string.h>

#include "../CFF.h"
#include "../CFF_Tables.h"

void survShortSrcFormatter(char *strBuffer);
void survLongSrcFormatter(short *consParams, char *strBuffer);
void survConstructCFF(int d, int t);

void addSurveyCFFs(CFF_Table *table);

#endif