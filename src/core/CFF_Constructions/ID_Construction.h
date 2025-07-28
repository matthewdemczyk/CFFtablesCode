#ifndef ID_CONSTRUCTION_HEADER
#define ID_CONSTRUCTION_HEADER
#include "../CFF_Tables.h"
#include "../CFF.h"

CoverFreeFamily* idMatrix(int n, int d);

void idShortSrcFormatter(char *strBuffer);
void idLongSrcFormatter(short *consParams, char *strBuffer);
void idConstructCFF(int d, int t);

// used in other files so must be extern
extern CFF_Construction_And_Name_Functions idConstructionFunctions;

#endif