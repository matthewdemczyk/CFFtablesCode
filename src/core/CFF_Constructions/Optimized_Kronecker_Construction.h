#ifndef OPT_KRONECKER_CONSTRUCTION_HEADER
#define OPT_KRONECKER_CONSTRUCTION_HEADER

#include "../CFF.h"
#include "../CFF_Tables.h"

void optKroneckerShortSrcFormatter(char *strBuffer);
void optKroneckerLongSrcFormatter(short *consParams, char *strBuffer);
void optKroneckerConstructCFF(int d, int t);

void applyPairConstructions(CFF_Table *table, CFF_Table *d_minus_one_table, int cff_d);

CoverFreeFamily* OptimizedKroneckerProduct
(
    CoverFreeFamily *kronecker_outter, // (d-1)-CFF(s,  n2)
    CoverFreeFamily *kronecker_inner,  //     d-CFF(t1, n1)
    CoverFreeFamily *bottom_cff        //     d-CFF(t2, n2)
);

#endif