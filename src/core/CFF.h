#ifndef CFF_HEADER
#define CFF_HEADER
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>

// a d-CFF(t, n):
// there are t subsets of the set n, stored as incidence matrix
typedef struct
{
    int d;
    int t;
    int n;
    int stride_bits; // helper for faster memory access
    uint8_t *matrix;
} CoverFreeFamily;

// getters/setters
void setCFFMatrixValue(CoverFreeFamily *cff, int r, int  c, int val);
int getCFFMatrixValue(CoverFreeFamily *cff, int r, int c);

// memory management
CoverFreeFamily* allocateCFF(int d, int t, int n);
void freeCFF(CoverFreeFamily *cff);

// allocates a d-CFF(t,n) and sets it from a 0-1 int matrix
CoverFreeFamily* initCoverFreeFamilyFromMatrix(int d, int t, int n, int matrix[t][n]);

// allocates and makes a copy of a cff
CoverFreeFamily* copyCoverFreeFamily(CoverFreeFamily *inCFF);

// print/write cff
void printCoverFreeFamily(CoverFreeFamily *cff);
void writeCFFtoFile(CoverFreeFamily *cff, FILE *fptr);

// verify if a CFF is valid
bool verifyCoverFreeFamily(CoverFreeFamily *cff);

#endif