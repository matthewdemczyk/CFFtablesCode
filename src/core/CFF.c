#include "CFF.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#include "CombinatorialAlgorithms.h"
#include "GeneralAlgorithms.h"

// if true, when running verifyCoverFreeFamily(), it will print out which rows in
// the CFF the d+1 rows of the ID matrix are found on
#define TESTMODE false

// allocates a d-CFF(t,n) filled with 0s
CoverFreeFamily* allocateCFF(int d, int t, int n)
{
    CoverFreeFamily* c = malloc(sizeof(CoverFreeFamily));
    c->d = d;
    c->t = t;
    c->n = n;
    c->stride_bits = ((n + 7) / 8) * 8;
    c->matrix = calloc(((n + 7) / 8) * t, sizeof(uint8_t));
    return c;
}

// free a CFF from memory
void freeCFF(CoverFreeFamily *cff)
{
    if (cff!=NULL) free(cff->matrix);
    free(cff);
}

// setter for the row "r" and column "c" for a value in the 0-1 CFF Matrix.
void setCFFMatrixValue(CoverFreeFamily *cff, int r, int  c, int val)
{
    int bitIndex = r * cff->stride_bits + c;
    if (val)
        cff->matrix[bitIndex / 8] |= (1 << (bitIndex % 8));
    else
        cff->matrix[bitIndex / 8] &= ~(1 << (bitIndex % 8));
}

// getter for the row "r" and column "c" for a value in the 0-1 CFF Matrix.
int getCFFMatrixValue(CoverFreeFamily *cff, int r, int c)
{
    int bitIndex = r * cff->stride_bits + c;
    return (cff->matrix[bitIndex / 8] >> (bitIndex % 8)) & 1;
}

// allocates memory for a d-CFF(t,n), sets its matrix to the passed matrix, then returns a pointer to it
// make sure to free the memory using freeCFF() when finished with the CFF.
CoverFreeFamily* initCoverFreeFamilyFromMatrix(int d, int t, int n, int matrix[t][n])
{
    CoverFreeFamily *cff = allocateCFF(d,t,n);
    for(int r  = 0; r < t; r++)
    {
        for(int c  = 0; c < n; c++)
        {
            setCFFMatrixValue(cff, r,c,matrix[r][c]);
        }
    }
    return cff;
}

CoverFreeFamily* copyCoverFreeFamily(CoverFreeFamily *inCFF)
{
    CoverFreeFamily *outCFF = allocateCFF(
        inCFF->d,
        inCFF->t,
        inCFF->n
    );
    for(int r = 0; r < inCFF->t; r++)
    {
        for(int c = 0; c < inCFF->n; c++)
        {
            setCFFMatrixValue(outCFF,r,c,
                getCFFMatrixValue(inCFF,r,c)
            );
        }
    }
    return outCFF;
}


// prints a CFF to console
void printCoverFreeFamily(CoverFreeFamily *cff)
{
    printf("%d-CFF(%d,%d):\n", cff->d, cff->t, cff->n);
    int value;
    for(int r = 0; r < cff->t; r++)
    {
        for(int c = 0; c < cff->n; c++)
        {
            value = getCFFMatrixValue(cff,r,c);
            if (value == 1)
            {
                printf("1 ");
            } else if (value == 0)
            {
                printf("- ");
            }
            else
            {   // This should never happen
                printf("? ");
            }
        }
        printf("\n");
    }
}

void writeCFFtoFile(CoverFreeFamily *cff, FILE *fptr)
{
    int value;
    for(int r = 0; r < cff->t; r++)
    {
        for(int c = 0; c < cff->n; c++)
        {
            value = getCFFMatrixValue(cff,r,c);
            if (value == 1)
            {
                fprintf(fptr, "1");
            } else if (value == 0)
            {
                fprintf(fptr, "0");
            }
            else
            {   // This should never happen
                fprintf(fptr, "?");
            }
        }
        printf("\n");
    }
}


// returns true if the CFF is a valid cover free family.
// returns false otherwise.
bool verifyCoverFreeFamily(CoverFreeFamily *cff)
{
    if (cff->d+1 > cff->n)
    { // return false if the parameters are invalid
        return false;
    }

    // cols will be an array of columns of size d+1 to test
    int k = cff->d + 1;
    int cols[k];
    for (int i = 0; i < k; i++)
    { // set cols to the smallest lexicographic ordering
        cols[i] = i;
    }

    do
    {
        if (TESTMODE)
        { //print out the current columns that are being tested
            printf("Testing cols:  ");
            for (int x = 0; x < k; x++)
            {
                printf("%d  ", cols[x]);
            }
            printf("| ID Matrix found on rows:  ");
        }

        // verify array "v". keeps track of the first i.d. row for every col in "cols"  array
        int v[k];
        for (int x = 0; x < k; x++)
        {
            v[x] = -1;
        }

        // s is sum of the current row, e will be the most recent column where a 1 was seen
        int s, e;
        // f is the number of identity rows found so far (for this particular subset of columns)
        int f = 0;
        // iterate over the rows in the cff
        for (int r = 0; r < cff->t; r++)
        {
            // iterate over the subset of columns we are currently testing
            s = 0;
            for (int c = 0; c < k; c++)
            {
                // if the current cell is 1:
                if (getCFFMatrixValue(cff, r, cols[c]) == 1)
                { // record the position of the column and increment row sum
                    e = c;
                    s++;
                }
            }
            // check if the row was the first identity row for it's column with a 1
            if (s == 1 && v[e] == -1)
            {
                if (TESTMODE) { printf("%d  ", r); }
                // record that c's identity row exists (by saving the row index)
                v[e] = r;
                f++;
                if (f == k)
                { // exit early if all identity rows are found
                    break;
                }

            }
        }

        // now return false if any of the identity rows were not there
        if (f != k)
        {
            if (TESTMODE) { printf(" Some ID rows not found! CFF is invalid\n"); }
            return false;
        }
        if (TESTMODE) { printf("\n"); }
    } while (nextCombination(cff->n, k, cols));
    return true;
}