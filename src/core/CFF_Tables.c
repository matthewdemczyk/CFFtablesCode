#include "CFF_Tables.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CFF_Constructions/ID_Construction.h"
#include "CFF_Constructions/STS_Construction.h"
#include "CFF_Constructions/Reed_Solomon_Construction.h"
#include "CFF_Constructions/Short_Reed_Solomon_Construction.h"
#include "CFF_Constructions/Sperner_Construction.h"
#include "CFF_Constructions/Doubling_Construction.h"
#include "CFF_Constructions/Optimized_Kronecker_Construction.h"
#include "CFF_Constructions/Kronecker_Construction.h"
#include "CFF_Constructions/Additive_Construction.h"
#include "CFF_Constructions/Ext_By_One_Construction.h"
#include "CFF_Constructions/Porat_Construction.h"
#include "CFF_Constructions/Fixed_CFFs.h"
#include "CombinatorialAlgorithms.h"
#include "GeneralAlgorithms.h"
#include "CFF.h"


unsigned t_max;
unsigned d_max;
CFF_n_type n_max;
CFF_Table **global_tables_array; //an array of pointers of size d_max

CoverFreeFamily* getByT(int d, int t)
{
    if (global_tables_array[d-1]->array[t].cff == NULL)
    {
        global_tables_array[d-1]->array[t].functions->constructionFunction(d, t);
    }
    return global_tables_array[d-1]->array[t].cff;
}

inline __attribute__((always_inline)) void updateTable(
//void updateTable(
    CFF_Table *table,
    int t,
    CFF_n_type n,
    CFF_Construction_And_Name_Functions *cffFunctions,
    int consParam0,
    int consParam1,
    int consParam2,
    int consParam3,
    int consParam4
)
{
    // if t fits, and the n the n is better: (it should go in the table)
    if ((t < table->numCFFs) && (n > table->array[t].n) && table->array[t].n != n_max)
    {
        // edge case for when the cff is bigger than n_max
        if ((n >= n_max))
        {
            table->array[t].n = n_max;
            table->numCFFs = t+1;
        }
        else
        {
            table->array[t].n = n;
        }
        table->array[t].consParams[0] = (short) consParam0;
        table->array[t].consParams[1] = (short) consParam1;
        table->array[t].consParams[2] = (short) consParam2;
        table->array[t].consParams[3] = (short) consParam3;
        table->array[t].consParams[4] = (short) consParam4;
        table->array[t].functions = cffFunctions;
        table->hasBeenChanged = true;
    }
}

int binarySearchTable(CFF_Table *table, CFF_n_type n) {
    int low = 0;
    int high = table->numCFFs - 1;

    if (table->array[high].n < n) {
        return -1; // n is not present in the range
    }

    while (low <= high) {
        // continue until low == high
        int mid = low + (high - low) / 2; // avoid overflow

        if (table->array[mid].n == n)
        {
            return mid;
        }

        if (table->array[mid].n > n) {
            high = mid - 1; // move high down
        } else {
            low = mid + 1; // move low up
        }
    }

    return low;
}

// helper used in makeTables()
CFF_Table* initializeTable(int numCFFs, int cff_d)
{
    CFF_Table *table = malloc(sizeof(CFF_Table));
    if (table == NULL)
    {
        printf("malloc fail'd table in initializeTable\n");
        exit(1);
    }
    table->hasBeenChanged = false;
    table->numCFFs = numCFFs;
    table->d = cff_d;
    // allocate space for the t=0 cff because this array should be indexed starting at 1
    table->array = malloc(sizeof(CFF_Table_Row)*(table->numCFFs+1));
    if (table->array == NULL)
    {
        printf("malloc fail'd table->array in initializeTable\n");
        exit(1);
    }
    // fill table with ID matices
    for (int t = 0; t < numCFFs + 1; t++)
    {
        table->array[t].n = t;
        table->array[t].consParams[0] = t;
        table->array[t].consParams[1] = 0;
        table->array[t].consParams[2] = 0;
        table->array[t].consParams[3] = 0;
        table->array[t].consParams[4] = 0;
        table->array[t].functions = &idConstructionFunctions;
        table->array[t].cff = NULL;
    }
    return table;
}

// helper used in makeTables()
CFF_Table* makeSpernerTable()
{
    CFF_Table *table = initializeTable(67, 1);
    for (int t = 4; t <= table->numCFFs; t++)
    {
        table->array[t].n = choose(t, t/2);
        table->array[t].consParams[0] = t;
        table->array[t].functions = &spernerConstructionFunctions;
    }
    return table;
}

void makeTables(unsigned d_maximum, unsigned t_maximum, CFF_n_type n_maximum, bool printLoops, bool useBinConstWeightCodes)
{
    // set globals for table bounds
    d_max = d_maximum;
    t_max = t_maximum;
    n_max = n_maximum;

    // save memory:
    if (t_max > n_max)
    {
        t_max = n_max;
    }

    // create an array of booleans to determine if numbers are prime/primepowers
    bool *prime_power_array = malloc(sizeof(bool)*t_max + 1);
    bool *prime_array = malloc(sizeof(bool)*t_max + 1);
    prime_power_sieve(t_max, prime_array, prime_power_array);

    // each table is stored in this array of pointers to structs, where each struct is one table
    CFF_Table **table_array = malloc(sizeof(CFF_Table*) * d_max);
    //CFF_Table* table_array[d_max];
    global_tables_array = table_array;

    // best 1-CFFs are sperner systems
    table_array[0] = makeSpernerTable();
    if (printLoops) printf("Finished d=1 table\n");

    // 2-CFFs have more constructions, so handle it seperately
    int c = 0;
    if (d_max > 1)
    {
        table_array[1] = initializeTable(t_max+1, 2);
        if (useBinConstWeightCodes)
        {
            addSurveyCFFs(table_array[1]);
        }
        addSTS(table_array[1]);
        addReedSolomonCodes(table_array[1], 2, prime_array);
        addPoratCodes(table_array[1], 2, prime_power_array);
        table_array[1]->hasBeenChanged = true;
        while(table_array[1]->hasBeenChanged)
        {
            c++;
            table_array[1]->hasBeenChanged = false;
            doublingConstructionFiller(table_array[1], table_array[0]);
            extendByOne(table_array[1], 2);
            applyPairConstructions(table_array[1], table_array[0], 2);
        }
        if (printLoops) printf("looped %d times for d=2\n", c);
    }

    //tables for d=3 ... d_max
    for (unsigned cff_d = 3; cff_d < d_max+1; cff_d++)
    {
        table_array[cff_d-1] = initializeTable(t_max+1, cff_d);
        addReedSolomonCodes(table_array[cff_d-1], cff_d, prime_array);
        addPoratCodes(table_array[cff_d-1], cff_d, prime_power_array);
        table_array[cff_d-1]->hasBeenChanged = true;
        c = 0;
        while(table_array[cff_d-1]->hasBeenChanged)
        {
            c++;
            table_array[cff_d-1]->hasBeenChanged = false;
            extendByOne(table_array[cff_d-1], cff_d);
            applyPairConstructions(table_array[cff_d-1], table_array[cff_d-2], cff_d);
        }
        if (printLoops) printf("looped %d times for d=%d\n", c, cff_d);
    }
    free(prime_array);
    free(prime_power_array);
}

// helper used in freeGlobalTableArray()
void freeTable(CFF_Table *table)
{
    for (int i = 0; i < table->numCFFs; i++)
    {
        freeCFF(table->array[i].cff);
    }
    free(table->array);
    free(table);
}

void freeGlobalTableArray()
{
    for (unsigned d = 0; d < d_max; d++)
    {
        freeTable(global_tables_array[d]);
    }
    free(global_tables_array);
}

void writeTableCSVs()
{
    for (unsigned i = 0; i < d_max; i++)
    {
        char filename[100];
        char shortSrc[100];
        char longSrc[100];
        FILE *fptr = NULL;
        sprintf(filename, "./tables/d_%d.csv", global_tables_array[i]->d);
        fptr = fopen(filename, "w");
        fprintf(fptr, "t,n,short source,long source");
        for (int t = 0; t < global_tables_array[i]->numCFFs; t++)
        {
            global_tables_array[i]->array[t].functions->shortSrcFormatter(shortSrc);
            global_tables_array[i]->array[t].functions->longSrcFormatter(global_tables_array[i]->array[t].consParams, longSrc);
            fprintf(
                fptr,
                "\n%d, %llu, %s, %s",
                t,
                global_tables_array[i]->array[t].n,
                shortSrc,
                longSrc
            );
        }
        fclose(fptr);
    }
}
