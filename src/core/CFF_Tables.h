#ifndef CFF_TABLES_HEADER
#define CFF_TABLES_HEADER

#include <stdbool.h>
#include <float.h>

#include "CFF.h"

#define MAX_CFF_D 25
#define MAX_CFF_T 20000 // 20 thousand
#define MAX_CFF_N 10000000000000 // 10 trillion

typedef unsigned long long CFF_n_type;

extern unsigned d_max;
extern unsigned t_max;
extern CFF_n_type n_max;

typedef void (*TableSourceShortSrcFormatter)(char*);
typedef void (*TableSourceLongSrcFormatter)(short*, char*);
typedef void (*ConstructionFunction)(int, int);

typedef struct
{
    TableSourceShortSrcFormatter shortSrcFormatter;
    TableSourceLongSrcFormatter longSrcFormatter;
    ConstructionFunction constructionFunction;
} CFF_Construction_And_Name_Functions;

typedef struct
{
    CFF_n_type n;
    short consParams[5];
    CFF_Construction_And_Name_Functions *functions;
    CoverFreeFamily *cff;
} CFF_Table_Row;

typedef struct
{
    bool hasBeenChanged;
    int d;
    int numCFFs;
    CFF_Table_Row *array; //an array of length numCFFs + 1
} CFF_Table;

// an array of pointers to tables
extern CFF_Table **global_tables_array;

// main recursive function to get a CFF after making the tables
CoverFreeFamily* getByT(int d, int t);

// called in loops to check constructions to update table with newly found CFFs
void updateTable(
    CFF_Table *table,
    int t,
    CFF_n_type n,
    CFF_Construction_And_Name_Functions *cffFunctions,
    int consParam0,
    int consParam1,
    int consParam2,
    int consParam3,
    int consParam4
);

// helper to search thru table
int binarySearchTable(CFF_Table *table, CFF_n_type n);

// initializes global_tables_array
void makeTables(unsigned d_maximum, unsigned t_maximum, CFF_n_type n_maximum, bool printLoops, bool useBinConstWeightCodes);

// frees global_tables_array's memory
void freeGlobalTableArray();

// makes csv files in ./tables/ from global_tables_array
void writeTableCSVs();

#endif