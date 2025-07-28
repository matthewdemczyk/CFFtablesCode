#ifndef COMBINATORIAL_ALGORITHMS_HEADER
#define COMBINATORIAL_ALGORITHMS_HEADER

#include <stdbool.h>

unsigned long long choose(int n, int k);
bool kSubsetLexSuccessor(int n, int k, int buffer[k]);
int kSubsetLexRank(int n, int k, int buffer[k]);
void kSubsetLexUnrank(int n, int k, int r, int buffer[k]);
bool kSubsetRevDoorSuccessor(int n, int k, int buffer[k+2]);
bool permLexSuccessor(int n, int buffer[n]);
long factorial(int n);
int permLexRank(int n, int buffer[n]);
void permLexUnrank(int n, int r, int buffer[n]);
bool grayCodeSuccessor(int n, int T[n]);
void setToFirstLexOrder(int k, int buffer[k]);
bool nextCombination(int n, int k, int buffer[k]);
int printKsubsetRevDoor(int n, int k);
bool nextPermutation(int n, int k, int buffer[k]);
bool nextLexicographicTuple(int n, int k, int buffer[k]);
unsigned long long printCombs(int n, int k);
unsigned long long printPerms(int n);
void setToAllZeros(int n, int T[n]);
bool subsetLexSuccessor(int k, int T[k]);
int subsetLexRank(int n, int T[n]);
void subsetLexUnrank(int n, int r, int T[n]);

#endif