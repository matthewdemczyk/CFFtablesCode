#include "CombinatorialAlgorithms.h"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "GeneralAlgorithms.h"


/*
    Returns nCk (n choose k) by generating Pascal's triangle.
    Returns 0 if overflow occurrs during addition in the function,
    or if the parameters are invalid (n < 1 OR n < k OR k < 0)
*/
unsigned long long choose(int n, int k)
{
    // invalid inputs
    if (n < k || n < 0 || k < 0)
    {
        return 0;
    }

    // only calculate left side
    if (n - k < k)
    {
        k = n - k;
    }

    // 1 cases
    if (n < 2 || k == 0)
    {
        return 1;
    }

    unsigned long long a[k+1]; // previous row of pascal's triangle
    unsigned long long b[k+1]; // current row of pascal's triangle

    unsigned long long *a_ptr, *b_ptr, *swap_ptr;
    a_ptr = a;
    b_ptr = b;

    // the leftmost cell of every row is one:
    a_ptr[0] = 1;
    b_ptr[0] = 1;
    // iterate over the n + 1 rows in pascal's triangle
    for (int row = 1; row < n + 1; row++)
    {
        int col; // (iniitialize outside for loop so its value can be accessed after)
        for (col = 1; col < row && col < k + 1; col++)
        {
            b_ptr[col] = a_ptr[col-1] + a_ptr[col];
            // returns 0 if overflow occurred during the addition in above line
            if (b_ptr[col] < a_ptr[col-1] || b_ptr[col] < a_ptr[col])
            {
                return 0;
            }
        }
        b_ptr[col] = 1;
        swap_ptr = a_ptr;
        a_ptr = b_ptr;
        b_ptr = swap_ptr;
    }
    return a_ptr[k];
}



/*
    Sets buffer to {0, 1, ..., k-1}
*/
void setToFirstLexOrder(int k, int buffer[k])
{
    for(int i = 0; i < k; i++)
    {
        buffer[i] = i;
    }
}


/*
    Replaces "buffer" with the next consecutive lexicographic subset of Zn (size k).
    "buffer" is assummed to be set to first lexicographic order on first call.
*/
bool kSubsetLexSuccessor(int n, int k, int buffer[k])
{
    for (int i = k-1; i > -1; i--)
    {
        if (buffer[i] != n-k+i)
        {
            buffer[i]++;
            for (int x = i+1; x < k; x++)
            {
                buffer[x] = buffer[i] + (x-i);
            }
            return true;
        }
    }
    return false;
}


bool nextCombination(int n, int k, int buffer[k])
{
    return kSubsetLexSuccessor(n,k,buffer);
}


int kSubsetLexRank(int n, int k, int buffer[k])
{
    int r = 0;
    int t[k+1];
    t[0] = 0;
    for (int a = 1; a < k+1; a++)
    {
        t[a] = buffer[a-1] + 1;
    }
    for (int i = 1; i <= k; i++)
    {
        if ((t[i-1] + 1) <= (t[i] - 1))
        {
            for (int j = t[i-1] + 1; j <= t[i] - 1; j++)
            {
                r = r + choose(n-j, k-i);
                //printf("i=%d, j=%d, r=%d\n", i, j, r);
            }
        }
    }
    return r;
}


//Page 45 combinatorial algorithms
void kSubsetLexUnrank(int n, int k, int r, int buffer[k])
{
    int x = 1;
    int algr = r;
    unsigned long long temp;
    for(int i = 1; i <= k; i++)
    {
        temp = choose(n-x, k-i);
        while((long long)temp <= algr)
        {
            algr = algr - temp;
            x++;
            temp = choose(n-x, k-i);
        }
        buffer[i-1] = x-1;
        x++;
    }
}


//Page 52 combinatorial algorithms
bool kSubsetRevDoorSuccessor(int n, int k, int t[k+2])
{
    t[0] = -1;
    t[k+1] = n+1;
    int j = 1;
    //return false if no successor:
    bool hasSuccessor = false;
    for (int i = 1; i < k; i++)
    {
        if (t[i] != i)
        {
            hasSuccessor = true;
        }
    }
    if (!hasSuccessor && t[k] == n)
    {
        return false;
    }

    //regular alg:
    while (j <= k && t[j] == j)
    {
        j++;
    }
    if (k % 2 != j % 2)
    {
        if (j == 1)
        {
            t[1]--;
        } else
        {
            t[j-1] = j;
            t[j-2] = j-1;
        }
    } else
    {
        if (t[j+1] != t[j] + 1)
        {
            t[j-1] = t[j];
            t[j]++;
        } else
        {
            t[j+1] = t[j];
            t[j] = j;
        }
    }
    return true;
}

int kSubsetRevDoorRank(int k, int t[k+2])
{
    int r = -1 * (k % 2);
    int s = 1;
    for (int i = k; i > 0; i--)
    {
        r += s * choose(t[i], i);
        s *= -1;
    }
    return r;
}

void kSubsetRevDoorUnrank(int n, int k, int r, int t[k])
{
    int x = n;
    for (int i = k; i > 0; i--)
    {
        while ((long long)choose(x, i) > r)
        {
            x--;
        }
        t[i-1] = x + 1;
        r = choose(x + 1, i) - r - 1;
    }
}


int printKsubsetRevDoor(int n, int k)
{
    int buffer[k+2];
    int unrankedBuffer[k];
    int counter = 0;
    //unsigned long long total = choose(n,k);
    //setToFirstLexOrder(k, buffer);
    for (int i = 0; i < k; i++)
    {
        buffer[i+1] = i + 1;
    }
    buffer[0] = 0;
    do
    {
        for (int i = 0; i <= k + 1; i++)
        {
             printf("%d ", buffer[i]);
        }
        //printf("\n");
        printf("[%d] ", kSubsetRevDoorRank(k, buffer));
        //printf("\n");
        kSubsetRevDoorUnrank(n, k, counter, unrankedBuffer);
        for (int i = 0; i < k; i++)
        {
             printf("%d ", unrankedBuffer[i]);
        }
        printf("\n");
        counter++;
        //if (counter == (long long) total) {break;}
    } while (kSubsetRevDoorSuccessor(n, k, buffer));
    return counter;
}


//From page 53 of combinatorial algorithms
bool permLexSuccessor(int n, int buffer[n])
{
    //Set first of buffer to 0 to ensure that the while loop terminates with 0<=i<=n-1
    //buffer[0] = -1;
    int t[n+1];
    t[0] = -1;
    for (int a = 1; a < n+1; a++)
    {
        t[a] = buffer[a-1];
    }
    int i = n - 1;
    while (t[i+1] < t[i])
    {
        i--;
    }
    //printf("i=%d ", i);
    if (i == 0)
    { //In this case buffer = [n, n-1, ... , 1] and is the last permuation lexicographically (->there is no successor)
        return false;
    }
    // J will be the position of the last element among buffer[i+1 : n] that is greater than all the rest from [j+1:n]
    int j = n;
    while (t[j] < t[i])
    {
        j--;
    }
    //printf("j=%d ", j);
    int swap = t[j];
    t[j] = t[i];
    t[i] = swap;
    int p[n];
    for (int h = i + 1; h <= n; h++)
    {
        p[h] = t[h];
    }
    for (int h = i + 1; h <= n; h++)
    {
        t[h] = p[n + i + 1 - h];
    }
    for (int a = 1; a < n+1; a++)
    {
        buffer[a-1] = t[a];
    }
    return true;
}


long factorial(int n)
{
    long result = 1;
    while (n > 0)
    {
        result = result * n;
        n--;
    }
    return result;
}


int permLexRank(int n, int buffer[n])
{
    int r = 0;
    int p[n+1];
    for (int a = 0; a < n+1; a++)
    {
        p[a] = buffer[a-1] + 1;
    }
    for (int j = 1; j <= n; j++)
    {
        r = r + ((p[j]-1)*factorial(n-j));
        for (int i = j + 1; i <=n; i++)
        {
            if (p[i] > p[j])
            {
                p[i]--;
            }
        }
    }
    return r;
}


void permLexUnrank(int n, int r, int buffer[n])
{
    buffer[n-1] = 1;
    for (int j = 1; j <= n - 1; j++)
    {
        int d = (r % factorial(j+1)) / factorial(j);
        r = r - (d * factorial(j));
        buffer[n-j-1] = d + 1;
        for (int i = n-j+1; i <= n; i++)
        {
            if (buffer[i-1] > d)
            {
                buffer[i-1]++;
            }
        }
    }
}


// This function should be called something like "next Lexicographic Tuple With Replacement" (of size k with max n)
// because it doesn't actually generate permutations
bool nextPermutation(int n, int k, int buffer[k])
{
    for (int i = k-1; i > -1; i--)
    {
        if (buffer[i] < n - 1)
        {
            buffer[i]++;
            for (int x = i+1; x < k; x++)
            {
                buffer[x] = 0;
            }
            return true;
        }
    }
    return false;
}


//rename but keep compatibility (compiler optimizes out the extra funciton call?)
bool nextLexicographicTuple(int n, int k, int buffer[k])
{
    return nextPermutation(n, k, buffer);
}


/*
    Prints all subsets of Zn (size k) in lexicographic order.
    Returns the number of subsets (which is n choose k).
*/
unsigned long long printCombs(int n, int k)
{
    // initialize "combinations" to an array of elements 0 ... k-1
    int combination[k];
    setToFirstLexOrder(k, combination);
    int ksubsetlexranktest[k];
    int r;
    unsigned long long c = 0;
    printf("All (%d choose %d) sets:\n", n, k);
    do
    {
        // print out the combination array
        for (int i = 0; i < k; i++)
        {
             printf("%d ", combination[i]);
        }
        r = kSubsetLexRank(n, k, combination);
        printf("[r=%d, c=%llu] ", r, c);
        kSubsetLexUnrank(n,k,(int)c,ksubsetlexranktest);
        for (int i = 0; i < k; i++)
        {
             printf("%d ", ksubsetlexranktest[i]);
        }
        printf("\n");
        c++;
    } while (kSubsetLexSuccessor(n, k, combination));
    return c;
}


/*
    Prints all permutations of size k from [0..n]
*/
unsigned long long printPerms(int n)
{
    int perms[n];
    //perms[0] = -1;
    for (int i = 0; i < n; i++)
    {
        perms[i] = i;
    }
    unsigned long long c = 0;
    int unRankedPerms[n];
    int r;
    do
    {
        c++;
        for (int i = 0; i < n; i++)
        {
             printf("%d ", perms[i]);
        }
        r = permLexRank(n, perms);
        printf("[%d] unranked:", r);
        permLexUnrank(n, r, unRankedPerms);
        for (int i = 0; i < n; i++)
        {
             printf("%d ", unRankedPerms[i]);
        }
        printf("\n");
    } while (permLexSuccessor(n, perms));
    return c;
}


/*
    Replaces T with the next successive subset in a minimal change ordering.
            (T is a 0-1 charactarisic vector)

    Returns true if there's a successor, false otherwise.
*/
bool grayCodeSuccessor(int n, int T[n])
{
    // Count the number of ones in the characteristic vector
    // Also keep track of where the rightmost "one" is
    int numOnes = 0;
    int positionOfRightmostOne = 0;
    for (int i = 0; i < n; i++)
    {
        if (T[i] == 1)
        {
            numOnes++;
            positionOfRightmostOne = i;
        }
    }

    if (numOnes == 1 && positionOfRightmostOne == 0)
    { // Return false if there is no successor
        return false;
    }

    if (numOnes % 2 == 0)
    { // Flip the last bit
        T[n - 1] = 1 - T[n - 1];
    } else
    { // Flip the bit left of the rightmost "one"
        T[positionOfRightmostOne - 1] = 1 - T[positionOfRightmostOne - 1];
    }

    return true;
}

/*
    Sets T to all zeros
*/
void setToAllZeros(int n, int T[n])
{
    for(int i = 0; i < n; i++)
    {
        T[i] = 0;
    }
}

bool subsetLexSuccessor(int k, int T[k])
{
    for (int i = k-1; i > -1; i--)
    {
        if (T[i] == 0)
        {
            T[i]++;
            for (int x = i+1; x < k; x++)
            {
                T[x] = 0;
            }
            return true;
        }
    }
    return false;
}

int subsetLexRank(int n, int T[n])
{
    int r = 0;
    for (int i = 0; i < n; i++)
    {
        if (T[i] == 1)
        {
            r += ipow(2, n-(i+1));
        }
    }
    return r;
}

void subsetLexUnrank(int n, int r, int T[n])
{
    for (int i = n; i > 0; i--)
    {
        if (r % 2 == 1)
        {
            T[i-1] = 1;
        } else
        {
            T[i-1] = 0;
        }
        r = r/2;
    }
}