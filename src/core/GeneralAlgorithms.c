#include "GeneralAlgorithms.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>

#include "CombinatorialAlgorithms.h"

/*
    qsort (from stdlib.h) int comparator

    example:
    qsort(array, array_length, sizeof(int), int_comparator)
*/
int int_comparator(const void *a, const void *b)
{
   const int *A = a, *B = b;
   return (*A > *B) - (*A < *B);
}


/*
    Integer exponentiation function, implemented instead of using
    pow from math.h because pow uses floating point numbers and may
    cause rounding errors / loss of precision errors
*/
int ipow(int base, int exp)
{
    if (exp == 0) {return 1;}
    int result = 1;
    for (int i = 0; i < exp; i++)
    {
        result = result  * base;
    }
    return result;
}


/*
    Struct used to rankify and sort an array at the same time in rankify_and_sort()
*/
typedef struct ValueRankPair {
    int value;
    int rank;
} ValueRankPair;


/*
    qsort comparator for a ValueRankPair array used in rankify_and_sort()
*/
int compare_value_rank_pair(const void *a, const void *b)
{
    ValueRankPair* vrp_a = (ValueRankPair *)a;
    ValueRankPair* vrp_b = (ValueRankPair *)b;
    int arg1 = vrp_a->value;
    int arg2 = vrp_b->value;

    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}


/*
    Copies the contents from input_unsorted_array into output_sorted_buffer
    (in sorted order) without changing input_unsorted_array. also puts its
    ranks into output_rankified_buffer. uses qsort from stdlib.h
*/
void rankify_and_sort
    (
        int array_length,
        int input_unsorted_array[array_length],
        int output_rankified_buffer[array_length],
        int output_sorted_buffer[array_length]
    )
{
    ValueRankPair vrp_array[array_length];
    for (int aa = 0; aa < array_length; aa++)
    {
        vrp_array[aa].value = input_unsorted_array[aa];
        vrp_array[aa].rank = aa;
    }
    qsort(vrp_array, array_length, sizeof(ValueRankPair), compare_value_rank_pair);
    for (int aa = 0; aa < array_length; aa++)
    {
        output_sorted_buffer[aa] = vrp_array[aa].value;
        output_rankified_buffer[aa] = vrp_array[aa].rank;
    }

}


/*
    Returns the greatest common divisor of a and b
    using the Euclidean algorithm

    https://en.wikipedia.org/wiki/Euclidean_algorithm#Implementations
*/
int euclid_GCD(int a, int b) //C is pass by value :)
{
    while(a != b)
    {
        if (a > b)
        {
            a -= b;
        } else
        {
            b -= a;
        }
    }
    return a;
}

/*
    https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
*/
void prime_power_sieve(int n, bool prime_array[n], bool prime_power_array[n])
{
    // set initial array values
    for (int i = 0; i < n; i++)
    {
        prime_power_array[i] = false;
        prime_array[i] = true;
    }

    prime_array[0] = false;
    prime_array[1] = false;

    // regular sieve to find all primes up to n
    for (int p = 2; p * p < n; p++)
    {
        if (prime_array[p])
        {
            for (int i = p * p; i < n; i += p)
            {
                prime_array[i] = false;
            }
        }
    }

    // for each prime number, mark all multiples of it and itself as prime powers
    for (int i = 2; i < n; i++)
    {
        if (prime_array[i])
        {
            prime_power_array[i] = true;

            // Use long long to prevent integer overflow
            long long prime_power = (long long)i * i;
            while (prime_power < n)
            {
                // Safely mark prime powers within array bounds
                prime_power_array[prime_power] = true;
                // Check if next multiplication would cause overflow or exceed n
                if (prime_power > (long long)n / i)
                    break;
                prime_power *= i;
            }
        }
    }
}

// the following code is not used anywhere, but the algorithm was included in the thesis as an example
// constants used to print colours to stdout (not used in algorithm)
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define BLUE    "\033[34m"
#define WHITE   "\033[37m"

// returns the weight function for a partial colouring
// colours is a 2d square array of integers of size n * n
double getWeight(int n, int (*colours)[n])
{
    double w = 0.0;
    int numRed, numBlue;
    // loop over every k4:
    for (int v1 = 0; v1 < n - 3; v1++)
    {
        for (int v2 = v1 + 1; v2 < n - 2; v2++)
        {
            for (int v3 = v2 + 1; v3 < n - 1; v3++)
            {
                for (int v4 = v3 + 1; v4 < n; v4++)
                {
                    numRed = 0;
                    numBlue = 0;
                    // count the number of red and blue edges in this K4
                    if (colours[v1][v2] == 1)
                    {
                        numRed++;
                    } else if (colours[v1][v2] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v1][v3] == 1)
                    {
                        numRed++;
                    } else if (colours[v1][v3] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v1][v4] == 1)
                    {
                        numRed++;
                    } else if (colours[v1][v4] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v2][v3] == 1)
                    {
                        numRed++;
                    } else if (colours[v2][v3] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v2][v4] == 1)
                    {
                        numRed++;
                    } else if (colours[v2][v4] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v3][v4] == 1)
                    {
                        numRed++;
                    } else if (colours[v3][v4] == 2)
                    {
                        numBlue++;
                    } else {}

                    // determine the weight for this particular K4 now that we know the colouring
                    if (numBlue == 0 && numRed == 0)
                    {
                        w += pow(2, -5);
                    } else if (numBlue == 0)
                    {
                        w += pow(2, numRed - 6);
                    } else if (numRed == 0)
                    {
                        w += pow(2, numBlue - 6);
                    }
                }
            }
        }
    }
    return w;
}

// RED = 1 BLUE = 2 uncoloured = 0
int findColouring(int n)
{
    int colours[n][n];
    // set every colour to nothing (uncoloured) at first
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            colours[i][j] = 0;
        }
    }

    // find colouring by looping over each edge
    // (loop over each pair of vertices to get each edge)
    double wRed, wBlue; // weights for the graph if edge {i, j} is coloured red or blue
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            // calculate the weight for this edge if it is blue
            colours[i][j] = 2;
            colours[j][i] = 2;
            wBlue = getWeight(n, colours);

            // calculate the weight for this edge if it is red
            colours[i][j] = 1;
            colours[j][i] = 1;
            wRed = getWeight(n, colours);

            // the edge is currently set to red
            // if blue was better weight, change it to blue
            if (wBlue < wRed)
            {
                colours[i][j] = 2;
                colours[j][i] = 2;
            }

        }
    }

    // print the edges colours
    printf("Colours for n=%d\n", n);
    for (int i = 0; i < n; i ++)
    {
        for (int j = 0; j < n; j++)
        {
            if (colours[i][j] == 1)
            {
                printf(RED "R " RESET);
            } else if (colours[i][j] == 2)
            {
                printf(BLUE "B " RESET);
            }
            else
            {
                printf(WHITE "- " RESET);
            }
        }
        printf("\n");
    }

    // count number of monochromatic K4
    int numMonochromatic = 0;
    int numRed, numBlue;
    for (int v1 = 0; v1 < n - 3; v1++)
    {
        for (int v2 = v1 + 1; v2 < n - 2; v2++)
        {
            for (int v3 = v2 + 1; v3 < n - 1; v3++)
            {
                for (int v4 = v3 + 1; v4 < n; v4++)
                {
                    // determine if this K4 is monochromatic
                    numRed = 0;
                    numBlue = 0;
                    if (colours[v1][v2] == 1)
                    {
                        numRed++;
                    } else if (colours[v1][v2] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v1][v3] == 1)
                    {
                        numRed++;
                    } else if (colours[v1][v3] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v1][v4] == 1)
                    {
                        numRed++;
                    } else if (colours[v1][v4] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v2][v3] == 1)
                    {
                        numRed++;
                    } else if (colours[v2][v3] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v2][v4] == 1)
                    {
                        numRed++;
                    } else if (colours[v2][v4] == 2)
                    {
                        numBlue++;
                    } else {}
                    if (colours[v3][v4] == 1)
                    {
                        numRed++;
                    } else if (colours[v3][v4] == 2)
                    {
                        numBlue++;
                    } else {}

                    if (numBlue == 0 || numRed == 0)
                    {
                        numMonochromatic++;
                    }
                }
            }
        }
    }

    return numMonochromatic;
}

// calculates the upper bound on number of monochromatic K4 in Kn that's guaranteed by prop 16.1.1
double getExpectedNumberOfMonochromatics(int n)
{
    return (
          ((double) choose(n, 4))
        * pow(2, -5)
    );
}