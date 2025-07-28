// small test to check if finite fields are precomputed properly
#include "flint/flint.h"
#include "flint/fq.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_factor.h"

#include "core/CFF_Constructions/Finite_Fields_Wrapper.h"
#include "core/CFF.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{

    if (argc != 3) {
        printf("give 2 parameters: p,k\n");
        return 1;
    }

    // Convert arguments to integers using atoi()
    int p = atoi(argv[1]);
    int k = atoi(argv[2]);

    int size = (int) pow(p,k);
    int add[size][size];
    int mult[size][size];
    int mult_inverses[size];
    int add_inverses[size];

    populateFiniteField(p,k,(int *)add,(int *)mult);
    populateAdditiveInverses(p, k,(int *) add, add_inverses);
    populateMultiplicativeInverses(p, k,(int *) mult, mult_inverses);
    printf("Addition table:\n");
    for(int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%d ", add[i][j]);
        }
        printf("\n");
    }
    printf("Multiplication table:\n");
    for(int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%d ", mult[i][j]);
        }
        printf("\n");
    }
    printf("Additive inverses:\n");
    for (int i = 0; i < size; i++)
    {
        printf("%d : %d\n", i, add_inverses[i]);
    }
    printf("Multiplicative inverses:\n");
    for (int i = 0; i < size; i++)
    {
        printf("%d : %d\n", i, mult_inverses[i]);
    }

    return 0;
}
