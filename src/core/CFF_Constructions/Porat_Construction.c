#include "Porat_Construction.h"

#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "../CFF_Tables.h"
#include "../CFF.h"
#include "../CombinatorialAlgorithms.h"
#include "ID_Construction.h"
#include "Finite_Fields_Wrapper.h"


void poratShortSrcFormatter(char *strBuffer)
{
    strcpy(strBuffer, "Porat and Rothschild");
}

void poratLongSrcFormatter(short *consParams, char *strBuffer)
{
    sprintf(strBuffer, "PR($%hd$,$%hd$,$%hd$,$%hd$)", consParams[0], consParams[1], consParams[2], consParams[3]);
}

void poratConstructCFF(int d, int t)
{
    // this function is not implemented since porat doesn't appear in any tables until
    // n is larger than unsigned long long max value (so this function will never be called)
    global_tables_array[d-1]->array[t].cff = idMatrix(t, d);
}

CoverFreeFamily* genmatrixToCFF(GeneratorMatrixStruct gs)
{
    CoverFreeFamily* cff = allocateCFF(
        (gs.m - 1 ) / (gs.m - (gs.minDistance)),// d
        gs.m*gs.q, // t
        (int)pow(gs.q, gs.k) //n
    );
    int (*code)[cff->n] = (int (*)[cff->n]) gs.code;
    for (int codewordIndex = 0; codewordIndex < cff->n; codewordIndex++)
    {
        for (int codewordPosition = 0; codewordPosition < gs.m; codewordPosition++)
        {
            int codewordLetter = code[codewordPosition][codewordIndex];
            setCFFMatrixValue(cff, (codewordPosition * gs.q) + codewordLetter, codewordIndex, 1);
        }
    }
    return cff;
}

double poratEntropyFunction(double q, double r)
{
    return (
        (
            (r - 1) / r
            * log( ((q - 1) * r) / (r - 1) )
            / log(q)
        )
        + (
            1 / r
            * log(r)
            / log(q)
        )
    );
}

CFF_Construction_And_Name_Functions poratConstructionFunctions = {
    .shortSrcFormatter = poratShortSrcFormatter,
    .longSrcFormatter = poratLongSrcFormatter,
    .constructionFunction = poratConstructCFF
};

void addPoratCodes(CFF_Table* table, int cff_d, bool* prime_power_array)
{
    CFF_n_type cff_n;
    unsigned cff_t, r, k, m;
    double Hq;
    r = cff_d + 1;
    for (unsigned q = 2 * r; q < 4 * r && q < t_max; q++)
    {
        if (prime_power_array[q])
        {
            k = 1;
            cff_n = q;
            Hq = poratEntropyFunction((double) q, (double) r);
            m = ceil( ((double) (k)) / (1.0 - Hq) );
            cff_t = m * q;
            while (cff_t <= t_max && cff_n != 0)
            {
                updateTable(table, cff_t, cff_n, &poratConstructionFunctions, q,k,m,r,0);
                // change array to test next
                k++;
                m = (unsigned) ceil( ((double) (k)) / (1.0 - Hq) );
                cff_n = pow(q, k);
                cff_t = m * q;
            }
        }
    }
}

// reverse (index revere, not total order reverse) lexicographic k-tuple successor
bool nextReverseLex(int n, int k, int buffer[k])
{
    for (int i = 0; i < k; i++)
    {
        if (buffer[i] < n - 1)
        {
            buffer[i]++;
            for (int x = i-1; x > -1; x--)
            {
                buffer[x] = 0;
            }
            return true;
        }
    }
    return false;
}

void freeGeneratorMatrixStruct(GeneratorMatrixStruct genMatrixStruct)
{
    free(genMatrixStruct.generatorMatrix);
    free(genMatrixStruct.code);
}

// construct Porat generator matrix, innefficient brute force version, only when q is prime (not prime power)
GeneratorMatrixStruct constructGeneratorMatrixBruteMod(int q, int k, int r, int m)
{
    if (q < 2 * r || q >= 4 * r)
    {
        printf("Parameters do not match thm requirments. q not in [2r, 4r)\n");
    }
    // determine other parameters from q,k,d (cff's d)
    double delta = (((double) r) - 1.0) / ((double) r);
    double Hq = poratEntropyFunction((double) q, (double) r);
    if (m == 0)
    {
        m = (int) ceil(((double) (k)) / (1.0 - Hq));
    }
    int D = (int) floor(delta*m);

    // store q^k, since it's in a for loop bound
    int q_to_the_k = (int) pow(q, k);

    printf("Starting porat cons with: q=%d k=%d r=%d m=%d Hq(δ)=%f δ=%f Distance=%d\n",q,k,r,m,Hq,delta,D);
    printf("This gives a: %d-CFF(%d,%d)\n", r-1, m*q, q_to_the_k);

    // allocate memory for the generator matrix
    int* g_pointer = malloc(sizeof(int) * k * m);
    // declare 2d array pointer ([m * k] size), call it g (to access the generator matrix easier) and assign it g_pointer's address
    int (*G)[k] = (int (*)[k]) g_pointer;

    // allocate memory for the full code
    int* code_pointer = malloc(sizeof(int) * q_to_the_k * m);
    // 2d ptr to access more easily:
    int (*code)[q_to_the_k] = (int (*)[q_to_the_k]) code_pointer;

    // set zero codeword
    for(int i = 0; i < m; i++)
    {
        code[i][0] = 0;
    }

    // A[x] is the amount of codeword letters fixed to zero for codeword x
    int A[q_to_the_k];
    // set a to all zeros, except for the zero codeword (first codeword)
    A[0] = m;
    for (int i = 1; i < q_to_the_k; i++)
    {
        A[i] = 0;
    }
    // y = all lexicographic k-tuples of Fq
    int y[q_to_the_k][k];
    // now fill y in with a the successor algorithm
    // (but first fill k_tuple_buffer with the 1st lexicographic element)
    int k_tuple_buffer[k];
    for (int i = 0; i < k; i++)
    {
        k_tuple_buffer[i] = 0;
    }
    int h = 0; //h is the current index of y to fill in the following loop
    //printf("y={\n");
    do
    {
        for (int i = 0; i < k ; i++)
        {
            y[h][i] = k_tuple_buffer[i];
            //printf("%d, ",y[h][i]);
        }
        //printf("\n");
        h++;
    } while (nextReverseLex(q, k, k_tuple_buffer));
    //printf("}\n\n");

    // loop over each cell in generator matrix, start filling it in
    int R, c;
    double m_i_choose_x, one_minus_power, one_over_power;
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= k; j++)
        {
            int v_best = -1;
            double E_best = DBL_MAX;
            for (int v = 0; v < q; v++)
            {
                G[i-1][j-1] = v;
                double E = 0.0;
                for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
                {
                    c = i - 1 - A[l];
                    R = 0;
                    for (int p = 1; p <= j; p++)
                    {
                        R += (G[i-1][p-1] * y[l][p-1]) % q;
                        R %= q;
                    }
                    if (R != 0)
                    {
                        c++;
                    }
                    //printf("{[l=%d, A[%d]=%d]",l,l,A[l]);
                    for (int x = 0; x < D - c; x++)
                    {
                        m_i_choose_x = (double) choose(m - i, x);
                        one_minus_power = pow(1.0 - (1.0 / ((double) q)), (double) x);
                        one_over_power = pow(1.0 / ((double) q), (double) (m - i - x));
                        E += m_i_choose_x * one_minus_power * one_over_power;
                        //E +=   ((double) choose(m - i, x))
                        //     * (pow(1.0 - (1.0 / ((double) q)), (double) x))
                        //     * (pow(1.0 / ((double) q), (double) (m - i - x)));
                        //printf("(E+=%f*%f*%f. E=%f, x=%d)",m_i_choose_x,one_minus_power,one_over_power,E,x);
                    }
                    //printf("} ");
                }
                //printf("E=%f (i=%d,j=%d,v=%d,R=%d,c=%d) --- ", E,i,j,v,R,c);
                if (E < E_best)
                {
                    E_best = E;
                    v_best = v;
                }
            }
            //printf("E_best=%f,v_best=%d\n\n",E_best, v_best);
            G[i-1][j-1] = v_best;
            for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
            {
                R = 0;
                for (int p = 1; p <= j; p++)
                {
                    R += (G[i-1][p-1] * y[l][p-1]) % q;
                    R %= q;
                }
                if (R == 0)
                {
                    A[l]++;
                }
                code[i-1][l] = R;
            }
        }
    }

    GeneratorMatrixStruct genMatrixStruct;
    genMatrixStruct.q = q;
    genMatrixStruct.k = k;
    genMatrixStruct.r = r;
    genMatrixStruct.m = m;

    genMatrixStruct.Hq = Hq;
    genMatrixStruct.delta = delta;

    genMatrixStruct.generatorMatrix = g_pointer;
    genMatrixStruct.code = code_pointer;

    genMatrixStruct.numCodewords = q_to_the_k;
    genMatrixStruct.minDistance = D;

    return genMatrixStruct;
}

GeneratorMatrixStruct constructGeneratorMatrixBruteFF(int p, int a, int k, int r, int m)
{
    int q = compute_field_size(p,a);
    int add[q][q];
    int mult[q][q];
    int mult_inverses[q];
    int add_inverses[q];

    populateFiniteField(p,a,(int *)add,(int *)mult);
    populateAdditiveInverses(p, a,(int *) add, add_inverses);
    populateMultiplicativeInverses(p, a,(int *) mult, mult_inverses);
    if (q < 2 * r || q >= 4 * r)
    {
        printf("Parameters do not match thm requirments. q not in [2r, 4r)\n");
    }
    // determine other parameters from q,k,d (cff's d)
    double delta = (((double) r) - 1.0) / ((double) r);
    double Hq = poratEntropyFunction((double) q, (double) r);
    if (m == 0)
    {
        m = (int) ceil(((double) (k)) / (1.0 - Hq));
    }
    int D = (int) floor(delta*m);

    // store q^k, since it's in a for loop bound
    int q_to_the_k = (int) pow(q, k);

    printf("Starting porat cons with: q=%d k=%d r=%d m=%d Hq(δ)=%f δ=%f Distance=%d\n",q,k,r,m,Hq,delta,D);
    printf("This gives a: %d-CFF(%d,%d)\n", r-1, m*q, q_to_the_k);

    // allocate memory for the generator matrix
    int* g_pointer = malloc(sizeof(int) * k * m);
    // declare 2d array pointer ([m * k] size), call it g (to access the generator matrix easier) and assign it g_pointer's address
    int (*G)[k] = (int (*)[k]) g_pointer;

    // allocate memory for the full code
    int* code_pointer = malloc(sizeof(int) * q_to_the_k * m);
    // 2d ptr to access more easily:
    int (*code)[q_to_the_k] = (int (*)[q_to_the_k]) code_pointer;

    // set zero codeword
    for(int i = 0; i < m; i++)
    {
        code[i][0] = 0;
    }

    // A[x] is the amount of codeword letters fixed to zero for codeword x
    int A[q_to_the_k];
    // set a to all zeros, except for the zero codeword (first codeword)
    A[0] = m;
    for (int i = 1; i < q_to_the_k; i++)
    {
        A[i] = 0;
    }
    // y = all lexicographic k-tuples of Fq
    int y[q_to_the_k][k];
    // now fill y in with a the successor algorithm
    // (but first fill k_tuple_buffer with the 1st lexicographic element)
    int k_tuple_buffer[k];
    for (int i = 0; i < k; i++)
    {
        k_tuple_buffer[i] = 0;
    }
    int h = 0; //h is the current index of y to fill in the following loop
    //printf("y={\n");
    do
    {
        for (int i = 0; i < k ; i++)
        {
            y[h][i] = k_tuple_buffer[i];
            //printf("%d, ",y[h][i]);
        }
        //printf("\n");
        h++;
    } while (nextReverseLex(q, k, k_tuple_buffer));
    //printf("}\n\n");

    // loop over each cell in generator matrix, start filling it in
    int R, c;
    double m_i_choose_x, one_minus_power, one_over_power;
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= k; j++)
        {
            int v_best = -1;
            double E_best = DBL_MAX;
            for (int v = 0; v < q; v++)
            {
                G[i-1][j-1] = v;
                double E = 0.0;
                for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
                {
                    c = i - 1 - A[l];
                    R = 0;
                    for (int p = 1; p <= j; p++)
                    {
                        //R += (G[i-1][p-1] * y[l][p-1]) % q;
                        //R %= q;
                        R = add[mult[G[i-1][p-1]][y[l][p-1]]][R];
                    }
                    if (R != 0)
                    {
                        c++;
                    }
                    //printf("{[l=%d, A[%d]=%d]",l,l,A[l]);
                    for (int x = 0; x < D - c; x++)
                    {
                        m_i_choose_x = (double) choose(m - i, x);
                        one_minus_power = pow(1.0 - (1.0 / ((double) q)), (double) x);
                        one_over_power = pow(1.0 / ((double) q), (double) (m - i - x));
                        E += m_i_choose_x * one_minus_power * one_over_power;
                        //E +=   ((double) choose(m - i, x))
                        //     * (pow(1.0 - (1.0 / ((double) q)), (double) x))
                        //     * (pow(1.0 / ((double) q), (double) (m - i - x)));
                        //printf("(E+=%f*%f*%f. E=%f, x=%d)",m_i_choose_x,one_minus_power,one_over_power,E,x);
                    }
                    //printf("} ");
                }
                //printf("E=%f (i=%d,j=%d,v=%d,R=%d,c=%d) --- ", E,i,j,v,R,c);
                if (E < E_best)
                {
                    E_best = E;
                    v_best = v;
                }
            }
            //printf("E_best=%f,v_best=%d\n\n",E_best, v_best);
            G[i-1][j-1] = v_best;
            for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
            {
                R = 0;
                for (int p = 1; p <= j; p++)
                {
                    //R += (G[i-1][p-1] * y[l][p-1]) % q;
                    //R %= q;
                    R = add[mult[G[i-1][p-1]][y[l][p-1]]][R];
                }
                if (R == 0)
                {
                    A[l]++;
                }
                code[i-1][l] = R;
            }
        }
    }

    GeneratorMatrixStruct genMatrixStruct;
    genMatrixStruct.q = q;
    genMatrixStruct.k = k;
    genMatrixStruct.r = r;
    genMatrixStruct.m = m;

    genMatrixStruct.Hq = Hq;
    genMatrixStruct.delta = delta;

    genMatrixStruct.generatorMatrix = g_pointer;
    genMatrixStruct.code = code_pointer;

    genMatrixStruct.numCodewords = q_to_the_k;
    genMatrixStruct.minDistance = D;

    return genMatrixStruct;
}


__attribute__((noinline)) GeneratorMatrixStruct constructGeneratorMatrixBestBrute(int p, int a, int k, int r, int m)
{
    int q = compute_field_size(p,a);
    int add[q][q];
    int mult[q][q];
    //int mult_inverses[q];
    //int add_inverses[q];

    populateFiniteField(p,a,(int *)add,(int *)mult);
    //populateAdditiveInverses(p, a,(int *) add, add_inverses);
    //populateMultiplicativeInverses(p, a,(int *) mult, mult_inverses);
    if (q < 2 * r || q >= 4 * r)
    {
        printf("Parameters do not match thm requirments. q not in [2r, 4r)\n");
    }
    // determine other parameters from q,k,d (cff's d)
    double delta = (((double) r) - 1.0) / ((double) r);
    double Hq = poratEntropyFunction((double) q, (double) r);
    if (m == 0)
    {
        m = (int) ceil(((double) (k)) / (1.0 - Hq));
    }
    int D = (int) floor(delta*m);

    // store q^k, since it's in a for loop bound
    int q_to_the_k = (int) pow(q, k);

    printf("Starting porat cons with: q=%d k=%d r=%d m=%d Hq(δ)=%f δ=%f Distance=%d\n",q,k,r,m,Hq,delta,D);
    printf("This gives a: %d-CFF(%d,%d)\n", r-1, m*q, q_to_the_k);

    // allocate memory for the generator matrix
    int* g_pointer = malloc(sizeof(int) * k * m);
    // declare 2d array pointer ([m * k] size), call it g (to access the generator matrix easier) and assign it g_pointer's address
    int (*G)[k] = (int (*)[k]) g_pointer;

    // allocate memory for the full code
    int* code_pointer = malloc(sizeof(int) * q_to_the_k * m);
    // 2d ptr to access more easily:
    int (*code)[q_to_the_k] = (int (*)[q_to_the_k]) code_pointer;

    // set zero codeword
    for(int i = 0; i < m; i++)
    {
        code[i][0] = 0;
    }

    // A[x] is the amount of codeword letters fixed to zero for codeword x
    int A[q_to_the_k];
    // set a to all zeros, except for the zero codeword (first codeword)
    A[0] = m;
    for (int i = 1; i < q_to_the_k; i++)
    {
        A[i] = 0;
    }
    // y = all lexicographic k-tuples of Fq
    int y[q_to_the_k][k];
    // now fill y in with a the successor algorithm
    // (but first fill k_tuple_buffer with the 1st lexicographic element)
    int k_tuple_buffer[k];
    for (int i = 0; i < k; i++)
    {
        k_tuple_buffer[i] = 0;
    }
    int h = 0; //h is the current index of y to fill in the following loop
    //printf("y={\n");
    do
    {
        for (int i = 0; i < k ; i++)
        {
            y[h][i] = k_tuple_buffer[i];
            //printf("%d, ",y[h][i]);
        }
        //printf("\n");
        h++;
    } while (nextReverseLex(q, k, k_tuple_buffer));
    //printf("}\n\n");

    // loop over each cell in generator matrix, start filling it in
    int R, c;
    double m_i_choose_x, one_minus_power, one_over_power;
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= k; j++)
        {
            printf("filling G[i=%d][j=%d]\n",i,j);
            int v_best = -1;
            double E_best = DBL_MAX;
            for (int v = 0; v < q; v++)
            {
                G[i-1][j-1] = v;
                double E = 0.0;
                for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
                {
                    c = i - 1 - A[l];
                    R = add
                        [mult
                            [G[i-1][j-1]]
                            [y[l][j-1]]
                        ]
                        [
                            code[i-1][(l%((int)pow(q,j-1)))]
                        ];
                    if (R != 0)
                    {
                        c++;
                    }
                    for (int x = 0; x < D - c; x++)
                    {
                        m_i_choose_x = (double) choose(m - i, x);
                        one_minus_power = pow(1.0 - (1.0 / ((double) q)), (double) x);
                        one_over_power = pow(1.0 / ((double) q), (double) (m - i - x));
                        E += m_i_choose_x * one_minus_power * one_over_power;
                    }
                }
                printf("\t Letter %d has E=%.17g\n",v,E);
                if (E < E_best)
                {
                    E_best = E;
                    v_best = v;
                }
            }
            G[i-1][j-1] = v_best;
            printf("\t Letter %d was chosen as best\n", v_best);
            for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
            {
                R = add
                        [mult
                            [G[i-1][j-1]]
                            [y[l][j-1]]
                        ]
                        [
                            code[i-1][(l%((int)pow(q,j-1)))]
                        ];
                if (R == 0)
                {
                    A[l]++;
                }
                code[i-1][l] = R;
            }
        }
    }

    GeneratorMatrixStruct genMatrixStruct;
    genMatrixStruct.q = q;
    genMatrixStruct.k = k;
    genMatrixStruct.r = r;
    genMatrixStruct.m = m;

    genMatrixStruct.Hq = Hq;
    genMatrixStruct.delta = delta;

    genMatrixStruct.generatorMatrix = g_pointer;
    genMatrixStruct.code = code_pointer;

    genMatrixStruct.numCodewords = q_to_the_k;
    genMatrixStruct.minDistance = D;

    return genMatrixStruct;
}

__attribute__((noinline)) GeneratorMatrixStruct constructGeneratorMatrixBestOptimized(int p, int a, int k, int r, int m)
{
    int q = compute_field_size(p,a);
    int add[q][q];
    int mult[q][q];
    int mult_inverses[q];
    int add_inverses[q];

    populateFiniteField(p,a,(int *)add,(int *)mult);
    populateAdditiveInverses(p, a,(int *) add, add_inverses); //remove this linear search
    populateMultiplicativeInverses(p, a,(int *) mult, mult_inverses);
    if (q < 2 * r || q >= 4 * r)
    {
        printf("Parameters do not match thm requirments. q not in [2r, 4r)\n");
    }
    // determine other parameters from q,k,d (cff's d)
    double delta = (((double) r) - 1.0) / ((double) r);
    double Hq = poratEntropyFunction((double) q, (double) r);
    if (m == 0)
    {
        m = (int) ceil(((double) (k)) / (1.0 - Hq));
    }
    int D = (int) floor(delta*m);

    // store q^k, since it's in a for loop bound
    int q_to_the_k = (int) pow(q, k);

    printf("Starting porat cons with: q=%d k=%d r=%d m=%d Hq(δ)=%f δ=%f Distance=%d\n",q,k,r,m,Hq,delta,D);
    printf("This gives a: %d-CFF(%d,%d)\n", r-1, m*q, q_to_the_k);

    // allocate memory for the generator matrix
    int* g_pointer = malloc(sizeof(int) * k * m);
    // declare 2d array pointer ([m * k] size), call it g (to access the generator matrix easier) and assign it g_pointer's address
    int (*G)[k] = (int (*)[k]) g_pointer;

    // allocate memory for the full code
    int* code_pointer = malloc(sizeof(int) * q_to_the_k * m);
    // 2d ptr to access more easily:
    int (*code)[q_to_the_k] = (int (*)[q_to_the_k]) code_pointer;

    // set zero codeword
    for(int i = 0; i < m; i++)
    {
        code[i][0] = 0;
    }

    // A[x] is the amount of codeword letters fixed to zero for codeword x
    int A[q_to_the_k];
    // set a to all zeros, except for the zero codeword (first codeword)
    A[0] = m;
    for (int i = 1; i < q_to_the_k; i++)
    {
        A[i] = 0;
    }
    // y = all lexicographic k-tuples of Fq
    int y[q_to_the_k][k];
    // now fill y in with a the successor algorithm
    // (but first fill k_tuple_buffer with the 1st lexicographic element)
    int k_tuple_buffer[k];
    for (int i = 0; i < k; i++)
    {
        k_tuple_buffer[i] = 0;
    }
    int h = 0; //h is the current index of y to fill in the following loop
    do
    {
        for (int i = 0; i < k ; i++)
        {
            y[h][i] = k_tuple_buffer[i];
        }
        h++;
    } while (nextReverseLex(q, k, k_tuple_buffer));

    // loop over each cell in generator matrix, start filling it in
    int c, v;
    double m_i_choose_x, one_minus_power, one_over_power;
    double E_best;
    double W[q];
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= k; j++)
        {
            printf("filling G[i=%d][j=%d]\n",i,j);
            for (int b = 0; b < q; b++)
            {
                W[b] = 0;
            }
            for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
            {
                c = i - A[l];
                v = mult
                    [code[i-1][(l%((int)pow(q,j-1)))]]
                    [add_inverses[mult_inverses[y[l][j-1]]]];
                if ( (D-c)>-1)
                {
                    m_i_choose_x = (double) choose(m - i, D-c);
                    one_minus_power = pow(1.0 - (1.0 / ((double) q)), (double) (D-c));
                    one_over_power = pow(1.0 / ((double) q), (double) ((m - i) - (D-c)));
                    W[v] -= m_i_choose_x * one_minus_power * one_over_power;
                }
            }
            E_best = W[0];
            int v_best = 0;
            for (int b = 0; b < q; b++)
            {
                printf("\tW[%d] = %.17g\n", b, W[b]);
                if (W[b] > E_best)
                {
                    v_best = b;
                    E_best = W[b];
                }
            }
            printf("\tLetter %d was chosen as best\n", v_best);
            G[i-1][j-1] = v_best;
            for (int l = (int) pow(q, j-1); l < (int) pow(q, j); l++)
            {
                v = add
                        [mult
                            [G[i-1][j-1]]
                            [y[l][j-1]]
                        ]
                        [
                            code[i-1][(l%((int)pow(q,j-1)))]
                        ];
                if (v == 0)
                {
                    A[l]++;
                }
                code[i-1][l] = v;
            }
        }
    }

    GeneratorMatrixStruct genMatrixStruct;
    genMatrixStruct.q = q;
    genMatrixStruct.k = k;
    genMatrixStruct.r = r;
    genMatrixStruct.m = m;

    genMatrixStruct.Hq = Hq;
    genMatrixStruct.delta = delta;

    genMatrixStruct.generatorMatrix = g_pointer;
    genMatrixStruct.code = code_pointer;

    genMatrixStruct.numCodewords = q_to_the_k;
    genMatrixStruct.minDistance = D;

    return genMatrixStruct;
}