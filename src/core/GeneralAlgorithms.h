#ifndef GENERAL_ALGORITHMS_HEADER
#define GENERAL_ALGORITHMS_HEADER

#include <stdbool.h>
int int_comparator(const void *a, const void *b);
int ipow(int base, int exp);
int compare_value_rank_pair(const void *a, const void *b);
int euclid_GCD(int a, int b);
void rankify_and_sort
    (
        int array_length,
        int input_unsorted_array[array_length],
        int output_rankified_buffer[array_length],
        int output_sorted_buffer[array_length]
    );
void prime_power_sieve(int n, bool prime_array[n], bool prime_power_array[n]);
#endif