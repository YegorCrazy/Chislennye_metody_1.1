#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv) {
    unsigned n;
    scanf("%u", &n);
    
    matrix *A = new_matrix(n, n);
    matrix *f = new_matrix(n, 1);
    FILE *input = fopen("input.txt", "r");
    read_matrix(input, A);
    read_matrix(input, f);
    fclose(input);
    
    printf("Determinant is %lf\n\n", determinant(A));
    
    double *arr = gauss_method(A, f, 0);
    for (unsigned i = 0; i < n; ++i) {
        printf("%lf\t", arr[i]);
    }
    printf("\n\n");
    free(arr);
    
    input = fopen("input.txt", "r");
    read_matrix(input, A);
    read_matrix(input, f);
    fclose(input);
    arr = gauss_method(A, f, 1);
    for (unsigned i = 0; i < n; ++i) {
        printf("%lf\t", arr[i]);
    }
    printf("\n");
    free(arr);
    
    delete_matrix(A);
    delete_matrix(f);
    
    return 0;
}
