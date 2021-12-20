#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv) {
    unsigned n;
    scanf("%u", &n);
    //n = 20;
    
    matrix *A = new_matrix(n, n);
    matrix *f = new_matrix(n, 1);
    
    FILE *input = fopen("input.txt", "r");
    read_matrix(input, A);
    read_matrix(input, f);
    fclose(input);
    //fill_matrix1(A, f, n);
    
    //print_matrix(A);
    printf("Starting f coefficents are:\n\n");
    print_matrix(f);
    
    /*matrix *new = copy_matrix(A);
    triangulate_matrix(new, NULL);
    print_matrix(new);
    delete_matrix(new);
    
    new = copy_matrix(A);
    unsigned *del = triangulate_matrix_lead(new, NULL);
    print_matrix(new);
    delete_matrix(new);
    free(del);*/
    
    long double det = determinant(A);
    
    printf("Determinant is %Lf (if it is 0, result can't be calculated well)\n\n", det);
    
    if (abs_d(det) < 0.00000001) {
        printf("Determinant is nearly zero, condition number can't be calculated\n\n");
    } else {
        matrix *reverse = reverse_matrix(A);
        long double norm = matrix_norm(A) * matrix_norm(reverse);
        printf("Condition number counted with column norm is %Lf\n\n", norm);
        delete_matrix(reverse);
    }
    
    long double *arr = gauss_method(A, f, 0);
    printf("Roots after default Gauss method are:\n\n");
    for (unsigned i = 0; i < n; ++i) {
        printf("x%u = %Lf\t", i + 1, arr[i]); //метод Гаусса без ведущего элемента
    }
    printf("\n\n");
    printf("Right part with higher roots counted is:\n\n");
    for (unsigned i = 0; i < n; ++i) { //m
        long double cur_res = 0;
        for (unsigned j = 0; j < n; ++j) { //n
            cur_res += ((A->elem)[j][i] * arr[j]);
        }
        printf("%Lf\n", cur_res);
    }
    printf("\n");
    free(arr);
    
    arr = gauss_method(A, f, 1);
    printf("Roots after Gauss method with leading element search are:\n\n");
    for (unsigned i = 0; i < n; ++i) {
        printf("x%u = %Lf\t", i + 1, arr[i]); // метод Гаусса с выбором ведущего элемента
    }
    printf("\n\n");
    printf("Right part with higher roots counted is:\n\n");
    for (unsigned i = 0; i < n; ++i) { //m
        long double cur_res = 0;
        for (unsigned j = 0; j < n; ++j) { //n
            cur_res += ((A->elem)[j][i] * arr[j]);
        }
        printf("%Lf\n", cur_res);
    }
    printf("\n");
    free(arr);
    
    delete_matrix(A);
    delete_matrix(f);
    
    return 0;
}
