#include <stdio.h>

typedef struct Matrix {
    unsigned m, n;
    long double **elem;
} matrix;

long double abs_d (long double x);

matrix *new_matrix (unsigned m, unsigned n);
matrix *copy_matrix (matrix *matr);
void delete_matrix (matrix *matr);
void print_matrix (matrix *matr);
void print_matrix_ext (matrix *A, matrix *f);
void read_matrix (FILE *input, matrix *matr);
void fill_matrix1 (matrix *A, matrix *f, unsigned n);
void fill_matrix2 (matrix *A, matrix *f, unsigned n);
long double q_M (long double M);
void fill_matrix3 (matrix *A, matrix *f, unsigned n);
void fill_matrix4 (matrix *A, matrix *f, unsigned n);
void fill_matrix5 (matrix *A, matrix *f, unsigned n);

unsigned leading_element (matrix *matr, unsigned n_st, unsigned m_aim);
void triangulate_matrix (matrix *matr, matrix *f, int *reverse);
unsigned *triangulate_matrix_lead (matrix *matr, matrix *f);

long double determinant (matrix *matr);
long double *gauss_method (matrix *A1, matrix *f1, int lead);

matrix *reverse_matrix (matrix *A2);
long double matrix_norm (matrix *A);
