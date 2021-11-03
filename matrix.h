#include <stdio.h>

typedef struct Matrix {
    unsigned m, n;
    double **elem;
} matrix;

double abs_d (double x);

matrix *new_matrix (unsigned m, unsigned n);
void delete_matrix (matrix *matr);
void print_matrix (matrix *matr);
void read_matrix (FILE *input, matrix *matr);

unsigned leading_element (matrix *matr, unsigned n_st, unsigned m_aim);
void triangulate_matrix (matrix *matr, matrix *f);
unsigned *triangulate_matrix_lead (matrix *matr, matrix *f);

double determinant (matrix *matr);
double *gauss_method (matrix *A, matrix *f, int flag);
