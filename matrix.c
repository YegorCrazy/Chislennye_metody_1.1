#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Matrix {
    unsigned m, n; //матрица, m строк, n столбцов
    double **elem;
} matrix;

double abs_d (double x) {
    if (x < 0) return -x;
    return x;
}

matrix *new_matrix (unsigned m, unsigned n) { //создание матрицы
    matrix *res = malloc(sizeof(matrix));
    if (res == NULL) return res;
    res->m = m;
    res->n = n;
    double **elem = malloc(sizeof(double *) * n);
    for (unsigned i = 0; i < n; ++i) {
        elem[i] = calloc(m, sizeof(double)); // память по столбцам для удобного выбора ведущего элемента
    }
    res->elem = elem;
    return res;
}

void delete_matrix (matrix *matr) { //удаление матрицы
    unsigned n = matr->n;
    for (unsigned i = 0; i < n; ++i) {
        free((matr->elem)[i]);
    }
    free(matr->elem);
    free(matr);
}

void print_matrix (matrix *matr) { //вывод матрицы
    unsigned m = matr->m;
    unsigned n = matr->n;
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            printf("%lf\t", (matr->elem)[j][i]);
        }
        printf("\n");
    }
}

void read_matrix (FILE *input, matrix *matr) { //чтение матрицы из файла
    unsigned m = matr->m;    
    unsigned n = matr->n;
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            fscanf(input, "%lf", &((matr->elem)[j][i]));
        }
    }
}

unsigned leading_element (matrix *matr, unsigned n_st, unsigned m_aim) { //менять местами столбцы для нахождения главного элемента
    double lead_el = 0;
    unsigned lead_indx = 0;
    unsigned n_fin = matr->n;
    for (unsigned i = n_st; i < n_fin; ++i) {
        if (abs_d((matr->elem)[i][m_aim]) > lead_el) {
            lead_el = abs((matr->elem)[i][m_aim]);
            lead_indx = i;
        }
    }
    double *tmp = (matr->elem)[n_st];
    (matr->elem)[n_st] = (matr->elem)[lead_indx];
    (matr->elem)[lead_indx] = tmp;
    return lead_indx;
}

void triangulate_matrix (matrix *matr, matrix *f) { //треугольный вид без выделения ведущего элемента
    unsigned m = matr->m;
    unsigned n = matr->n;
    double **elem = matr->elem;
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = i + 1; j < m; ++j) {
            double koef = elem[i][j] / elem[i][i];
            for (unsigned k = i; k < n; ++k) {
                elem[k][j] -= elem[k][i] * koef;
            }
            if (f != NULL) {
                (f->elem)[0][j] -= (f->elem)[0][i] * koef;
            }
        }
    }
}
    
unsigned *triangulate_matrix_lead (matrix *matr, matrix *f) { //треугольный вид с выделением ведущего элемента
    unsigned m = matr->m;
    unsigned n = matr->n;
    double **elem = matr->elem;
    unsigned swp;
    unsigned *arr = malloc(n * sizeof(unsigned));
    for (unsigned i = 0; i < n; ++i) {
        arr[i] = i;
    }
    for (unsigned i = 0; i < m; ++i) {
        swp = leading_element(matr, i, i);
        unsigned tmp = arr[i];
        arr[i] = arr[swp];
        arr[swp] = tmp;
        for (unsigned j = i + 1; j < m; ++j) {
            double koef = elem[i][j] / elem[i][i];
            for (unsigned k = i; k < n; ++k) {
                elem[k][j] -= elem[k][i] * koef;
            }
            if (f != NULL) {
                (f->elem)[0][j] -= (f->elem)[0][i] * koef;
            }
        }
    }
    return arr;
}

double determinant (matrix *matr) { //поиск определителя
    unsigned m = matr->m;
    unsigned n = matr->n;
    matrix *new = new_matrix(matr->m, matr->n);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            (new->elem)[i][j] = (matr->elem)[i][j];
        }
    }
    triangulate_matrix(new, NULL);
    double res = 1.0;
    for (unsigned i = 0; i < n; ++i) {
        res *= (new->elem)[i][i];
    }
    delete_matrix(new);
    return res;
}

double *gauss_method (matrix *A, matrix *f, int lead) { //метод Гаусса
    unsigned *arr = NULL;
    if (lead == 0) {
        triangulate_matrix(A, f);
    } else {
        arr = triangulate_matrix_lead(A, f);
    }
    unsigned m = A->m;
    unsigned n = A->n;
    double *roots_raw = calloc(n, sizeof(double)); //тут немного разнятся коэффиценты, но при m = n все хорошо
    for (unsigned i = m - 1; ; --i) {
        roots_raw[i] = ((f->elem)[0][i]);
        for (unsigned j = i + 1; j < n; ++j) {
            roots_raw[i] -= roots_raw[j] * (A->elem)[j][i];
        }
        roots_raw[i] /= (A->elem)[i][i];
        if (i == 0) break;
    }
    if (lead == 0) return roots_raw;
    double *roots = calloc(n, sizeof(double));
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (arr[j] == i) {
                roots[i] = roots_raw[j];
                break;
            }
        }
    }
    free(arr);
    free(roots_raw);
    return roots;
}
