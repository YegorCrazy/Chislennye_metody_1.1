#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Matrix {
    unsigned m, n; //матрица, m строк, n столбцов
    long double **elem;
} matrix;

long double abs_d (long double x) {
    if (x < 0) return -x;
    return x;
}

matrix *new_matrix (unsigned m, unsigned n) { //создание матрицы
    matrix *res = malloc(sizeof(matrix));
    if (res == NULL) return res;
    res->m = m;
    res->n = n;
    long double **elem = malloc(sizeof(long double *) * n);
    for (unsigned i = 0; i < n; ++i) {
        elem[i] = calloc(m, sizeof(long double)); // память по столбцам для удобного выбора ведущего элемента
    }
    res->elem = elem;
    return res;
}

matrix *copy_matrix (matrix *matr) {
    unsigned m = matr->m;
    unsigned n = matr->n;
    matrix *new = new_matrix(matr->m, matr->n);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            (new->elem)[i][j] = (matr->elem)[i][j];
        }
    }
    return new;
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
            printf("%Lf\t", (matr->elem)[j][i]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix_ext (matrix *A, matrix *f) { //вывод матрицы
    unsigned m = A->m;
    unsigned n = A->n;
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            printf("%Lf\t", (A->elem)[j][i]);
        }
        printf("%Lf\n", (f->elem)[0][i]);
    }
    printf("\n");
}

void read_matrix (FILE *input, matrix *matr) { //чтение матрицы из файла
    unsigned m = matr->m;    
    unsigned n = matr->n;
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            fscanf(input, "%Lf", &((matr->elem)[j][i]));
        }
    }
}

void fill_matrix1 (matrix *A, matrix *f, unsigned n) {
    if (n != 20) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double m1 = 8, n1 = 20;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = (i + j + 2) / (m1 + n1);
            } else {
                (A->elem)[i][j] = n1 + (m1 * m1) + ((j + 1) / m1) + ((i + 1) / n1);
            }
        }
    }
    for (unsigned i = 0; i < n; ++i) {
        (f->elem)[0][i] = 200 + (50 * (i + 1));
    }
}

void fill_matrix2 (matrix *A, matrix *f, unsigned n) {
    if (n != 25) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double m1 = 10, n1 = 25;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = (i + j + 2) / (m1 + n1);
            } else {
                (A->elem)[i][j] = n1 + (m1 * m1) + ((j + 1) / m1) + ((i + 1) / n1);
            }
        }
    }
    for (unsigned i = 0; i < n; ++i) {
        (f->elem)[0][i] = (i + 1) * (i + 1) - n1;
    }
}

long double q_M (long double M) {
    return 1.001 - (2 * M * 0.001);
}

void fill_matrix3 (matrix *A, matrix *f, unsigned n) {
    if (n != 100) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double M1 = 4, n1 = 100;
    long double qM = q_M(M1);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = pow(qM, (double)(i + j + 2)) + 0.1 * ((double)j - (double)i);
            } else {
                (A->elem)[i][j] = pow(abs_d(qM - 1), (double)(i + j + 2));
                if ((i + j) % 2 == 1) (A->elem)[i][j] *= (-1);
            }
        }
    }
    for (unsigned i = 0; i < n; ++i) {
        (f->elem)[0][i] = n1 * exp(M1 / (i + 1)) * cos(M1);
    }
}

void fill_matrix4 (matrix *A, matrix *f, unsigned n) {
    if (n != 100) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double M1 = 5, n1 = 100;
    long double qM = q_M(M1);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = pow(qM, (double)(i + j + 2)) + 0.1 * ((double)j - (double)i);
            } else {
                (A->elem)[i][j] = pow(abs_d(qM - 1), (double)(i + j + 2));
                if ((i + j) % 2 == 1) (A->elem)[i][j] *= (-1);
            }
        }
    }
    for (unsigned i = 0; i < n; ++i) {
        (f->elem)[0][i] = abs_d(M1 - (n1 / 10)) * i * sin(M1);
    }
}

void fill_matrix5 (matrix *A, matrix *f, unsigned n) {
    if (n != 100) {
        printf("ERROR ERROR ERROR\n");
        return;
    }
    long double M1 = 6;
    long double qM = q_M(M1);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (i != j) {
                (A->elem)[i][j] = pow(qM, (double)(i + j + 2)) + 0.1 * ((double)j - (double)i);
            } else {
                (A->elem)[i][j] = pow(abs_d(qM - 1), (double)(i + j + 2));
                if ((i + j) % 2 == 1) (A->elem)[i][j] *= (-1);
            }
        }
    }
    for (unsigned i = 0; i < n; ++i) {
        (f->elem)[0][i] = M1 * exp(M1 / (i + 1)) * cos(M1 / (i + 1));
    }
}

unsigned leading_element (matrix *matr, unsigned n_st, unsigned m_aim) { //менять местами столбцы для нахождения главного элемента
    long double lead_el = 0;
    unsigned lead_indx = 0;
    unsigned n_fin = matr->n;
    for (unsigned i = n_st; i < n_fin; ++i) {
        if (abs_d((matr->elem)[i][m_aim]) > lead_el) {
            lead_el = abs((matr->elem)[i][m_aim]);
            lead_indx = i;
        }
    }
    long double *tmp = (matr->elem)[n_st];
    (matr->elem)[n_st] = (matr->elem)[lead_indx];
    (matr->elem)[lead_indx] = tmp;
    return lead_indx;
}

void triangulate_matrix (matrix *matr, matrix *f) { //треугольный вид без выделения ведущего элемента
    unsigned m = matr->m;
    unsigned n = matr->n;
    long double **elem = matr->elem;
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = i + 1; j < m; ++j) {
            long double koef = elem[i][j] / elem[i][i];
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
    long double **elem = matr->elem;
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
            long double koef = elem[i][j] / elem[i][i];
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

long double determinant (matrix *matr) { //поиск определителя
    unsigned n = matr->n;
    matrix *new = copy_matrix(matr);
    triangulate_matrix(new, NULL);
    long double res = 1.0;
    for (unsigned i = 0; i < n; ++i) {
        res *= (new->elem)[i][i];
    }
    delete_matrix(new);
    return res;
}

long double *gauss_method (matrix *A1, matrix *f1, int lead) { //метод Гаусса
    matrix *A = copy_matrix(A1);
    matrix *f = copy_matrix(f1);
    unsigned *arr = NULL;
    if (lead == 0) {
        triangulate_matrix(A, f);
    } else {
        arr = triangulate_matrix_lead(A, f);
    }
    unsigned m = A->m;
    unsigned n = A->n;
    long double *roots_raw = calloc(n, sizeof(long double)); //тут немного разнятся коэффиценты, но при m = n все хорошо
    for (unsigned i = m - 1; ; --i) {
        roots_raw[i] = ((f->elem)[0][i]);
        for (unsigned j = i + 1; j < n; ++j) {
            roots_raw[i] -= roots_raw[j] * (A->elem)[j][i];
        }
        roots_raw[i] /= (A->elem)[i][i];
        if (i == 0) break;
    }
    delete_matrix(A);
    delete_matrix(f);
    if (lead == 0) return roots_raw;
    long double *roots = calloc(n, sizeof(long double));
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
