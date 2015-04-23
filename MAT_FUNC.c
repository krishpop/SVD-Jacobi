#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>

void swap_rows(double * i_mat, int row1, int row2, int n) {
    i_mat[row1 * n + row1] = 0.0;
    i_mat[row1 * n + row2] = 1.0;
    i_mat[row2 * n + row2] = 0.0;
    i_mat[row2 * n + row1] = 1.0;
}

void multiply(double * m1, double * m2, double * new_matrix, int n) {
    int new_row, new_column, old_x;

    for (new_row = 0; new_row < n; new_row++) {

        for (new_column = 0; new_column < n; new_column++) {
            new_matrix[new_row*n + new_column] = 0;

            for (old_x = 0; old_x < n; old_x++) {
                new_matrix[new_row*n + new_column] += m1[new_row*n + old_x] * m2[old_x*n + new_column];
            }
        }
    }
}

// Copies a to b
void copy_matrix(double * a, double * b, int n) {
    for (int i = 0; i < n*n; i++) {
        b[i] = a[i];
    }
}

void identity_matrix(double ** m, int n) {
    int i;
    *m = (double*)calloc(sizeof(double), n*n);
    for (i=0; i<n; i++) {
        (*m)[i*n+i] = 1.0;
    }
}


void transpose(double * a, int n) {
    double temp;
    int row, column;

    for (row = 0; row < n; row++) {
        for (column = (row + 1); column < n; column++) {
            if (row != column) {
              temp = a[(row * n) + column];
              a[(row * n) + column] = a[(column * n) + row];
              a[(column * n) + row] = temp;
            }
        }
    }
}

double* gen_matrix(int n) {
    double* a = malloc(n * n * sizeof(double));
        // int i, j;
        // double i_p_1, j_p_1;
    a[0] = 42.0;
    a[1] = 84.0;
    a[2] = 19.0;
    a[3] = 65.0;
    a[4] = 45.0;
    a[5] = 75.0;
    a[6] = 6.0;
    a[7] = 8.0;
    a[8] = 50.0;
    a[9] = 72.0;
    a[10] = 38.0;
    a[11] = 70.0;
    a[12] = 87.0;
    a[13] = 83.0;
    a[14] = 48.0;
    a[15] = 5.0;

    // for (i = 0; i < n; ++i)
    // {
    //     for (j = 0; j < n; ++j)
    //     {
    //         i_p_1 = ((double) i) + 1;
    //         j_p_1 = ((double) j) + 1;
    //         a[i * n + j] = sqrt(i_p_1 * i_p_1 + j_p_1 * j_p_1);
    //     }
    // }
    return a;
}

void print_matrix(double *m, int rows, int columns) {
    int row, column;
    for (row=0; row<rows; row++) {
        for (column=0; column<columns; column++) {
            printf("[%d, %d]: %.17f\t", row, column, m[row*rows+column]);
        }
        printf("\n");
    }
    printf("\n");
}
