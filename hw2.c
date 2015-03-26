//
//  CPSC 440 - hw2.c
//  Jacobi Polynomials
//
//  Created by Krishnan Srinivasan on 3/23/15.
//

// SVD using Jacobi Algorithm from:
// Brent, Luk, Loan, "Computation of the Singular Value Decomposition Using
// Mesh-Connected Processors"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>

#define ERROR(msg) exit (fprintf (stderr, "%s\n", msg))
#define LEFT 0
#define RIGHT 1

// matrix functions
void copy_matrix(double * a, double * b, int n);
void identity_matrix(double ** m, int n);
void transpose(double *matrix, double ** transposed_matrix, int n);
void multiply(double *m1, double *m2, double * new_matrix, int n);
void print_matrix(double *m, int rows, int columns);
double* gen_matrix(int n);
void swap_rows(double * i_mat, int row1, int row2, int n);

// Jacobi method functions
void jacobi(double * a, int n, double * s, double * u, double * v);
void generate_composite_matrix(double ** matrix, int p, int q, int n, double angle, int side);
void rotate(double * a, int n, double * u, double * v, int p, int q);
int not_converged(double * a, int n);
void order_a(double * a, int n, double * u, double * v);
void generate_s(double * a, double * s, int n);

void jacobi(double * a, int n, double * s, double * u, double * v) {
    // Creates a copy of A, make all changes to A in A
    double * a_copy = (double *)malloc(sizeof(double) * n * n);
    copy_matrix(a, a_copy, n);
    
    // Initializes U and V as identity matrices
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                u[i*n + j] = 1.0;
                v[i*n + j] = 1.0;
            } else {
                u[i*n + j] = 0.0;
                v[i*n + j] = 0.0;
            }
        }
    }
    int sweep = 0;
    int max_sweeps = INT_MAX;

    while (not_converged(a, n)) { // epsilon is 10^-15
        for (int p = 0; p < n-1; p++) {
            for (int q = p+1; q < n; q ++) {
                rotate(a, n, u, v, p, q);
            }
        }
        sweep++;
        if (sweep > max_sweeps) ERROR("failing to converge");
    }

    order_a(a, n, u, v);
    generate_s(a, s, n);
    // TODO: order a, u, and v. s must have non-neg decreasing values
}

void rotate(double * a, int n, double * u, double * v, int p, int q) {
    double theta, phi;
    double a_pp = a[p*n + p];
    double a_pq = a[p*n + q];
    double a_qp = a[q*n + p];
    double a_qq = a[q*n + q];
    // calculate left and right rotation angles, see Brent, et. al
    double theta_phi = atan((a_qp + a_pq) / (a_qq - a_pp));
    double theta_phi_d = atan((a_qp - a_pq) / (a_qq + a_pp));
    theta = (theta_phi - theta_phi_d) / 2;
    phi = (theta_phi + theta_phi_d) / 2;
    double * l_matrix, * r_matrix;
    generate_composite_matrix(&l_matrix, p, q, n, theta, LEFT);
    generate_composite_matrix(&r_matrix, p, q, n, phi, RIGHT);
    
    double * temp = (double *)malloc(sizeof(double) * n * n);
    // Rotate a to the left and right
    multiply(l_matrix, a, temp, n);
    multiply(temp, r_matrix, a, n);
    
    // Apply rotation on U
    multiply(l_matrix, u, temp, n);
    copy_matrix(temp, u, n);
    
    // Apply rotation on V

    multiply(v, r_matrix, temp, n);
    copy_matrix(temp, v, n);
    free(temp); free(l_matrix); free(r_matrix);
}

void order_a(double * a, int n, double * u, double * v) {
    int d_i, sort_index; // d_i is diagonal element index, sort_index is index of sort
    double * permutation_matrix;
    double * temp_matrix = (double *)malloc(sizeof(double) * n * n);
    identity_matrix(&permutation_matrix, n);

    for (sort_index = 0; sort_index < n; sort_index++) {
        double max_d = 0.0;   // maximum diagonalized element
        int max_d_index = 0;
        for (d_i = sort_index; d_i < n; d_i++) {
            if (fabs(a[d_i*n + d_i]) > fabs(max_d)) {
                max_d = a[d_i*n + d_i];
                max_d_index = d_i;
            }
        }
        if (max_d_index > 0) {
            // reordering values in a
            double temp_d = a[sort_index * n + sort_index];
            a[sort_index * n + sort_index] = max_d;
            a[max_d_index * n + max_d_index] = temp_d;
            swap_rows(permutation_matrix, sort_index, max_d_index, n);
            // reorder u and v wrt swapped values in a
            multiply(permutation_matrix, u, temp_matrix, n);
            copy_matrix(temp_matrix, u, n);
            multiply(v, permutation_matrix, temp_matrix, n);
            copy_matrix(temp_matrix, v, n);
            // set permutation matrix back to identity
            free(permutation_matrix);
            identity_matrix(&permutation_matrix, n);
        }
    }
    for (d_i = 0; d_i < n; d_i ++) {
        // if diagonal element less than 0, make it positive
        if (a[d_i*n + d_i] < 0.0) {
            permutation_matrix[d_i*n + d_i] = -1.0;
        }
    }
    multiply(permutation_matrix, a, temp_matrix, n);
    copy_matrix(temp_matrix, a, n);
    multiply(u, permutation_matrix, temp_matrix, n);
    copy_matrix(temp_matrix, u, n);
    free(permutation_matrix); free(temp_matrix);
}

// swapping rows of an identity matrix
void swap_rows(double * i_mat, int row1, int row2, int n) {
    i_mat[row1*n + row1] = 0.0;
    i_mat[row1*n + row2] = 1.0;
    i_mat[row2*n + row2] = 0.0;
    i_mat[row2*n + row1] = 1.0;
}

void generate_s(double * a, double * s, int n) {
    for (int d = 0; d < n; d ++) {
        s[d] = a[d*n + d];
    }
}

// TODO: optimize, just multiply col p and q, then row p and q
void multiply(double * m1, double * m2, double * new_matrix, int n) {
    int new_row, new_column, old_x;

    for (new_row = 0; new_row < n; new_row++) {
        
        for (new_column = 0; new_column < n; new_column++) {
            new_matrix[new_row*n + new_column] = 0;
            
            for (old_x = 0; old_x < n; old_x++) {
                new_matrix[new_row*n + new_column] += m1[new_row*n + old_x] \
                    * m2[old_x*n + new_column];
            }
        }
    }
}

void generate_composite_matrix(double ** matrix, int p, int q, int n, double angle, int side) {
    identity_matrix(matrix, n);
    (*matrix)[n*p + p] = (*matrix)[n*q + q] = cos(angle);
    if (side == LEFT) {
        (*matrix)[n*p + q] = -sin(angle);
        (*matrix)[n*q + p] = sin(angle);
    }
    else {
        (*matrix)[n*p + q] = sin(angle);
        (*matrix)[n*q + p] = -sin(angle);
    }
}

int not_converged(double *a, int n) {
    double epsilon = 1e-15;
    int row, col;
    for(row = 0; row < n; row ++) {
        for(col = 0; col < n; col ++) {
            if (row == col) continue;
            else {
                if (fabs(a[row*n + col]) > epsilon) // not yet converged, error > epsilon
                    return 1; 
            }
        }
    }
    return 0;
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

void transpose(double * a, double ** transposed_matrix, int n) {
    int row, column;
    *transposed_matrix = (double*)malloc(sizeof(double)*n*n);
    
    for (row=0; row < n; row++) {
        for (column=0; column < n; column++) {
            (*transposed_matrix)[row*n+column] = a[column*n+row]; 
            (*transposed_matrix)[column*n+row] = a[row*n+column];
        }
    }
}

double* gen_matrix(int n) {
    double* a = malloc(n * n * sizeof(double));
    int i, j;
    double i_p_1, j_p_1;
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            i_p_1 = ((double) i) + 1;
            j_p_1 = ((double) j) + 1;
            a[i * n + j] = sqrt(i_p_1 * i_p_1 + j_p_1 * j_p_1);
        }
    }
    return a;
}

void print_matrix(double *m, int rows, int columns) {
    int row, column;
    for (row=0; row<rows; row++) {
        for (column=0; column<columns; column++) {
            printf("[%d, %d]: %f\t", row, column, m[row*rows+column]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char const *argv[]) {
    if (argc != 2) {
        ERROR("Invoke as ./hw2 [n] \n");
    }
    int n = atoi(argv[1]);
    double *a, *u, *v, *s;
    u = (double *)malloc(sizeof(double) * n * n);
    v = (double *)malloc(sizeof(double) * n * n);
    s = malloc(sizeof(double) * n);
    a = gen_matrix(n);
    jacobi(a, n, s, u, v);

    printf("\nFINAL a, s, u, v\n");
    
    printf("\na = \n");
    print_matrix(a, n, n);
    
    printf("\ns = \n");
    print_matrix(s, 1, n);

    printf("\nu = \n");
    print_matrix(u, n, n);
    
    printf("\nv = \n");
    print_matrix(v, n, n);

    free(s); free(u); free(v); free(a);
    return 0;
}