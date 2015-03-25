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
#include <stdint.h>
#include <math.h>

#define LEFT 0
#define RIGHT 1
// matrix functions
void copy_matrix(double * a, double ** b, int n);
void identity_matrix(double ** m, int n);
void transpose(double *matrix, double ** transposed_matrix, int n);
void multiply(double *m1, double *m2, double ** new_matrix, int p, int q, int n);

// Jacobi method functions
void jacobi(double * a, int n, double * s, double * u, double * v);
void compute_theta(double * a, double ** left_matrix);
void generate_composite_matrix(double ** matrix, int p, int q, double angle, int side);

void jacobi(double * a, int n, double * s, double * u, double * v) {
    // Creates a copy of A, make all changes to A in A
    double * a_copy;
    copy_matrix(a, &a_copy, n);
    
    // Initializes U and V as identity matrices
    identity_matrix(&v, n);
    identity_matrix(&u, n);
    int sweep = 0;

    double epsilon = 1e-15;

    while (not_converged(a)) {
        for (int p = 0; p < n-1; p++) {
            for (int q = q+1; q < n; q ++) {
                double theta, phi;
                double a_pp = a[p*n + p];
                double a_pq = a[p*n + q];
                double a_qp = a[q*n + p];
                double a_qq = a[q*n + q];
                double theta_phi = atan((a_qp + a_pq) / (a_qq - a_pp));
                double theta_phi_d = atan((a_qp - a_pq) / (a_qq + a_pp));
                theta = (theta_phi - theta_phi_d) / 2;
                phi = (theta_phi + theta_phi_d) / 2;
                double * l_matrix, * r_matrix;
                generate_composite_matrix(&l_matrix, p, q, theta, LEFT);
                generate_composite_matrix(&r_matrix, p, q, phi, RIGHT);
                double * temp;
                copy_matrix(a, &temp, n);
            }
        }
        sweep++;
    }
}

void multiply(double *m1, double *m2, double ** new_matrix, int p, int q, int n) {
    int new_row, new_column, old_x;

    for (new_row = 0; new_row < n; new_row++) {
        
        for (new_column = 0; new_column < n; new_column++) {
            (*new_matrix)[new_row*n + new_column] = 0;
            
            for (old_x = 0; old_x < n; old_x++) {
                (*new_matrix)[new_row*n + new_column] += m1[new_row*n + old_x] \
                    * m2[old_x*n + new_column];
            }
        }
    }
}

void generate_composite_matrix(double ** matrix, int p, int q, double angle, int side) {
    identity_matrix(matrix);
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

void copy_matrix(double * a, double ** b, int n) {
    *b = (double *)malloc(sizeof(double)*n*n);
    for (int i = 0; i < n*n; i++) {
        (*b)[i] = a[i];
    }
}

void identity_matrix(double **m, int n) {
    int i;
    *m = (double*)calloc(sizeof(double), n*n);
    for (i=0; i<n; i++) {
        (*m)[i*n+i] = (double)1;
    }
}

void transpose(double *a, double ** transposed_matrix, int n) {
    int row, column;
    *transposed_matrix = (double*)malloc(sizeof(double)*n*n);
    
    for (row=0; row < n; row++) {
        for (column=0; column < n; column++) {
            (*transposed_matrix)[row*n+column] = a[column*n+row]; 
            (*transposed_matrix)[column*n+row] = a[row*n+column];
        }
    }
}

int main () {

    return 0;
}