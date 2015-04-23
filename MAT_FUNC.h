#include "MAT_FUNC.c"

// matrix functions
void copy_matrix(double * a, double * b, int n);
void identity_matrix(double ** m, int n);
void transpose(double *matrix, int n);
void multiply(double *m1, double *m2, double * new_matrix, int n);
void print_matrix(double *m, int rows, int columns);
double* gen_matrix(int n);
void swap_rows(double * i_mat, int row1, int row2, int n);
