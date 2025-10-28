#pragma once

#include "matrix.h"

// Do the Jacobi rotation algorithm for computing the eigenvalues and eigenvectors of matrix mat
int linalg_eigen_decomp(matrix mat, double** eigenvalues, matrix* eigenvectors);

// Compute the X = mat^(-1/2) that have the property X^ mat X = I
matrix linalg_inv_sqrt(matrix mat);
