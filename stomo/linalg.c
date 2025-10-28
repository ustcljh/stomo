#include "linalg.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 100
#define EPS 1e-10

int linalg_eigen_decomp(matrix mat, double** eigenvalues, matrix* eigenvectors)
{
	if (mat_size1(mat) != mat_size2(mat))
	{
		fprintf(stderr, "linalg_eigen_decomp: not a n*n matrix.\n");
		abort();
	}
	int n = mat_size1(mat);

	// Initialize results
	*eigenvalues = malloc(sizeof(double) * n);
	if (*eigenvalues == NULL)
	{
		fprintf(stderr, "linalg_eigen_decomp: failed to allocate %d units for eigenvalues.", n);
		abort();
	}

	*eigenvectors = mat_alloc(n, n);

	// Initial to In for eigenvectors
	for (int i = 0; i < n; ++i)
	{
		(*eigenvectors)[i][i] = 1;
	}

	// Temporary matrix for diagonalize
	matrix temp = mat_copy(mat);

	for (int iter = 0; iter < MAX_ITER; iter++) {
		// Find largest off-diagonal element
		int p = 0, q = 1;
		double max = fabs(temp[p][q]);
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double val = fabs(temp[i][j]);
				if (val > max)
				{
					max = val;
					p = i;
					q = j;
				}
			}
		}

		if (max < EPS)
			break; // converged

		double app = temp[p][p];
		double aqq = temp[q][q];
		double apq = temp[p][q];

		double phi = 0.5 * atan2(2.0 * apq, aqq - app);
		double c = cos(phi);
		double s = sin(phi);

		// Rotate matrix A
		for (int k = 0; k < n; k++) {
			double aik = temp[p][k];
			double aqk = temp[q][k];
			temp[p][k] = c * aik - s * aqk;
			temp[q][k] = s * aik + c * aqk;
		}
		for (int k = 0; k < n; k++) {
			double akp = temp[k][p];
			double akq = temp[k][q];
			temp[k][p] = c * akp - s * akq;
			temp[k][q] = s * akp + c * akq;
		}

		// Update diagonal
		temp[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
		temp[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
		temp[p][q] = 0.0;
		temp[q][p] = 0.0;

		// Update eigenvectors
		for (int k = 0; k < n; k++) {
			double vip = (*eigenvectors)[k][p];
			double viq = (*eigenvectors)[k][q];
			(*eigenvectors)[k][p] = c * vip - s * viq;
			(*eigenvectors)[k][q] = s * vip + c * viq;
		}
	}

	// Diagonal entries are eigenvalues
	for (int i = 0; i < n; i++)
		(*eigenvalues)[i] = temp[i][i];

	mat_free(temp);

	return n;
}

matrix linalg_inv_sqrt(matrix mat)
{
	if (mat_size1(mat) != mat_size2(mat))
	{
		fprintf(stderr, "linalg_inv_sqrt: not a n*n matrix.\n");
		abort();
	}
	int n = mat_size1(mat);

	// Do the diagonalize
	double* eigenval;
	matrix eigenvec;
	linalg_eigen_decomp(mat, &eigenval, &eigenvec);

	matrix dhalf = mat_alloc(n, n);
	for (int i = 0; i < n; i++)
	{
		if (eigenval[i] < 1e-12)
		{
			eigenval[i] = 1e-12;
		}
		dhalf[i][i] = 1.0 / sqrt(eigenval[i]);
	}

	matrix temp = mat_matmul(eigenvec, dhalf);
	matrix eigenvec_transpose = mat_transpose(eigenvec);
	matrix result = mat_matmul(temp, eigenvec_transpose);

	mat_free(eigenvec);
	mat_free(temp);
	mat_free(eigenvec_transpose);
	free(eigenval);

	return result;
}
