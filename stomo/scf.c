#include "scf.h"

#include "orbital.h"
#include "matrix.h"
#include "linalg.h"
#include "integral.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define MAX_SCF_ITER 16

matrix mat_mul_3(matrix a, matrix b, matrix c)
{
	matrix temp = mat_matmul(a, b);
	matrix result = mat_matmul(temp, c);
	mat_free(temp);

	return result;
}

void scf(int nocc)
{
	// Compute the overlap integral matrix S
	matrix S = mat_alloc(norbitals, norbitals);
	for (int u = 0; u < norbitals; ++u)
	{
		for (int v = 0; v < norbitals; ++v)
		{
			S[u][v] = integral_s(u, v);
		}
	}

	// Compute the orthogonalizer X = S^(-1/2) and its transposed form Xt
	matrix X = linalg_inv_sqrt(S);
	matrix Xt = mat_transpose(X);

	// Compute Hcore
	matrix Hcore = mat_alloc(norbitals, norbitals);
	for (int u = 0; u < norbitals; ++u)
	{
		for (int v = 0; v < norbitals; ++v)
		{
			Hcore[u][v] = integral_t(u, v) + integral_v(u, v);
		}
	}

	// Do the initial guess
	// Build up the initial Fock matrix Fp = Xt * Hcore * X
	matrix Fp = mat_mul_3(Xt, Hcore, X);

	// Solve eigenproblem F * Cp = Cp * e
	double* e;
	matrix Cp;
	linalg_eigen_decomp(Fp, &e, &Cp);

	// Give back to C by C = X * Cp
	matrix C = mat_matmul(X, Cp);

	// Build the initial density matrix P with Puv = sum(2 * Cui * Cvi).
	matrix P = mat_alloc(norbitals, norbitals);
	for (int u = 0; u < norbitals; ++u)
	{
		for (int v = 0; v < norbitals; ++v)
		{
			P[u][v] = 0;
			for (int i = 0; i < nocc; ++i)
			{
				P[u][v] += 2 * C[u][i] * C[v][i];
				printf("Puv (%.3lf) += 2 * %.3lf * %.3lf\n", P[u][v], C[u][i], C[v][i]);
			}
		}
	}

	printf("--- C matrix ---\n");
	mat_print(C);
	printf("--- P matrix ---\n");
	mat_print(P);

	double E_old = 0;

	// This is the main SCF loop
	for (int iter = 0; iter < MAX_SCF_ITER; ++iter)
	{
		printf("--- C matrix ---\n");
		mat_print(C);
		printf("--- P matrix ---\n");
		mat_print(P);

		// Build the Fock matrix
		matrix F = mat_alloc(norbitals, norbitals);

		for (int u = 0; u < norbitals; ++u)
		{
			for (int v = 0; v < norbitals; ++v)
			{
				// Compute coulomb + exchange contrib
				double Guv = 0;
				for (int a = 0; a < norbitals; ++a)
				{
					for (int b = 0; b < norbitals; ++b)
					{
						Guv += P[a][b] * (integral_four(u, v, a, b) - 0.5 * integral_four(u, a, v, b));
					}
				}

				// printf("G%d%d = %.8lf\n", u, v, Guv);

				F[u][v] = Hcore[u][v] + Guv;
			}
		}

		printf("--- F matrix ---\n");
		mat_print(F);

		// Release old Fp, e, Cp, C
		mat_free(Fp);
		free(e);
		mat_free(Cp);
		mat_free(C);

		// Transform Fp to orthonormal basis
		Fp = mat_mul_3(Xt, F, X);

		// Solve the eigenproblem
		linalg_eigen_decomp(Fp, &e, &Cp);

		// Back transform C' to C
		C = mat_matmul(X, Cp);

		// Free the old P
		mat_free(P);

		// Form the new density matrix P
		P = mat_alloc(norbitals, norbitals);
		for (int u = 0; u < norbitals; ++u)
		{
			for (int v = 0; v < norbitals; ++v)
			{
				P[u][v] = 0;
				for (int i = 0; i < nocc; ++i)
				{
					P[u][v] += 2 * C[u][i] * C[v][i];
				}
			}
		}

		// Compute the electron energy
		double Eelec = 0;
		for (int u = 0; u < norbitals; ++u)
		{
			for (int v = 0; v < norbitals; ++v)
			{
				Eelec += 0.5 * P[u][v] * (Hcore[u][v] + F[u][v]);
			}
		}

		printf("SCF iter %d, Energy %.8lf\n", iter, Eelec);

		if (fabs(E_old - Eelec) <= 1e-6)
		{
			printf("SCF converged.\n");
			break;
		}

		E_old = Eelec;

		mat_free(F);
	}
}
