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
	matrix T = mat_alloc(norbitals, norbitals);
	matrix V = mat_alloc(norbitals, norbitals);
	for (int u = 0; u < norbitals; ++u)
	{
		for (int v = 0; v < norbitals; ++v)
		{
			T[u][v] = integral_t(u, v);
			V[u][v] = integral_v(u, v);
		}
	}

	matrix Hcore = mat_alloc(norbitals, norbitals);
	for (int u = 0; u < norbitals; ++u)
	{
		for (int v = 0; v < norbitals; ++v)
		{
			Hcore[u][v] = T[u][v] + V[u][v];
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
			}
		}
	}

	double E_old = 0;

	// This is the main SCF loop
	for (int iter = 0; iter < MAX_SCF_ITER; ++iter)
	{
		// Build the Fock matrix
		matrix F = mat_alloc(norbitals, norbitals);

		for (int u = 0; u < norbitals; ++u)
		{
			for (int v = 0; v < norbitals; ++v)
			{
				double Guv = 0;
				for (int a = 0; a < norbitals; ++a)
				{
					for (int b = 0; b < norbitals; ++b)
					{
						double Juvab = integral_four(u, v, a, b);
						double Kuvab = integral_four(u, a, v, b);
						Guv += P[a][b] * (Juvab - 0.5 * Kuvab);
					}
				}
				F[u][v] = Hcore[u][v] + Guv;
			}
		}

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

		mat_free(F);

		if (fabs(E_old - Eelec) <= 1e-10)
		{
			printf("SCF converged.\n");

			printf("--- The MO coefficients: ---\n");
			mat_print(C);

			for (int i = 0; i < norbitals; ++i)
			{
				printf("Orbital %d energy: %.10f\n", i, e[i]);
			}

			break;
		}

		E_old = Eelec;
	}

	mat_free(S);
	mat_free(X);
	mat_free(Xt);
	mat_free(T);
	mat_free(V);
	mat_free(Hcore);
	mat_free(Fp);
	mat_free(Cp);
	mat_free(C);
	mat_free(P);

	free(e);
}
