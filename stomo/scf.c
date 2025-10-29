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

/* Call this right after you detect convergence (before freeing F, C, e, Cp, etc.) */

void debug_print_orbitals(matrix F, matrix C, matrix S, matrix Fp, double* e, matrix Cp, int norb)
{
	printf("\n--- DEBUG ORBITAL CHECKS ---\n");

	/* 1) Print eigenvalues returned from orthonormal eigenproblem (these are canonical) */
	printf("Eigenvalues from Fp diagonalization (e[i]):\n");
	for (int i = 0; i < norb; ++i) printf(" e[%d] = %.12f\n", i, e[i]);

	/* 2) Rayleigh quotient: (C_i^T F C_i) / (C_i^T S C_i) */
	printf("\nRayleigh quotients (C_i^T F C_i) / (C_i^T S C_i):\n");
	for (int i = 0; i < norb; ++i)
	{
		double num = 0.0, den = 0.0;
		for (int u = 0; u < norb; ++u)
		{
			for (int v = 0; v < norb; ++v)
			{
				num += C[u][i] * F[u][v] * C[v][i];
				den += C[u][i] * S[u][v] * C[v][i];
			}
		}
		printf(" Rayleigh[%d] = %.12f   (denominator = %.12f)\n", i, num / den, den);
	}

	/* 3) Diagonal of C^T F C and diagonal of C^T S C */
	matrix Ct = mat_transpose(C);
	matrix CtF = mat_matmul(Ct, F);
	matrix CtFC = mat_matmul(CtF, C);
	matrix CtS = mat_matmul(Ct, S);
	matrix CtSC = mat_matmul(CtS, C);

	printf("\nDiagonal of C^T F C:\n");
	for (int i = 0; i < norb; ++i) printf(" CtFC[%d,%d] = %.12f\n", i, i, CtFC[i][i]);

	printf("\nMatrix C^T S C (should be identity):\n");
	mat_print(CtSC);

	/* 4) Compare eigenvalues e[i] to Rayleigh and CtFC diagonals */
	printf("\nDifferences (e[i] - Rayleigh) and (e[i] - CtFC_diag):\n");
	for (int i = 0; i < norb; ++i)
	{
		double ray = 0.0, ctfcdiag = 0.0;
		double den = 0.0;
		for (int u = 0; u < norb; ++u)
			for (int v = 0; v < norb; ++v)
			{
				ray += C[u][i] * F[u][v] * C[v][i];
				den += C[u][i] * S[u][v] * C[v][i];
			}
		ray /= den;
		ctfcdiag = CtFC[i][i];
		printf(" i=%d: e=%.12f  ray=%.12f  ctfc=%.12f   e-ray=%.5e  e-ctfc=%.5e\n",
			i, e[i], ray, ctfcdiag, e[i] - ray, e[i] - ctfcdiag);
	}

	/* free temporaries */
	mat_free(Ct);
	mat_free(CtF);
	mat_free(CtFC);
	mat_free(CtS);
	mat_free(CtSC);

	printf("Cp:\n");
	mat_print(Cp);
	printf("Fp:\n");
	mat_print(Fp);

	printf("--- END DEBUG ---\n\n");
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

	/* printf("--- S matrix ---\n");
	mat_print(S); */

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

	printf("--- T matrix ---\n");
	mat_print(T);
	printf("--- V matrix ---\n");
	mat_print(V);

	matrix Hcore = mat_alloc(norbitals, norbitals);
	for (int u = 0; u < norbitals; ++u)
	{
		for (int v = 0; v < norbitals; ++v)
		{
			Hcore[u][v] = T[u][v] + V[u][v];
		}
	}

	printf("--- Hcore matrix ---\n");
	mat_print(Hcore);

	// Do the initial guess
	// Build up the initial Fock matrix Fp = Xt * Hcore * X
	matrix Fp = mat_mul_3(Xt, Hcore, X);

	// Solve eigenproblem F * Cp = Cp * e
	double* e;
	matrix Cp;
	linalg_eigen_decomp(Fp, &e, &Cp);

	/* printf("--- C' matrix ---\n");
	mat_print(Cp); */

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

	/* for (int i = 0; i < norbitals; ++i)
	{
		printf("Energy level %d: %.8lf\n", i, e[i]);
	} */

	double E_old = 0;

	// This is the main SCF loop
	for (int iter = 0; iter < MAX_SCF_ITER; ++iter)
	{
		printf("--- C matrix ---\n");
		mat_print(C);
		printf("--- C' matrix ---\n");
		mat_print(Cp);
		printf("--- P matrix ---\n");
		mat_print(P);

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

		printf("--- F matrix ---\n");
		mat_print(F);

		// Release old Fp, e, Cp, C
		/* mat_free(Fp);
		free(e);
		mat_free(Cp);
		mat_free(C); */

		// Transform Fp to orthonormal basis
		Fp = mat_mul_3(Xt, F, X);

		// Solve the eigenproblem
		linalg_eigen_decomp(Fp, &e, &Cp);

		printf("--- Fp matrix ---\n");
		mat_print(Fp);

		// Back transform C' to C
		C = mat_matmul(X, Cp);

		// Free the old P
		// mat_free(P);

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

		// Check P * S
		/* matrix PS = mat_matmul(P, S);
		printf("--- Matrix P*S ---\n");
		mat_print(PS); */

		// Ct S C
		/* matrix CtSC = mat_mul_3(mat_transpose(C), S, C);
		printf("--- Matrix Ct*S*C ---\n");
		mat_print(CtSC); */

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

		debug_print_orbitals(F, C, S, Fp, e, Cp, norbitals);

		if (fabs(E_old - Eelec) <= 1e-10)
		{
			printf("SCF converged.\n");

			printf("--- The MO coefficients: ---\n");
			mat_print(C);

			matrix temp = mat_matmul(mat_transpose(C), F);
			matrix eps_mat = mat_matmul(temp, C);
			// mat_free(temp);

			mat_print(eps_mat);

			for (int i = 0; i < norbitals; ++i) {
				printf("Orbital %d energy: %.10f\n", i, eps_mat[i][i]);
			}

			break;
		}

		E_old = Eelec;

		// mat_free(F);
	}
}
