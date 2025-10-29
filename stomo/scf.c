#include "scf.h"

#include "orbital.h"
#include "matrix.h"
#include "linalg.h"
#include "integral.h"

#include <stdlib.h>

void scf()
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
	matrix F0 = mat_matmul(Xt, Hcore);
	matrix Fp = mat_matmul(F0, X);

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
			for (int i = 0; i < norbitals; ++i)
			{
				P[u][v] += 2 * C[u][i] * C[v][i];
			}
		}
	}
}
