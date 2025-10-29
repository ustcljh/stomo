#include <stdio.h>

#include "matrix.h"
#include "basis.h"
#include "structure.h"
#include "orbital.h"
#include "integral.h"
#include "linalg.h"
#include "scf.h"

int main(void)
{
	basis_load("D:\\dev\\stomo\\data\\basis\\bse.txt");
	// basis_print();
	structure_load("D:\\dev\\stomo\\data\\struct.txt");
	orbitals_init();
	// orbitals_print();

	/* printf("---- S matrix ----\n");
	matrix sint = mat_alloc(norbitals, norbitals);

	for (int i = 0; i < norbitals; ++i)
	{
		for (int j = 0; j < norbitals; ++j)
		{
			sint[i][j] = integral_s(i, j);
		}
	}

	mat_print(sint);

	printf("---- T matrix ----\n");
	matrix tint = mat_alloc(norbitals, norbitals);

	for (int i = 0; i < norbitals; ++i)
	{
		for (int j = 0; j < norbitals; ++j)
		{
			tint[i][j] = integral_t(i, j);
		}
	}

	mat_print(tint);

	printf("---- V matrix ----\n");
	matrix vint = mat_alloc(norbitals, norbitals);

	for (int i = 0; i < norbitals; ++i)
	{
		for (int j = 0; j < norbitals; ++j)
		{
			vint[i][j] = integral_v(i, j);
		}
	}

	mat_print(vint);

	printf("---- mat matrix ----\n");
	matrix mat = mat_alloc(3, 3);
	mat[0][0] = 1;
	mat[0][1] = 2;
	mat[0][2] = 3;
	mat[1][0] = 2;
	mat[1][1] = 5;
	mat[1][2] = 6;
	mat[2][0] = 3;
	mat[2][1] = 6;
	mat[2][2] = 10;
	mat_print(mat);

	printf("---- mat^(-1/2) matrix ----\n");
	matrix invsqrt = linalg_inv_sqrt(mat);
	mat_print(invsqrt);

	matrix invsqrt_T = mat_transpose(invsqrt);
	matrix xm = mat_matmul(invsqrt_T, mat);
	matrix xmx = mat_matmul(xm, invsqrt);

	mat_print(xmx); */

	scf(1);

	/* matrix m = mat_alloc(5, 5);
	mat_print(m);

	printf("\n");

	m = mat_resize1(m, 6);
	m[5][1] = 1958;

	mat_print(m);

	printf("\n");

	m = mat_resize2(m, 8);
	m[5][7] = 2025.12345678;

	mat_print(m);

	printf("\n");

	m = mat_resize1(m, 2);
	m = mat_resize2(m, 3);

	mat_print(m);

	printf("\n");

	mat_free(m); */
}
