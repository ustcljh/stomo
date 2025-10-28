#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

int mat_alloc_counter = 0;

matrix_cap* _matcap(matrix mat)
{
	return (matrix_cap*)((uint8_t*)mat - sizeof(matrix_cap));
}

matrix _matdat(matrix_cap* cap)
{
	return (matrix)((uint8_t*)cap + sizeof(matrix_cap));
}

matrix mat_alloc(int size1, int size2)
{
	if (size1 <= 0)
	{
		fprintf(stderr, "error: matalloc cannot allocate matrix with invalid size1 = %d\n", size1);
		abort();
	}

	if (size2 <= 0)
	{
		fprintf(stderr, "error: matalloc cannot allocate matrix with invalid size2 = %d\n", size2);
		abort();
	}

	/* The allocated memory is like
	*
	* pcap____      pdata__
	*        |            |
	*  ______v____________v_________________________
	* |  cap   ...  |  double*  |  double*  |  ...
	* |             |    [0]    |    [1]    |  ...
	* |_____________|___________|___________|_______
	*    (12 byte)        |           |_____________________________
	* (possible 16b)      |                                         |
	*                _____v_________________________           _____v_________________________
	*               |   double  |   double  |  ...            |   double  |   double  |  ...
	*               |   [0][0]  |   [0][1]  |  ...            |   [1][0]  |   [1][1]  |  ...
	*               |___________|___________|_______          |_______________________________
	*/

	// Compute the first dimension size
	int bsz1 = size1 * sizeof(double*) + sizeof(matrix_cap);

	// Allocate the first part (dimension 1 index + cap)
	matrix_cap* pcap = malloc(bsz1);
	if (pcap == NULL)
	{
		fprintf(stderr, "error: matalloc failed to allocate %d bytes of memory for part I matrix.\n", bsz1);
		abort();
	}

	pcap->size1 = size1;
	pcap->size2 = size2;

	// Compute where the data region really lies, skip the cap part
	matrix pdata = _matdat(pcap);

	// Allocate each sub array
	for (int i = 0; i < size1; ++i)
	{
		int bsz2 = size2 * sizeof(double);

		// Do the allocation
		double* psubdata = malloc(bsz2);
		if (psubdata == NULL)
		{
			fprintf(stderr, "error: matalloc failed to allocate %d bytes of memory for part II sub array in matrix.\n", bsz2);
			abort();
		}

		// Clean to zero for the newly allocated part
		memset(psubdata, 0, bsz2);

		// Save the allocated address
		pdata[i] = psubdata;
	}

	return pdata;
}

void mat_free(matrix mat)
{
	if (mat == NULL)
	{
		fprintf(stderr, "error: matfree got NULL matrix.\n");
		abort();
	}

	for (int i = 0; i < mat_size1(mat); ++i)
	{
		if (mat[i] != NULL)
		{
			free(mat[i]);
		}
		else
		{
			fprintf(stderr, "error: matfree found NULL dimension 2 at [%d] = %p.\n", i, mat[i]);
			abort();
		}
	}

	free(_matcap(mat));
}

int mat_size1(matrix mat)
{
	return _matcap(mat)->size1;
}

int mat_size2(matrix mat)
{
	return _matcap(mat)->size2;
}

matrix mat_resize1(matrix mat, int new_size)
{
	if (mat == NULL)
	{
		fprintf(stderr, "error: matresize1 got NULL matrix.\n");
		abort();
	}

	// Use current size if new size is -1
	if (new_size == -1)
	{
		new_size = mat_size1(mat);
	}

	if (new_size <= 0)
	{
		fprintf(stderr, "error: matresize1 cannot allocate matrix with invalid size = %d\n", new_size);
		abort();
	}

	matrix_cap* pcap = _matcap(mat);

	// Check and resize dimension 1
	if (new_size != mat_size1(mat))
	{
		// Compute dimension 1 new size in byte
		int bsz = new_size * sizeof(double*) + sizeof(matrix_cap);

		// If need to shrink, free the corresponding dimension II arrays
		if (new_size < mat_size1(mat))
		{
			for (int i = new_size; i < mat_size1(mat); ++i)
			{
				if (mat[i] != NULL)
				{
					free(mat[i]);
				}
				else
				{
					fprintf(stderr, "error: matresize1 found NULL dimension 2 at [%d] = %p.\n", i, mat[i]);
					abort();
				}
			}
		}

		// Do the reallocation
		matrix_cap* new_pcap = realloc(pcap, bsz);

		if (new_pcap == NULL)
		{
			fprintf(stderr, "error: matresize1 failed to allocate %d bytes of memory for part I matrix.\n", bsz);
			abort();
		}

		// Update the pcap pointer, and mat associated to it
		pcap = new_pcap;
		mat = _matdat(pcap);

		// If need to expand, allocate dimension 2 arrays for the newly expanded region
		if (new_size > mat_size1(mat))
		{
			for (int i = mat_size1(mat); i < new_size; ++i)
			{
				int bsz2 = mat_size2(mat) * sizeof(double);

				// Do the allocation
				double* psubdata = malloc(bsz2);
				if (psubdata == NULL)
				{
					fprintf(stderr, "error: matresize1 failed to allocate %d bytes of memory for part II sub array in matrix.\n", bsz2);
					abort();
				}

				// Clean to zero for the newly allocated part
				memset(psubdata, 0, bsz2);

				// Save the allocated address
				mat[i] = psubdata;
			}
		}

		// Update the size1 metadata
		pcap->size1 = new_size;
	}

	return mat;
}

matrix mat_resize2(matrix mat, int new_size)
{
	if (mat == NULL)
	{
		fprintf(stderr, "error: matresize2 got NULL matrix.\n");
		abort();
	}

	// Use current size if new size is -1
	if (new_size == -1)
	{
		new_size = mat_size2(mat);
	}

	if (new_size <= 0)
	{
		fprintf(stderr, "error: matresize2 cannot allocate matrix with invalid size = %d\n", new_size);
		abort();
	}

	matrix_cap* pcap = _matcap(mat);

	// Check and resize dimension 2
	if (new_size != mat_size2(mat))
	{
		for (int i = 0; i < mat_size1(mat); ++i)
		{
			// Compute dimension II size
			int bsz = new_size * sizeof(double);

			if (mat[i] == NULL)
			{
				fprintf(stderr, "error: matresize2 got NULL dimension 2 array.\n");
				abort();
			}

			// Do the reallocation
			double* new_subarr = realloc(mat[i], bsz);

			if (new_subarr == NULL)
			{
				fprintf(stderr, "error: matresize2 failed to allocate %d bytes of memory for part I matrix.\n", bsz);
				abort();
			}

			// Update the subarray pointer
			mat[i] = new_subarr;

			// Clean to zero for the expanded part
			int sizediff = new_size - mat_size2(mat);
			if (sizediff > 0)
			{
				memset(mat[i] + mat_size2(mat), 0, sizediff * sizeof(double));
			}
		}

		// Update size2 metadata
		pcap->size2 = new_size;
	}

	return mat;
}

matrix mat_matmul(matrix a, matrix b)
{
	if (mat_size2(a) != mat_size1(b))
	{
		fprintf(stderr, "mat_matmul: different size of dim2 (%d) in A and dim1 (%d) in B.\n", mat_size2(a), mat_size1(b));
		abort();
	}

	matrix s = mat_alloc(mat_size1(a), mat_size2(b));

	for (int i = 0; i < mat_size1(a); ++i)
	{
		for (int j = 0; j < mat_size2(b); ++j)
		{
			double sum = 0;

			for (int k = 0; k < mat_size2(a); ++k)
			{
				sum += a[i][k] * b[k][j];
			}

			s[i][j] = sum;
		}
	}

	return s;
}

matrix mat_copy(matrix source)
{
	matrix dest = mat_alloc(mat_size1(source), mat_size2(source));
	for (int i = 0; i < mat_size1(source); ++i)
	{
		for (int j = 0; j < mat_size2(source); ++j)
		{
			dest[i][j] = source[i][j];
		}
	}

	return dest;
}

matrix mat_transpose(matrix mat)
{
	matrix result = mat_alloc(mat_size2(mat), mat_size1(mat));

	for (int i = 0; i < mat_size1(mat); ++i)
	{
		for (int j = 0; j < mat_size2(mat); ++j)
		{
			result[j][i] = mat[i][j];
		}
	}

	return result;
}

void mat_print(matrix mat)
{
	for (int i = 0; i < mat_size1(mat); ++i)
	{
		for (int j = 0; j < mat_size2(mat); ++j)
		{
			printf("%8g ", mat[i][j]);
		}
		printf("\n");
	}
}
