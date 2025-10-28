#pragma once

// The real matrix
typedef double** matrix;

// The cap part before the data region of a matrix
typedef struct
{
	int size1;
	int size2;
	int alloc_count;
} matrix_cap;

// Matrix allocation counter
extern int mat_alloc_counter;

// Allocate a matrix with given size
matrix mat_alloc(int size1, int size2);

// Free an allocated matrix
void mat_free(matrix mat);

// Get the dimensional size
int mat_size1(matrix mat);
int mat_size2(matrix mat);

// Resize the matrix, use -1 to keep current size
matrix mat_resize1(matrix mat, int new_size);
matrix mat_resize2(matrix mat, int new_size);

// Matrix multiplication
matrix mat_matmul(matrix a, matrix b);

// Matrix copy
matrix mat_copy(matrix source);

// Matrix transpose
matrix mat_transpose(matrix mat);

// Pretty print the matrix
void mat_print(matrix mat);
