#pragma once

typedef struct
{
	double x;
	double y;
	double z;
} vector;

// result = a + b;
vector vector_add(vector a, vector b);

// result = a - b;
vector vector_sub(vector a, vector b);

// result = c * a;
vector vector_linmul(double c, vector a);

// result = |a|
double vector_len(vector a);

// result = |a - b|
double vector_dist(vector a, vector b);
