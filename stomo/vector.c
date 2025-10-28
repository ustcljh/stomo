#include "vector.h"

#include <math.h>

vector vector_add(vector a, vector b)
{
	vector result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;

	return result;
}

vector vector_sub(vector a, vector b)
{
	vector result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;

	return result;
}

vector vector_linmul(double c, vector a)
{
	vector result;
	result.x = a.x * c;
	result.y = a.y * c;
	result.z = a.z * c;

	return result;
}

double vector_len(vector a)
{
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double vector_dist(vector a, vector b)
{
	return vector_len(vector_sub(a, b));
}
