#include "plot.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "orbital.h"

double plot_wfunc_orb_value(orbital* orb, vector v)
{
	double val = 0.0;

	// printf("X %.3lf Y %.3lf Orb center %.3lf\n", v.x, v.y, orb->atom->atom_pos.x);

	if (orb->orbital->type == BASIS_ORBITAL_TYPE_S)
	{
		basis_gto* gtos = orb->orbital->gtos;

		for (int i = 0; i < orb->orbital->ngtos; ++i)
		{
			val += s_normalize_factor(gtos[i].exponent) * gtos[i].coefficient * exp(-gtos[i].exponent * pow(vector_dist(v, orb->atom->atom_pos), 2));
		}
	}

	return val;
}

double plot_wfunc_value(matrix C, int index, vector v)
{
	double sum = 0.0;
	for (int i = 0; i < norbitals; ++i)
	{
		sum += plot_wfunc_orb_value(&(orbitals[i]), v) * C[i][index];
	}

	return sum;
}

void plot_wfunc_gnuplot(const char* filename, matrix C, int index)
{
	FILE* f = fopen(filename, "w");

	if (f == NULL)
	{
		fprintf(stderr, "plot_wfunc_gnuplot: failed to open output file %s.\n", filename);
		abort();
	}

	int xstep = 100;
	int ystep = 100;
	double xfrom = -1.0;
	double xto = 3.0;
	double yfrom = -2.0;
	double yto = 2.0;
	double z = 0;

	for (int i = 0; i < xstep; ++i)
	{
		for (int j = 0; j < ystep; ++j)
		{
			double x = (xto - xfrom) / (xstep + 1) * i + xfrom;
			double y = (yto - yfrom) / (ystep + 1) * j + yfrom;

			vector v;
			v.x = x;
			v.y = y;
			v.z = z;

			fprintf(f, "%.8lf %.8lf %.8lf\n", x, y, plot_wfunc_value(C, index, v));
		}

		fprintf(f, "\n");
	}

	fclose(f);
}
