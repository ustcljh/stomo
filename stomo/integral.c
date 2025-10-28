#include "integral.h"
#include "global.h"
#include "structure.h"
#include "orbital.h"

#define _USE_MATH_DEFINES
#include <math.h>

#define orb(id, gto) (orbitals[(id)].orbital->gtos[(gto)])
#define atom(id) (orbitals[(id)].atom)

double integral_normalize_factor(double exponent)
{
	return pow(2.0 * exponent / M_PI, 0.75);
}

// Do the overlap integral on primitive gaussian (u.p, v.q)
double integral_primitive_s(int u, int p, int v, int q)
{
	return orb(u, p).coefficient *
		orb(v, q).coefficient *
		integral_normalize_factor(orb(u, p).exponent) *
		integral_normalize_factor(orb(v, q).exponent) *
		pow(M_PI / (orb(u, p).exponent + orb(v, q).exponent), 1.5) *
		exp(-((orb(u, p).exponent * orb(v, q).exponent) / (orb(u, p).exponent + orb(v, q).exponent)
			* pow(vector_dist(atom(u)->atom_pos, atom(v)->atom_pos), 2)));
}

double integral_s(int u, int v)
{
	double sum = 0;
	for (int p = 0; p < orbitals[u].orbital->ngtos; ++p)
	{
		for (int q = 0; q < orbitals[v].orbital->ngtos; ++q)
		{
			sum += integral_primitive_s(u, p, v, q);
		}
	}

	return sum;
}

// Do the kinetic integral on primitive gaussian (u.p, v.q)
double integral_primitive_t(int u, int p, int v, int q)
{
	// The (a * b) / (a + b) where a and b are the two exponents of u.p and v.q
	double ab_over_aplusb = (orb(u, p).exponent * orb(v, q).exponent) / (orb(u, p).exponent + orb(v, q).exponent);

	return integral_primitive_s(u, p, v, q) * ab_over_aplusb *
		(3 - 2 * ab_over_aplusb * pow(vector_dist(atom(u)->atom_pos, atom(v)->atom_pos), 2));
}

double integral_t(int u, int v)
{
	double sum = 0;
	for (int p = 0; p < orbitals[u].orbital->ngtos; ++p)
	{
		for (int q = 0; q < orbitals[v].orbital->ngtos; ++q)
		{
			sum += integral_primitive_t(u, p, v, q);
		}
	}

	return sum;
}

// The Boys function of order 0
double integral_fboys0(double t)
{
	if (t < 1e-10)
	{
		// F0(0) = 1, as the limit says. Do not go on since pi/t would probably overflow.
		return 1;
	}

	return 0.5 * sqrt(M_PI / t) * erf(sqrt(t));
}

// Do the potential integral with primitive gaussian (u.p, v.q) and core c
double integral_primitive_v(int u, int p, int v, int q, int c)
{
	// Two exponents
	double a = orb(u, p).exponent;
	double b = orb(v, q).exponent;

	// The core c
	atom* cc = &(strut.atoms[c]);

	// Exponent weight-ed point P of the two cores atom(u) and atom(v)
	vector pp = vector_add(vector_linmul(a / (a + b), atom(u)->atom_pos),
		vector_linmul(b / (a + b), atom(v)->atom_pos));

	return -(cc->atom_number) *
		orb(u, p).coefficient *
		orb(v, q).coefficient *
		integral_normalize_factor(a) *
		integral_normalize_factor(b) *
		(2 * M_PI / (a + b)) * exp(-a * b / (a + b) * pow(vector_dist(atom(u)->atom_pos, atom(v)->atom_pos), 2)) *
		integral_fboys0((a + b) * pow(vector_dist(pp, cc->atom_pos), 2));
}

double integral_v(int u, int v)
{
	double sum = 0;

	// Sum over all cores
	for (int core = 0; core < strut.natoms; ++core)
	{
		// Sum over all primitive gaussians
		for (int p = 0; p < orbitals[u].orbital->ngtos; ++p)
		{
			for (int q = 0; q < orbitals[v].orbital->ngtos; ++q)
			{
				sum += integral_primitive_v(u, p, v, q, core);
			}
		}
	}

	return sum;
}

// Do the four-center integrals, (pq|rs) over four primitive gaussians
double integral_four_primitive(int c1, int p, int c2, int q, int c3, int r, int c4, int s)
{
	// Four exponents
	double a = orb(c1, p).exponent;
	double b = orb(c2, q).exponent;
	double c = orb(c3, r).exponent;
	double d = orb(c4, s).exponent;

	// Exponent weighted centers P and Q
	vector pp = vector_add(vector_linmul(a / (a + b), atom(c1)->atom_pos),
		vector_linmul(b / (a + b), atom(c2)->atom_pos));
	vector qq = vector_add(vector_linmul(c / (c + d), atom(c3)->atom_pos),
		vector_linmul(d / (c + d), atom(c4)->atom_pos));

	// Distance of |PQ|, |AB| and |CD|
	double pq = vector_dist(pp, qq);
	double ab = vector_dist(atom(c1)->atom_pos, atom(c2)->atom_pos);
	double cd = vector_dist(atom(c3)->atom_pos, atom(c4)->atom_pos);

	// Factor for pq and rs
	double facpq = exp(-a * b / (a + b) * pow(ab, 2));
	double facrs = exp(-c * d / (c + d) * pow(cd, 2));

	return orb(c1, p).coefficient * orb(c2, q).coefficient *
		orb(c3, r).coefficient * orb(c4, s).coefficient *
		(2 * pow(M_PI, 2.5) / ((a + b) * (c + d) * sqrt(a + b + c + d))) *
		facpq * facrs *
		integral_fboys0(((a + b) * (c + d) / (a + b + c + d)) * pow(pq, 2));
}

double integral_four(int c1, int c2, int c3, int c4)
{
	double sum = 0;

	for (int p = 0; p < orbitals[c1].orbital->ngtos; ++p)
	{
		for (int q = 0; q < orbitals[c2].orbital->ngtos; ++q)
		{
			for (int r = 0; r < orbitals[c3].orbital->ngtos; ++r)
			{
				for (int s = 0; s < orbitals[c4].orbital->ngtos; ++s)
				{
					sum += integral_four_primitive(c1, p, c2, q, c3, r, c4, s);
				}
			}
		}
	}

	return sum;
}
