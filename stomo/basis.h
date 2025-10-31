#pragma once

#define BASIS_ORBITAL_TYPE_S 0
#define BASIS_ORBITAL_TYPE_P 1

// Gaussian function term in basis
typedef struct
{
	double exponent;
	double coefficient;
} basis_gto;

// Orbital in basis
typedef struct
{
	int type;
	int shell;

	int amx, amy, amz;

	int ngtos;
	basis_gto* gtos;
} basis_orbital;

// Basis for an atom type
typedef struct
{
	int norbitals;
	basis_orbital* orbitals;
} basis_atom;

extern int basis_size;
extern basis_atom* basis;

// Load from basis file
void basis_load(const char* basis_filename);

// Pretty print the basis information
void basis_print();

// Clean loaded basis
void basis_clean();
