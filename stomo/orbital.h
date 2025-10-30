#pragma once

#include "structure.h"
#include "basis.h"

// Structure of one MO
typedef struct
{
	atom* atom;
	basis_orbital* orbital;
} orbital;

extern int norbitals;
extern orbital* orbitals;

// Initialize the orbitals from structure and basis
void orbitals_init();

// Pretty print the orbitals
void orbitals_print();

// Normalization factors
double s_normalize_factor(double exponent);
