#include "orbital.h"

#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

int norbitals;
orbital* orbitals;

double s_normalize_factor(double exponent)
{
	return pow(2.0 * exponent / M_PI, 0.75);
}

void orbitals_init()
{
	if (strut.natoms == 0)
	{
		fprintf(stderr, "error: orbitals_init found no atoms in structure.\n");
		abort();
	}

	norbitals = 0;

	// Allocate the orbitals
	orbitals = NULL;

	// Set orbital infos
	for (int i = 0; i < strut.natoms; ++i)
	{
		// The atom number
		int atom_number = strut.atoms[i].atom_number;

		// Orbital number for this atom
		int nAO = basis[atom_number].norbitals;

		// Reallocate the orbitals array
		orbital* new_orbitals = realloc(orbitals, (norbitals + nAO) * sizeof(orbital));
		if (new_orbitals == NULL)
		{
			fprintf(stderr, "error: orbitals_init cannot allocate %d MOs.\n", norbitals + nAO);
			abort();
		}
		orbitals = new_orbitals;

		// Read in the AOs
		for (int j = norbitals, orb = 0; j < norbitals + nAO; ++j, ++orb)
		{
			orbitals[j].atom = &strut.atoms[i];
			orbitals[j].orbital = &basis[atom_number].orbitals[orb];
		}

		norbitals += nAO;
	}
}

void orbitals_print()
{
	for (int i = 0; i < norbitals; ++i)
	{
		printf("Orbital %d: Atom %d, Center %.8lf %.8lf %.8lf, Shell %d, Type %d, %d GTOs\n",
			i,
			orbitals[i].atom->atom_number,
			orbitals[i].atom->atom_pos.x,
			orbitals[i].atom->atom_pos.y,
			orbitals[i].atom->atom_pos.z,
			orbitals[i].orbital->shell,
			orbitals[i].orbital->type,
			orbitals[i].orbital->ngtos);
	}
}
