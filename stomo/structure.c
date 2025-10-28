#include "structure.h"

#include <stdio.h>
#include <stdlib.h>

structure strut;

void structure_load(const char* filename)
{
	// Open the structure file.
	FILE* file_structure = fopen(filename, "r");
	if (file_structure == NULL)
	{
		fprintf(stderr, "error: load_structure cannot open structure file: %s.\n", filename);
		abort();
	}

	// Read number of electrons
	if (fscanf(file_structure, "%d%d", &strut.natoms, &strut.nelec) != 2)
	{
		fprintf(stderr, "error: load_structure cannot read total atom number and/or total electron number.\n");
		abort();
	}

	// Allocate the required entries
	strut.atoms = malloc(strut.natoms * sizeof(atom));
	if (strut.atoms == NULL)
	{
		fprintf(stderr, "error: load_structure cannot allocate %d elements for atom.\n", strut.natoms);
		abort();
	}

	for (int i = 0; i < strut.natoms; ++i)
	{
		// Do the read in
		if (fscanf(file_structure, "%d%lf%lf%lf",
			&strut.atoms[i].atom_number,
			&strut.atoms[i].atom_pos.x,
			&strut.atoms[i].atom_pos.y,
			&strut.atoms[i].atom_pos.z) != 4)
		{
			fprintf(stderr, "error: load structure cannot read atom number or XYZ coordinate for atom %d.\n", i);
		}
	}
}
