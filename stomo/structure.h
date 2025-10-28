#pragma once

#include "vector.h"

typedef struct
{
	vector atom_pos;
	int atom_number;
} atom;

typedef struct
{
	int natoms;
	atom* atoms;

	int nelec;
} structure;

extern structure strut;

void structure_load(const char* filename);
