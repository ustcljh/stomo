#include "atom.h"

#include <string.h>
#include <stdbool.h>

const char* atom_names[] = { "?", "H", "He", "Li", "Be", "B",
	"C", "N", "O", "F", "Ne",
	"Na", "Mg", "Al", "Si", "P",
	"S", "Cl", "Ar", "\0" };

int atom_lookup(const char* name)
{
	int atom_number = 0;

	while (true)
	{
		if (atom_names[atom_number][0] == '\0')
		{
			// End of atom name table: not found
			return -1;
		}

		if (strcmp(name, atom_names[atom_number]) == 0)
		{
			// Matched
			return atom_number;
		}

		atom_number++;
	}
}
