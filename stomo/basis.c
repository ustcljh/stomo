#include "basis.h"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "global.h"
#include "atom.h"

basis_atom* basis = NULL;
int basis_size = 0;

#define BASIS_LINE_SZBUF 1000
#define BASIS_ATOM_NAME_SZBUF 10
#define BASIS_ORBITAL_TYPE_SZBUF 10

int count_terms(const char* str)
{
	int count = 0;
	while (*str != '\0')
	{
		while (isspace(*str) && *str != '\0')
			++str;

		if (*str == '\0')
			break;

		++count;

		while (!isspace(*str) && *str != '\0')
			++str;
	}

	return count;
}

void basis_load(const char* basis_filename)
{
	// Open the basis file.
	FILE* basis_file = fopen(basis_filename, "r");

	if (basis_file == NULL)
	{
		fprintf(stderr, "error: basis_load cannot open basis file for read in.\n");
		abort();
	}

	// Thw current processing atom
	basis_atom* current_process_atom = NULL;
	int nshell = 0;

	while (true)
	{
		// Read in the line content into a buffer
		char linebuf[BASIS_LINE_SZBUF];
		fgets(linebuf, BASIS_LINE_SZBUF, basis_file);

		if (feof(basis_file))
		{
			break;
		}

		// Remove leading whitespaces
		char* plinebuf = linebuf;
		while (*plinebuf == ' ' ||
			*plinebuf == '\t')
		{
			++plinebuf;
		}

		// Skip the empty lines, comment lines, atom separator lines etc.
		if (*plinebuf == '!' ||
			*plinebuf == '\r' ||
			*plinebuf == '\n' ||
			*plinebuf == '\0')
		{
			continue;
		}

		if (*plinebuf == '*')
		{
			// There is **** after an atom basis. Reset the current record.
			current_process_atom = NULL;
			nshell = 0;
			continue;
		}

		// Check the line type
		int term_number = count_terms(plinebuf);
		if (term_number == 3)
		{
			// This is an orbital description line
			// like: S    3   1.00

			if (current_process_atom == NULL)
			{
				fprintf(stderr, "error: basis_load unexpected orbital specification line without atom: %s.\n", plinebuf);
				abort();
			}

			char orbital_type[BASIS_ORBITAL_TYPE_SZBUF];
			int number_of_gto;

			// Parse the input line
			if (sscanf(plinebuf, "%s%d", orbital_type, &number_of_gto) == 2)
			{
				int n_orbitals;

				if (strcmp(orbital_type, "S") == 0)
				{
					n_orbitals = 1;
				}
				else if (strcmp(orbital_type, "SP") == 0)
				{
					n_orbitals = 2;
				}
				else
				{
					fprintf(stderr, "error: basis_load unexpected orbital type: %s.\n", orbital_type);
					abort();
				}

				// Allocate memory for the orbitals
				basis_orbital* new_orbitals = realloc(current_process_atom->orbitals, (current_process_atom->norbitals + n_orbitals) * sizeof(basis_orbital));
				if (new_orbitals == NULL)
				{
					fprintf(stderr, "error: basis_load failed to allocate %d orbitals.\n", current_process_atom->norbitals);
					abort();
				}
				current_process_atom->orbitals = new_orbitals;

				for (int i = current_process_atom->norbitals, type = BASIS_ORBITAL_TYPE_S; i < current_process_atom->norbitals + n_orbitals; ++i, ++type)
				{
					// Set the orbital info
					current_process_atom->orbitals[i].ngtos = number_of_gto;
					current_process_atom->orbitals[i].shell = nshell;

					// n_orbitals is just the orbital type
					current_process_atom->orbitals[i].type = type;

					// Allocate the GTO terms
					current_process_atom->orbitals[i].gtos = malloc(number_of_gto * sizeof(basis_gto));
					if (current_process_atom->orbitals[i].gtos == NULL)
					{
						fprintf(stderr, "error: basis_load failed to allocate %d GTOs for %d GTOs.\n", number_of_gto, i);
						abort();
					}
					current_process_atom->orbitals[i].ngtos = number_of_gto;
				}

				// Read in the basis data
				for (int i = 0; i < number_of_gto; ++i)
				{
					// Read in the exponent (shared between GTOs)
					double exponent;
					if (fscanf(basis_file, "%lf", &exponent) != 1)
					{
						fprintf(stderr, "error: basis_load failed to parse basis: cannot read exponent for GTO %d.\n", i);
						abort();
					}

					// Read in the coefficients
					for (int j = current_process_atom->norbitals; j < current_process_atom->norbitals + n_orbitals; ++j)
					{
						double coefficient;
						if (fscanf(basis_file, "%lf", &coefficient) != 1)
						{
							fprintf(stderr, "error: basis_load failed to parse basis: cannot read coefficient for GTO %d.\n", j);
							abort();
						}

						// Set coefficients and exponents
						current_process_atom->orbitals[j].gtos[i].exponent = exponent;
						current_process_atom->orbitals[j].gtos[i].coefficient = coefficient;
					}
				}

				current_process_atom->norbitals += n_orbitals;
			}
			else
			{
				fprintf(stderr, "error: basis_load orbital specification line with wrong format: %s.\n", plinebuf);
				abort();
			}

			++nshell;
		}
		else if (term_number == 2)
		{
			// This is an atom description line
			// like: He     0
			char atom_name[BASIS_ATOM_NAME_SZBUF];
			take(sscanf(plinebuf, "%s", atom_name));

			int atom_number = atom_lookup(atom_name);
			if (atom_number == -1)
			{
				// Not found
				fprintf(stderr, "error: basis_load found invalid or not supported atom name: %s.\n", atom_name);
				abort();
			}

			if (basis_size < atom_number + 1)
			{
				// Expand the basis storage
				basis_atom* new_basis = realloc(basis, (atom_number + 1) * sizeof(basis_atom));
				if (new_basis == NULL)
				{
					// Unable to allocate more memory
					fprintf(stderr, "error: basis_load cannot allocate %d units for atom %s.\n", atom_number + 1, atom_name);
					abort();
				}

				basis = new_basis;
				basis_size = atom_number + 1;
			}

			// Remember the current atom
			nshell = 0;
			current_process_atom = basis + atom_number;
			current_process_atom->norbitals = 0;
			current_process_atom->orbitals = NULL;
		}
		else
		{
			fprintf(stderr, "error: basis_load found unexpected line in basis file: %s.\n", linebuf);
			abort();
		}
	}
}

void basis_print()
{
	for (int i = 0; i < basis_size; ++i)
	{
		printf("Basis for atom %d (%s):\n", i, atom_names[i]);

		for (int j = 0; j < basis[i].norbitals; ++j)
		{
			printf("  Orbital %d: Shell %d, Type %d\n", j, basis[i].orbitals[j].shell, basis[i].orbitals[j].type);

			for (int k = 0; k < basis[i].orbitals[j].ngtos; ++k)
			{
				printf("    GTO %d: Exp %.8lf Coeff %.8lf\n", k, basis[i].orbitals[j].gtos[k].exponent, basis[i].orbitals[j].gtos[k].coefficient);
			}
		}
	}
}

void basis_clean()
{
}
