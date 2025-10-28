#pragma once

// Names of the atoms, [1] = "H" etc.
extern const char* atom_names[];

// Get atom number from name
int atom_lookup(const char* name);
