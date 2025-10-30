#pragma once

#include "matrix.h"

// Do the plotting of wavefunction, exporting the data to GNUPlot-ready text. Coefficient matrix is C and it takes the [index] of that
void plot_wfunc_gnuplot(const char* filename, matrix C, int index);
