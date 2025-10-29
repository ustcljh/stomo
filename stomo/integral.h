#pragma once

// Do the S_uv overlap integral. u and v are ids in orbitals
double integral_s(int u, int v);

// Do the T_uv kinetic integral. u and v are ids in orbitals
double integral_t(int u, int v);

// Do the V_uv potential integral. u and v are ids in orbitals
double integral_v(int u, int v);

// Do the (c1c2|c3c4) four center integral.
double integral_four(int c1, int c2, int c3, int c4);
