// NK: Numerical Recipes for Lattice Boltzmann

#ifndef _NRLB_H
#define _NRLB_H

#include "Globals.h"

void lubksb(Real **a, int n, int *indx, Real b[]);
void ludcmp(Real **a, int n, int *indx, Real *d);

#endif
