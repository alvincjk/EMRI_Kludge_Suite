// NK: Numerical Recipes for IEKG

#ifndef _NRIEKG_H
#define _NRIEKG_H

#include "Globals.h"

void lubksb(Real **a, int n, int *indx, Real b[]);
void ludcmp(Real **a, int n, int *indx, Real *d);

#endif
