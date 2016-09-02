// NK: Numerical Recipes for Mismatch

#ifndef _NRMM_H
#define _NRMM_H

#include "Globals.h"

void tqli(Real d[], Real e[], const int n, Real **z);
void svdcmp(Real **a, int m, int n, Real *w, Real **v);
void tred2(Real **a, const int n, Real d[], Real e[]);
void gaussj(Real **a, int n, Real **b, int m);
void ludcmp(Real **a, int n, int *indx, Real *d);
Real ran2(long *idum);

#endif
