// NK: Numerical Recipes for CKG

#ifndef _NRCKG_H
#define _NRCKG_H

#include "Globals.h"

Real elle(const Real phi, const Real ak);
Real ellf(const Real phi, const Real ak);
Real ellpi(const Real phi, const Real en, const Real ak);

Real rd(const Real x, const Real y, const Real z);
Real rf(const Real x, const Real y, const Real z);
Real rj(const Real x, const Real y, const Real z, const Real p);

#endif
