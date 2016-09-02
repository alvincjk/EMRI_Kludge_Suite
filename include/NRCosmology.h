// NK: Numerical Recipes for Cosmology

#ifndef _NRCosmology_H
#define _NRCosmology_H

#include "Globals.h"

Real elle(const Real phi, const Real ak);
Real ellf(const Real phi, const Real ak);
Real ellpi(const Real phi, const Real en, const Real ak);

Real rd(const Real x, const Real y, const Real z);
Real rf(const Real x, const Real y, const Real z);
Real rj(const Real x, const Real y, const Real z, const Real p);

void sncndn(Real uu, Real emmc, Real *sn, Real *cn, Real *dn);

#endif
