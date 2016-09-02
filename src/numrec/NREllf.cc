#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

Real ellf(const Real phi, const Real ak)
{
  Real rf(const Real x, const Real y, const Real z);
  const Real s = sin(phi);

  if (phi > M_PI/2.)
    return (-1.*s*rf(SQR(cos(phi)), (1.0-s*ak)*(1.0+s*ak), 1.0) + 2.*rf(0.,(1.0-ak)*(1.0+ak), 1.0));
  else
    return s*rf(SQR(cos(phi)), (1.0-s*ak)*(1.0+s*ak), 1.0);
}
