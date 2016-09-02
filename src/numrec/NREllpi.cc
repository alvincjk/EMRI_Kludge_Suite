#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

Real ellpi(const Real phi, const Real en, const Real ak)
{
  Real rf(const Real x, const Real y, const Real z);
  Real rj(const Real x, const Real y, const Real z,
	  const Real p);

  const Real cc = SQR(cos(phi));
  const Real s = sin(phi);
  const Real enss = en*s*s;
  const Real q = (1.0-s*ak)*(1.0+s*ak);

  return s*(rf(cc, q, 1.0)-enss*rj(cc, q, 1.0, 1.0+enss)/3.0);
}
