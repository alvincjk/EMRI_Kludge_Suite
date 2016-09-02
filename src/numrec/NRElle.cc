#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

Real elle(const Real phi, const Real ak)
{
  Real rd(const Real x, const Real y, const Real z);
  Real rf(const Real x, const Real y, const Real z);
  
  const Real s=sin(phi);
  const Real cc=SQR(cos(phi));
  const Real q=(1.0-s*ak)*(1.0+s*ak);

  return s*(rf(cc, q, 1.0)-(SQR(s*ak))*rd(cc, q, 1.0)/3.0);
}
