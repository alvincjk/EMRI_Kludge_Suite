#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

Real pythag(const Real a, const Real b)
{
  const Real absa = fabs(a);
  const Real absb = fabs(b);
  if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}
