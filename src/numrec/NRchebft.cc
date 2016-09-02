#include "Globals.h"
#include "NRUtil.h"

void chebft(Real a, Real b, Real c[], int n, Real (*func)(Real))
{
  int k, j;
  Real fac, bpa, bma, *f;
  
  f = Realvector(0, n - 1);
  bma = 0.5*(b - a);
  bpa = 0.5*(b + a);
  for (k = 0; k < n; k++) {
    Real y = cos(M_PI*(k + 0.5)/n);
    f[k] = (*func)(y*bma + bpa);
  }
  fac = 2.0/n;
  for (j = 0; j < n; j++) {
    double sum = 0.0;
    for (k = 0; k < n; k++)
      sum += f[k]*cos(M_PI*j*(k + 0.5)/n);
    c[j] = fac*sum;
  }
  free_Realvector(f, 0, n - 1);
}
