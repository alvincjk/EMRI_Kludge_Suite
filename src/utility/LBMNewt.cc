#include <math.h>
#include "Globals.h"
#include "LB.h"
#include "NRLB.h"
#include "NRUtil.h"

void LB::mnewt(int ntrial, Real x[], int n, Real tolx, Real tolf)
{

  int k, i, *indx;
  Real errx, errf, d, *fvec, **fjac, *p;
  
  indx = ivector(1, n);
  p = Realvector(1, n);
  fvec = Realvector(1, n);
  fjac = Realmatrix(1, n, 1, n);
  for (k = 1; k <= ntrial; k++) {
    func_and_jac(x, n, fvec, fjac);
    errf = 0.0;
    for (i = 1; i <= n; i++) errf += fabs(fvec[i]);
    if (errf <= tolf) {
      free_Realmatrix(fjac, 1, n, 1, n);
      free_Realvector(fvec, 1, n);
      free_Realvector(p, 1, n);
      free_ivector(indx, 1, n);
      return;
    }
    for (i = 1; i <= n; i++) p[i] = -fvec[i];
    ludcmp(fjac, n, indx, &d);
    lubksb(fjac, n, indx, p);
    errx = 0.0;
    for (i = 1; i <= n; i++) {
      errx += fabs(p[i]);
      x[i] += p[i];
    }
    if (errx <= tolx) {
      free_Realmatrix(fjac, 1, n, 1, n);
      free_Realvector(fvec, 1, n);
      free_Realvector(p, 1, n);
      free_ivector(indx, 1, n);
      return;
    }
  }
  free_Realmatrix(fjac, 1, n, 1, n);
  free_Realvector(fvec, 1, n);
  free_Realvector(p, 1, n);
  free_ivector(indx, 1, n);
  return;
}
