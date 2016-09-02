#include <math.h>
#include "Globals.h"
#include "LB.h"

void LB::func_and_jac(Real *x, int n, Real *fvec, Real **fjac)
{
  const Real E = x[1];
  const Real Lz = x[2];
  const Real Q = x[3];

  fvec[1] = gkg->RFunc(r, a, 0., E, Lz, Q);
  fvec[2] = gkg->dr_RFunc(r, a, 0., E, Lz, Q);
  fvec[3] = gkg->ddr_RFunc(r, a, 0., E, Lz, Q);
  //
  // No separate calls for the Jacobian elements ...
  // they're fairly simple, though.
  //
  const Real tmp1 = r*r + a*a;
  const Real tmp2 = Lz - a*E;

  fjac[1][1] = 2.*r*(E*r*tmp1 - 2.*a*tmp2);
  fjac[1][2] = 2.*r*(2.*tmp2 - Lz*r);
  fjac[1][3] = 2.*r  - tmp1;
  fjac[2][1] = 4.*(E*r*(a*a + 2.*r*r) - a*tmp2);
  fjac[2][2] = 4.*(tmp2 - Lz*r);
  fjac[2][3] = 2.*(1. - r);
  fjac[3][1] = 4.*E*(a*a + 6.*r*r);
  fjac[3][2] = -4.*Lz;
  fjac[3][3] = -2.;
}
