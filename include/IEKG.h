// NK: Inclined, eccentric Kerr geodesics

#ifndef _IEKG_H
#define _IEKG_H

#include <math.h>
#include "Globals.h"
#include "GKG.h"
#include "NRRoot.h"

void rzextr(int iest, Real xest, Real yest[], Real yz[], Real dy[],
	    Real **d, Real *x, int nv);

class IEKG {
public:
  IEKG(const Real semilat, const Real eee, const Real costilt,
       const Real spin);
  ~IEKG();

  void Constants();
  void Rroots(const Real Ehere, const Real Lzhere, const Real Qhere,
	      Real &arr1, Real &arr2, Real &arr3, Real &arr4,
	      int &stable);
  void mnewt(int ntrial, Real x[], int n, Real tolx, Real tolf);
  void TakeAStep(Real x[], const Real t1, const Real t2);
  Real rfunc(const Real psi);
  Real costhetafunc(const Real chi);
  void Frequencies(Real *Omega);

  Real p, ecc, cosiota;
  Real peri, ap, tansqiota, taniota;
  Real p3, p4;

  GKG *gkg;

  Real a;
  Real E, Lz, Q;
  Real alpha_wilkins, beta_wilkins;
  Real zedplus, zedminus;
  Real betazedplus;
  Real thetamin, thetamax;
  Real r1, r2, r3, r4;
  //
  // Stable ==  1: stable orbit.
  // Stable ==  0: marginally stable orbit.
  // Stable == -1: unstable orbit.
  //
  int Stable;
  void IEKG_derivs(Real t, Real *x, Real *dxdt);               

private:
  void odeint(Real ystart[], const int nvar, const Real x1, const Real x2,
	      const Real eps, const Real h1, const Real hmin);
  void bsstep(Real y[], Real dydx[], int nv, Real *xx, Real htry, Real eps,
	      Real yscal[], Real *hdid, Real *hnext);
  void mmid(Real y[], Real dydx[], int nvar, Real xs, Real htot, int nstep,
	    Real yout[]);
  void func_and_jac(Real *x, int n, Real *fvec, Real **fjac);
};
#endif
