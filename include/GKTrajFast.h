// AAK: Fast generic Kerr inspirals

#ifndef _GKTRAJFAST_H
#define _GKTRAJFAST_H

#include <math.h>
#include "Globals.h"
#include "IEKG.h"
#include "GKR.h"
#include "CKG.h"

#define ECC_ONE 0.9995

//
// All times that are used in this routine are in "radiation reaction
// time".  To convert to "physical time" (the time in which waveforms
// are measured) multiply by M^2/\mu.
//
class GKTrajFast {
public:
  GKTrajFast(const Real cosiotastart, const Real spin);

  int Eccentric(const Real dtstart, TrajData traj[],
		int StepsMax, int & StepsTaken);
  int Circular(const Real dtstart, TrajData traj[],
	       int StepsMax, int & StepsTaken);

  void TakeAStep(Real x[], const Real t1);

  Real p, ecc, cosiota, t;

  Real a, dt, dt_orig;

  Real **p_lso, *e_lso, *csi_lso;

  void Param_derivs(Real t, Real *x, Real *dxdt);
  Real P_lso(const Real csi, const Real e);

private:
  int LSO_Max_Index,LSO_csi_Max_Index;

  int odeint(Real ystart[], const int nvar, const Real x1, const Real x2,
	     const Real eps, const Real h1, const Real hmin);
  int bsstep(Real y[], Real dydx[], int nv, Real *xx, Real htry, Real eps,
	     Real yscal[], Real *hdid, Real *hnext);  
  int mmid(Real y[], Real dydx[], int nvar, Real xs, Real htot, int nstep,
	   Real yout[]);

  void Rroots(const Real Ehere, const Real Lzhere, const Real Qhere,
	      Real &arr1, Real &arr2, Real &arr3, Real &arr4);

  IEKG *iekg_tmp;
  CKG *ckg_tmp;
};
#endif
