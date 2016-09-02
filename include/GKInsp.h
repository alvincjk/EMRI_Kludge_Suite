// NK: Inclined, eccentric Kerr inspirals

#ifndef _GKINSP_H
#define _GKINSP_H

#include <math.h>
#include "GKG.h"
#include "Globals.h"

void rzextr(int iest, Real xest, Real yest[], Real yz[], Real dy[],
	    Real **d, Real *x, int nv);

class GKInsp {
public:
  bool errencount;
  GKInsp(TrajData EccTraj[], const int ne,
	 TrajData CircTraj[], const int nc,
	 const Real Mass, const Real Massratio, const Real spin,
	 Real dt_est, Real pinit, Real pfin);
  ~GKInsp();
  //
  Real LinInterp(const Real t, Real timearr[],
		 Real dataarr[], const int N);
  //
  void TakeAStep(Real x[], const Real t1, const Real t2);
  void GKInsp_QuadrupoleMoments(Real t, Real *x, Real *Quad);
  //
  Real *x; // Holder for integrated data
  //
  void getorbitalelements(const Real t, Real *elts);
  void getorbitalspeeds(const Real t, Real *speeds, Real psi, Real chi);
  void getorbitalconsts(const Real t, Real *elts);
  Real rfunc(const Real t, const Real psi);
  Real rstarfunc(const Real t, const Real psi);
  Real psifunc(const Real t, const Real psi, const Real phi);
  Real costhetafunc(const Real t, const Real chi);
  //
  // Trajectory data, re-stored to make interpolations work better.
  // Note that the time is stored in *physical* time --- we multiply
  // by M^2/\mu as soon as it is passed in.
  //
  int Necc, Ncirc;
  Real *t_ecc;
  Real *p_ecc, *e_ecc, *cosi_ecc;
  Real *E_ecc, *Lz_ecc, *Q_ecc;
  Real *p3_ecc, *p4_ecc;
  Real *beta_wilkins_ecc, *zedminus_ecc, *betazedplus_ecc;
  Real *t_circ;
  Real *p_circ, *cosi_circ;
  Real *E_circ, *Lz_circ, *Q_circ;
  Real *p3_circ, *p4_circ;
  Real *beta_wilkins_circ, *zedminus_circ, *betazedplus_circ;
  //
  // The following are used for searching data arrays.  The time
  // variables are all in physical time!
  //
  //
  Real deltat_est;
  Real tecc_init, tecc_fin;
  Real tcirc_init, tcirc_fin;
  int ihi;
  //
  int DoingEcc; // Distinguish eccentric versus circular inspiral.
  int sign;
  void GKInsp_derivs(Real t, Real *x, Real *dxdt);

private:
  Real M, Movermu, a;

  void odeint(Real ystart[], const int nvar, const Real x1, const Real x2,
	      const Real eps, const Real h1, const Real hmin);
  void bsstep(Real y[], Real dydx[], int nv, Real *xx, Real htry, Real eps,
	      Real yscal[], Real *hdid, Real *hnext);
  void mmid(Real y[], Real dydx[], int nvar, Real xs, Real htot, int nstep,
	    Real yout[]);
};

#endif
