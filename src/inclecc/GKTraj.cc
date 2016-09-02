#include <math.h>
#include <stdio.h>
#include "Globals.h"
#include "GKTraj.h"
#include "IEKG.h"
#include "CKG.h"
#include "LB.h"
#include "NRUtil.h"

#define MAXSTP 10000000
#define TINY 1.0e-30
#define EPSILON 1.e-15
#define KMAXX 8
#define IMAXX (KMAXX + 1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define SCALMX 0.1

GKTraj::GKTraj(const Real cosiotastart, const Real spin)
  : cosiota(cosiotastart), a(spin)
{
  //
  // The main thing the constructor does is store the coordinates of
  // the LSO at this inclination angle.
  //
  LSO_Max_Index = 200;
  LSO_csi_Max_Index=100;
  p_lso = Realmatrix(1, LSO_csi_Max_Index,1,LSO_Max_Index);
  e_lso = Realvector(1, LSO_Max_Index);
  csi_lso = Realvector(1, LSO_csi_Max_Index);
  //
  const Real r_ipro = Kerr::isco_pro(a);
  const Real r_iret = Kerr::isco_ret(a);
  const Real r_horz = Kerr::rplus(a);
  //
  Real p_l, p_h, p_g;
  int i,j;
  //
  //cout << "In GKTraj, cos i = " << cosiota << ", 1.-cosi = " << 1.-cosiota << endl;
  for (j=1; j <= LSO_csi_Max_Index; j++) {
    csi_lso[j] = 1.-((Real)(j - 1))*2./((Real)(LSO_csi_Max_Index - 1));    
  for (i = 1; i <= LSO_Max_Index; i++) {
    e_lso[i] = ((Real)(i - 1))*ECC_ONE/((Real)(LSO_Max_Index - 1));
    if (e_lso[i] == 0.) {
      //
      // If circular and equatorial, it's pretty trivial.
      //
      /* NB This can cause problems if the offset in the ORBIT_I direction is small and we are centred on an equatorial orbit.
	    The computer then takes a positive and negative offset to be the same as the zero offset, so we do not get 
	    off-diagonal elements that are 2*gxx, but actually just gxx. Must be wary of this in the future.... */

      if (fabs(csi_lso[j] - 1.0) < 1.e-38) {
	//cout << "Correcting cos i " << endl;
	p_l = p_h = r_ipro;
	csi_lso[j] = 1.0;
      } else if (fabs(csi_lso[j] + 1.0) < 1.e-38) {
	p_l = p_h = r_iret;
	csi_lso[j] = -1.0;
      } else {
	//
	// If circular and inclined, use the circular LB code.
	//
	// Inner bound: prograde orbit
	//
	p_l = r_ipro;
	Real *x_l; x_l = Realvector(1, 3);
	x_l[1] = Kerr::Eeqpro(p_l, a);
	x_l[2] = Kerr::Lzeqpro(p_l, a);
	x_l[3] = 0.;
	//
	// Outer bound: retrograde orbit
	//
	p_h = r_iret;
	Real *x_h; x_h = Realvector(1, 3);
	x_h[1] = Kerr::Eeqret(p_h, a);
	x_h[2] = Kerr::Lzeqret(p_h, a);
	x_h[3] = 0.;
	//
	Real *x_g; x_g = Realvector(1, 3);
	while (p_h - p_l > 1.e-12) {
	  p_g = 0.5*(p_h + p_l);
	  x_g[1] = 0.5*(x_l[1] + x_h[1]);
	  x_g[2] = 0.5*(x_l[2] + x_h[2]);
	  x_g[3] = 0.5*(x_l[3] + x_h[3]);
	  LB lb(p_g, a);
	  const Real tolx = 1.e-14, tolf = 1.e-14;
	  lb.mnewt(100, x_g, 3, tolx, tolf);
	  Real cosiota_g = x_g[2]/sqrt(x_g[2]*x_g[2] + x_g[3]);
	  if (csi_lso[j] < cosiota_g) { // Move farther out (larger iota)
	    p_l = p_g;
	    x_l[1] = x_g[1];
	    x_l[2] = x_g[2];
	    x_l[3] = x_g[3];
	  } else {
	    p_h = p_g;
	    x_h[1] = x_g[1];
	    x_h[2] = x_g[2];
	    x_h[3] = x_g[3];
	  }
	}
	p_g = 0.5*(p_h + p_l);
	p_l = p_h = p_g;
      }
    } else {
      // Set p_l to slightly below previous answer, p_h to above.
      p_l = p_g - 0.01;
      p_h = p_g + 0.1;
    }
    // Now, bisect
    while (p_h - p_l > 1.e-14) {
      p_g = 0.5*(p_h + p_l);
      IEKG iekg(p_g, e_lso[i], csi_lso[j], a);
      if (iekg.Stable == 1) {
	// Guess is too far out: move p_h in
	p_h = p_g;
      } else if (iekg.Stable == -1) {
	// Guess is too far in: move p_l out
	p_l = p_g;
      } else if (iekg.Stable == 0) {
	// Land right on it ... it could happen!
	p_h = p_l = p_g;
      }
    }
    p_g = 0.5*(p_h + p_l);
    p_lso[j][i] = p_g;
   
  }
  }
}

//
// This routine computes the eccentric part of an inspiral.
//
// Possible return values:
//     2: traj reaches maximum allowed steps before plunging or
//        circularizing.
//     1: traj circularizes before plunging.
//     0: traj plunges before circularizing.
//    -1: traj has unstable initial conditions.
//    -2: traj has gone beyond validity of the hybrid approximation.
//
// Note that dtstart is in radiation reaction time units!
//
int GKTraj::Eccentric(const Real dtstart, TrajData traj[],
		      int StepsMax, int & StepsTaken)
{
  //cout << "Called GKTraj Eccentric" << endl;
  dt = dtstart;
  //
  // Let's keep initial value of ecc as a reference point.
  //
  Real ecc_ref = ecc;
  //
  // x[1] = p; x[2] = e; x[3] = cos(iota).
  //
  Real *x; x = Realvector(1, 3);
  x[1] = p; x[2] = ecc; x[3] = cosiota;
  t = 0;
  //
  // This will be used in some testing code in the main loop.
  //
  Real *dxdt_est; dxdt_est = Realvector(1, 3);
  //
  int returnflag, DONE = 0;
  //
  // Test whether the initial conditions are OK.
  //
  if (p < P_lso(cosiota,ecc)) {
    //cout << "In GKTraj, e = " << ecc << ", cos(iota) = " << cosiota << ", p_LSO = " << P_lso(cosiota,ecc) << endl;
    DONE = 1;
    returnflag = -1;
  }
  //
  // If the trajectory wants ecc = 0, the two equations [R(peri) = 0,
  // R(ap) = 0] don't provide separate conditions.  If so, return 1.
  //
  if (ecc == 0) {
    DONE = 1;
    returnflag = 1;
  }
  //
  // Begin main inspiral loop.
  //
  int i = 0;
  while (!DONE) {
    //cout << "Starting loop" << endl;
    iekg_tmp = new IEKG(x[1], x[2], x[3], a);
    //
    // Find roots of R(r)
    //
    Real r1, r2, r3, r4;
    Rroots(iekg_tmp->E, iekg_tmp->Lz, iekg_tmp->Q, r1, r2, r3, r4);
    const Real p3 = r3*(1 - x[2]);
    const Real p4 = r4*(1 + x[2]);
    //
    // Quantities used to determine the range in theta
    //
    const Real alpha_wilkins = iekg_tmp->Q + iekg_tmp->Lz*iekg_tmp->Lz;
    const Real beta_wilkins = a*a*(1. - iekg_tmp->E*iekg_tmp->E);
    const Real tmp = alpha_wilkins + beta_wilkins;
    Real zedplus, zedminus, betazedplus;
    if (beta_wilkins > EPSILON) {
      zedplus = (tmp + sqrt(tmp*tmp - 4.*iekg_tmp->Q*beta_wilkins));
      zedplus /= 2.*beta_wilkins;
      betazedplus = (tmp + sqrt(tmp*tmp - 4.*iekg_tmp->Q*beta_wilkins))/2.;
      zedminus = (tmp - sqrt(tmp*tmp - 4.*iekg_tmp->Q*beta_wilkins));
      zedminus /= 2.*beta_wilkins;
      if (fabs(zedminus) < 1.e-15) zedminus = 0.0;
    } else {
      zedplus = tmp/beta_wilkins - (iekg_tmp->Q/tmp)*
	(1. - 2.*beta_wilkins/(tmp*tmp));
      betazedplus = tmp - (iekg_tmp->Q*beta_wilkins/tmp)*
	(1. - 2.*beta_wilkins/(tmp*tmp));
      zedminus = (iekg_tmp->Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
    }
    //cout << "In GKTraj, zedminus = " << zedminus << endl;
    //
    // Increment step count, store data, check whether we've reach
    // maxsteps.
    //
    i++;
    traj[i].t = t;
    traj[i].p = x[1];
    traj[i].ecc = x[2];
    traj[i].cosiota = x[3];
    traj[i].E = iekg_tmp->E;
    traj[i].Lz = iekg_tmp->Lz;
    traj[i].Q = iekg_tmp->Q;
    traj[i].p3 = p3;
    traj[i].p4 = p4;
    traj[i].beta_wilkins = beta_wilkins;
    traj[i].zedminus = zedminus;
    traj[i].betazedplus = betazedplus;
    if (i == StepsMax) {
      DONE = 1;
      returnflag = 2;
    }
    //
    // Save p and ecc values prior to stepping --- will use this to
    // see whether they are evolving in a sensical way.
    //
    Real p_prev = x[1], ecc_prev = x[2];
    //
    // Step.
    //
    //cout << "Trying to take a step" << endl;
    TakeAStep(x, t);
    //cout << x[1] << " " << x[2] << " " << x[3] << endl;
    //cout << "GKTraj has taken a step" << endl;
    //
    // If p is growing, the approx is no good.  Exit cleanly. But only if going forward in time!!
    //
    if (dt*(x[1] - p_prev) > 0.) {
      DONE = 1;
      returnflag = -2;
      cout << "pgrowth" << endl;
    }
    //
    // If ecc is heading towards 1, the approx is no good. Exit cleanly. Likewise!!
    //
    if (x[2] > 0.98 && dt*(x[2] - ecc_prev) > 0.) {
      DONE = 1;
      returnflag = -2;
      cout << "e growth" << endl;
    }
    //
    // Check whether the next step is likely to take us past the LSO.
    // If yes, shrink the timestep.  Note that if we're essentially
    // already on the unstable point this can lead to an infinite
    // loop, so make sure that dt doesn't get too small while we're at
    // it.
    //
    Param_derivs(t, x, dxdt_est);
    while (x[1] + dt*dxdt_est[1] < P_lso(x[3],x[2]) && fabs(dt) > 1.e-9) {
      dt /= 2.;
    }
    //
    // Also check that we don't step into negative eccentricity.
    //
    while (x[2] + dt*dxdt_est[2] < 0 && fabs(dt) > 1.e-9) {
      dt /= 2.;
    }
    //
    // If the eccentricity has shrunk to 1/2 of the reference value,
    // we are possibly circularizing.  Decrease the timestep and
    // update the reference val.
    //
    if (x[2] < 0.5*ecc_ref) {
      ecc_ref = x[2];
      dt /= 2.;
    }
    //
    // If the timestep or the eccentricity has shrunk too much, we're
    // out of here.
    //
    if (fabs(dt) < 1.e-8 || x[2] < 1.e-6) {
      DONE = 1;
      if (x[2] < 1.e-3 && iekg_tmp->Stable == 1) {
	returnflag = 1;
      } else {
	returnflag = 0;
      }
    }
    delete iekg_tmp;
    t += dt;
  }
  StepsTaken = i;
  //
  // Copy over the state of things at the end of it all.
  //
  p = x[1];
  ecc = x[2];
  cosiota = x[3];
  free_Realvector(x, 1, 3);
  free_Realvector(dxdt_est, 1, 3);
  return(returnflag);
}

//
// This routine computes the circular part of an inspiral.
//
// Possible return values:
//     2: traj reaches maximum allowed steps before plunging.
//     1: traj runs until plunging.
//    -1: traj has unstable initial conditions.
//
// Note that dtstart is in radiation reaction time units!
//
int GKTraj::Circular(const Real dtstart, TrajData traj[],
		     int StepsMax, int & StepsTaken)
{
  dt = dtstart;
  //
  // x[1] = p; x[2] = e; x[3] = cos(iota).
  //
  Real *x; x = Realvector(1, 3);
  //
  // We will use this to estimate whether we're about to step over the
  // LSO or not.
  //
  Real *dxdt_est; dxdt_est = Realvector(1, 3);
  x[1] = p; x[2] = 0.0; x[3] = cosiota;
  //
  int returnflag, DONE = 0;
  //
  Real r_LSO = P_lso(cosiota,0.);
  //
  // Sanity check: we'd best be outside the LSO right now!
  //
  if (x[1] <= r_LSO) {
    DONE = 1;
    returnflag = -1;
  }
  //
  // Begin main inspiral loop.
  //
  int i = 0;
  while (!DONE) {
    //
    ckg_tmp = new CKG(2, 0, x[3], USING_cosiota, x[1], a);
    //
    // Find roots of R(r)
    //
    Real r1, r2, r3, r4;
    Rroots(ckg_tmp->E, ckg_tmp->Lz, ckg_tmp->Q, r1, r2, r3, r4);
    const Real p3 = r3;
    const Real p4 = r4;
    //
    // Quantities used to determine the range in theta
    //
    const Real alpha_wilkins = ckg_tmp->Q + ckg_tmp->Lz*ckg_tmp->Lz;
    const Real beta_wilkins = a*a*(1. - ckg_tmp->E*ckg_tmp->E);
    const Real tmp = alpha_wilkins + beta_wilkins;
    Real zedplus, zedminus, betazedplus;
    if (beta_wilkins > EPSILON) {
      zedplus = (tmp + sqrt(tmp*tmp - 4.*ckg_tmp->Q*beta_wilkins));
      zedplus /= 2.*beta_wilkins;
      betazedplus = (tmp + sqrt(tmp*tmp - 4.*ckg_tmp->Q*beta_wilkins))/2.;
      zedminus = (tmp - sqrt(tmp*tmp - 4.*ckg_tmp->Q*beta_wilkins));
      zedminus /= 2.*beta_wilkins;
      if (fabs(zedminus) < 1.e-15) zedminus = 0.0;
    } else {
      zedplus = tmp/beta_wilkins - (ckg_tmp->Q/tmp)*
	(1. - 2.*beta_wilkins/(tmp*tmp));
      betazedplus = tmp - (ckg_tmp->Q*beta_wilkins/tmp)*
	(1. - 2.*beta_wilkins/(tmp*tmp));
      zedminus = (ckg_tmp->Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
    }
    //
    // Increment step count, store data, check whether we've reach
    // maxsteps.
    //
    i++;
    traj[i].t = t;
    traj[i].p = x[1];
    traj[i].ecc = 0.;
    traj[i].cosiota = x[3];
    traj[i].E = ckg_tmp->E;
    traj[i].Lz = ckg_tmp->Lz;
    traj[i].Q = ckg_tmp->Q;
    traj[i].p3 = p3;
    traj[i].p4 = p4;
    traj[i].beta_wilkins = beta_wilkins;
    traj[i].zedminus = zedminus;
    traj[i].betazedplus = betazedplus;
    if (i == StepsMax) {
      DONE = 1;
      returnflag = 2;
    }
    TakeAStep(x, t);
    t += dt;
    //
    // Code to test ending conditions.  If the next step is likely
    // to take us beyond the LSO, we need to decrease the timestep.
    //
    Param_derivs(t, x, dxdt_est);
    while(x[1] + dt*dxdt_est[1] < r_LSO && dt > 1.e-9)
      dt /= 2.;
    //
    // When the timestep shrinks too much, we're out of here.
    //
    if (dt < 1.e-8) {
      DONE = 1;
      returnflag = 1;
    }
    delete ckg_tmp;
  }
  StepsTaken = i;
  free_Realvector(x, 1, 3);
  free_Realvector(dxdt_est, 1, 3);
  return(returnflag);
}

//
// Computes the p coordinate of the lso at a particular value of e;
// uses linear interpolation.
//
Real GKTraj::P_lso(const Real csi, const Real e)
{
  //
  // Bracket.
  //
  int csihi = 2;
  while (csi_lso[csihi] > csi) csihi++;  
  int ehi=2;
  while (e_lso[ehi] < e) ehi++;
  //
  // Interpolate
  //
  Real efrac = (e_lso[ehi] - e)/(e_lso[ehi] - e_lso[ehi - 1]);
  Real csifrac = (csi_lso[csihi] - csi)/(csi_lso[csihi] - csi_lso[csihi - 1]);

  return(csifrac*(efrac*p_lso[csihi-1][ehi-1]+(1.-efrac)*p_lso[csihi-1][ehi]) + 
(1.-csifrac)*(efrac*p_lso[csihi][ehi-1]+(1.-efrac)*p_lso[csihi][ehi]));
  //return(p_lso[ihi - 1]*frac + p_lso[ihi]*(1. - frac));
}

void GKTraj::TakeAStep(Real x[], const Real t1)
{
  const Real EPS = EPSILON*1000.;
  int DONE = 0;
  int integrate_value;
  while (!DONE) {
    integrate_value = odeint(x, 3, t1, t1 + dt, EPS, dt, EPS);
    if (integrate_value == 1) {
      //
      // We're out of here...
      //
      DONE = 1;
    } else if (integrate_value == 0) {
      //
      // Need to reduce timestep and try again.
      //
      dt /= 10.;
    } else if (integrate_value == -1) {
      //
      // We've done the best possible.  Set dt = 0 and get out of here.
      //
      dt = 0.0;
      DONE = 1;
    }
  }
}

//
// Returns 1 if all is OK
// Returns 0 if we need to reduce the timestep and try again.
// Returns -1 if the timestep goes to zero and we should bail.
//
int GKTraj::odeint(Real ystart[], const int nvar, const Real x1,
		   const Real x2, const Real eps, const Real h1,
		   const Real hmin)
{
  int nstp, i;
  Real xsav, x, hnext, hdid, h;
  Real *yscal;
  Real *y, *dydx;
  
  yscal = Realvector(1, nvar);
  y = Realvector(1, nvar);
  dydx = Realvector(1, nvar);
  x = x1;
  h = SIGN(h1, x2 - x1);
  for (i = 1; i <= nvar; i++) y[i] = ystart[i];
  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    Param_derivs(x, y, dydx);
    for (i = 1; i <= nvar; i++)
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
    if ((x + h - x2)*(x + h - x1) > 0.0) h = x2 - x;
    int bsstep_val = bsstep(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext);
    if (bsstep_val == -1) {
      //
      // bsstep can't resolve any better.  Most likely, we've integrated
      // in to zero eccentricity.
      //
      free_Realvector(dydx, 1, nvar);
      free_Realvector(y, 1, nvar);
      free_Realvector(yscal, 1, nvar);
      return(-1);
    }
    if ((x - x2)*(x2 - x1) >= 0.0) {
      for (i = 1; i <= nvar; i++) ystart[i] = y[i];
      free_Realvector(dydx, 1, nvar);
      free_Realvector(y, 1, nvar);
      free_Realvector(yscal, 1, nvar);
      return(1);
    }
    if (fabs(hnext) <= hmin || bsstep_val == -2) {
      free_Realvector(dydx, 1, nvar);
      free_Realvector(y, 1, nvar);
      free_Realvector(yscal, 1, nvar);
      return(0);
    }
    h = hnext;
  }
  return(0);
}

//
// Returns 1 if all ends OK.
// Returns -1 if stepsize shrinks too much.
// Returns -2 if mmid gives it trouble (indicating odeint needs
//  to reduce the stepsize).
//
int GKTraj::bsstep(Real y[], Real dydx[], int nv, Real *xx,
		   Real htry, Real eps, Real yscal[], Real *hdid,
		   Real *hnext)
{
  int i, iq, k, kk, km;
  static int first = 1, kmax, kopt;
  static Real epsold = -1.0, xnew;
  Real **d;
  Real *x;
  Real eps1, errmax, fact, h, red, scale, work, wrkmin, xest;
  Real *err;
  Real *yerr, *ysav, *yseq;
  static Real a[IMAXX + 1];
  static Real alf[KMAXX + 1][KMAXX + 1];
  static int nseq[IMAXX + 1] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
  int reduct, exitflag = 0;

  d = Realmatrix(1, KMAXX, 1, KMAXX);
  err = Realvector(1, KMAXX);
  x = Realvector(1, KMAXX);
  yerr = Realvector(1, nv);
  ysav = Realvector(1, nv);
  yseq = Realvector(1, nv);
  if (eps != epsold) {
    *hnext = xnew = -1.0e29;
    eps1 = SAFE1*eps;
    a[1] = nseq[1] + 1;
    for (k = 1; k <= KMAXX; k++) a[k + 1] = a[k] + nseq[k + 1];
    for (iq = 2; iq <= KMAXX; iq++) {
      for (k = 1; k < iq; k++)
	alf[k][iq] = pow(eps1, (a[k + 1] - a[iq + 1])/
			 ((a[iq + 1] - a[1] + 1.0)*(2*k + 1)));
    }
    epsold = eps;
    for (kopt = 2; kopt < KMAXX; kopt++)
      if (a[kopt + 1] > a[kopt]*alf[kopt - 1][kopt]) break;
    kmax = kopt;
  }
  h = htry;
  for (i = 1; i <= nv; i++) ysav[i] = y[i];
  if (*xx != xnew || h != (*hnext)) {
    first = 1;
    kopt = kmax;
  }
  reduct = 0;
  for (;;) {
    for (k = 1; k <= kmax; k++) {
      xnew = (*xx) + h;
      if (xnew == (*xx)) {
	free_Realvector(yseq, 1, nv);
	free_Realvector(ysav, 1, nv);
	free_Realvector(yerr, 1, nv);
	free_Realvector(x, 1, KMAXX);
	free_Realvector(err, 1, KMAXX);
	free_Realmatrix(d, 1, KMAXX, 1, KMAXX);
	return(-1);
      }
      int mmid_val = mmid(ysav, dydx, nv, *xx, h, nseq[k], yseq);
      if (mmid_val == -1) return(-2);
      xest = SQR(h/nseq[k]);
      rzextr(k, xest, yseq, y, yerr, d, x, nv);
      if (k != 1) {
	errmax = TINY;
	for (i = 1; i <= nv; i++)
	  errmax = FMAX(errmax, fabs(yerr[i]/yscal[i]));
	errmax /= eps;
	km = k - 1;
	err[km] = pow(errmax/SAFE1, 1.0/(2*km + 1));
      }
      if (k != 1 && (k >= kopt - 1 || first)) {
	if (errmax < 1.0) {
	  exitflag = 1;
	  break;
	}
	if (k == kmax || k == kopt + 1) {
	  red = SAFE2/err[km];
	  break;
	}
	else if (k == kopt && alf[kopt - 1][kopt] < err[km]) {
	  red = 1.0/err[km];
	  break;
	}
	else if (kopt == kmax && alf[km][kmax - 1] < err[km]) {
	  red = alf[km][kmax - 1]*SAFE2/err[km];
	  break;
	}
	else if (alf[km][kopt] < err[km]) {
	  red = alf[km][kopt - 1]/err[km];
	  break;
	}
      }
    }
    if (exitflag) break;
    red = FMIN(red, REDMIN);
    red = FMAX(red, REDMAX);
    h *= red;
    reduct = 1;
  }
  *xx = xnew;
  *hdid = h;
  first = 0;
  wrkmin = 1.0e35;
  for (kk = 1; kk <= km; kk++) {
    fact = FMAX(err[kk], SCALMX);
    work = fact*a[kk + 1];
    if (work < wrkmin) {
      scale = fact;
      wrkmin = work;
      kopt = kk + 1;
    }
  }
  *hnext = h/scale;
  if (kopt >= k && kopt != kmax && !reduct) {
    fact = FMAX(scale/alf[kopt - 1][kopt], SCALMX);
    if (a[kopt + 1]*fact <= wrkmin) {
      *hnext = h/fact;
      kopt++;
    }
  }
  free_Realvector(yseq, 1, nv);
  free_Realvector(ysav, 1, nv);
  free_Realvector(yerr, 1, nv);
  free_Realvector(x, 1, KMAXX);
  free_Realvector(err, 1, KMAXX);
  free_Realmatrix(d, 1, KMAXX, 1, KMAXX);
  return(1);
}

//
// The parameter derivatives.
//
void GKTraj::Param_derivs(Real t, Real *x, Real *dxdt)
{
  const Real ppp = x[1];
  const Real eee = x[2];
  const Real cosincl = x[3];
  //
  // These are only need to map (p, e, i) -> (E, Lz, Q)
  //
  IEKG *inclined;
  CKG *circular;
  //
  Real Ehere, Lzhere, Qhere;
  if (eee != 0.0) {
    inclined = new IEKG(ppp, eee, cosincl, a);
    Ehere = inclined->E;
    Lzhere = inclined->Lz;
    Qhere = inclined->Q;
    delete inclined;
  } else {
    circular = new CKG(2, 0, cosincl, USING_cosiota, ppp, a);
    Ehere = circular->E;
    Lzhere = circular->Lz;
    Qhere = circular->Q;
    delete circular;
  }
  //
  GKR gkr(ppp, eee, Ehere, Lzhere, Qhere, a);
  //
  dxdt[1] = gkr.pdot;
  //
  dxdt[2] = gkr.edot;
  //
  dxdt[3] = gkr.cosidot;
}

//
// Returns 1 if all is OK.
// Returns -1 if the code tries to put the semilatus rectum negative.
//
int GKTraj::mmid(Real y[], Real dydx[], int nvar, Real xs,
		 Real htot, int nstep, Real yout[])
{
  int n, i;
  Real swap, *ym, *yn;
  Real x, h2, h;
  
  ym = Realvector(1, nvar);
  yn = Realvector(1, nvar);
  h = htot/nstep;
  for (i = 1; i <= nvar; i++) {
    ym[i] = y[i];
    yn[i] = y[i] + h*dydx[i];
  }
  x = xs + h;
  Param_derivs(x, yn, yout);

  h2 = 2.0*h;
  for (n = 2; n <= nstep; n++) {
    for (i = 1; i <= nvar; i++) {
      swap = ym[i] + h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    //
    // If at any point in this process, yn[1] gets significantly
    // inside the LSO --- say, 5% beyond --- we're stepping too far
    // and finding crap.  We need to decrease the stepsize.  Test this
    // by comparing with the LSO for the circular orbit.
    //
    if (yn[1] - P_lso(yn[3],0.) < -0.05*P_lso(yn[3],0.)) return(-1);
    //
    // Likewise, if the eccentricity has decided to become negative,
    // we're stepping into the land of land of serpents and dragons.
    // Quoth Chef: "Let's get the fudge out, children".
    //
    if (yn[2] < 0.) return(-1);
    x += h;
    Param_derivs(x, yn, yout);
  }
  for (i = 1; i<= nvar; i++)
    yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);
  free_Realvector(yn, 1, nvar);
  free_Realvector(ym, 1, nvar);
  return(1);
}

void GKTraj::Rroots(const Real Ehere, const Real Lzhere, const Real Qhere,
		    Real &arr1, Real &arr2, Real &arr3, Real &arr4)
{
  int N_R = 4;
  Complex *R_cofs, *R_roots;
  Real ome2 = (1. - Ehere)*(1. + Ehere);
  //
  // Store coefficients of polynomial for R in array, find roots.
  //
  R_cofs = Complexvector(0, N_R);
  R_cofs[0] = -a*a*Qhere/ome2;
  R_cofs[1] = 2.*(Qhere + (a*Ehere - Lzhere)*(a*Ehere - Lzhere))/ome2;
  R_cofs[2] = -(a*a*ome2 + Lzhere*Lzhere + Qhere)/ome2;
  R_cofs[3] = 2./ome2;
  R_cofs[4] = -1.;
  R_roots = Complexvector(1, N_R);
  zroots(R_cofs, N_R, R_roots, 1, EPSILON);
  free_Complexvector(R_cofs, 0, N_R);
  //
  // Remap the roots to what I've got in my notes.  Note that with
  // the current scheme for inputting parameters we NEVER get an
  // imaginary part to the roots (well, maybe a small one due to
  // roundoff error).
  //
  arr1 = R_roots[4].real();
  arr2 = R_roots[3].real();
  arr3 = R_roots[2].real();
  arr4 = R_roots[1].real();
  //
  // Stability should not be an issue in this code, so cut out all
  // of that crap.
  //
  free_Complexvector(R_roots, 1, N_R);
}

#undef MAXSTP
#undef TINY
#undef EPSILON
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef SCALMX
