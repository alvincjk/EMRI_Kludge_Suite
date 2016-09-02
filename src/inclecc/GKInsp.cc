#include <math.h>
#include <stdio.h>
#include "Globals.h"
#include "GKInsp.h"
#include "IEKG.h"
#include "CKG.h"
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

//
// Note that deltat_est is in physical time units!
//
GKInsp::GKInsp(TrajData EccTraj[], const int ne,
	       TrajData CircTraj[], const int nc,
	       const Real Mass, const Real Massratio, const Real spin,
	       Real dt_est, Real pinit, Real pfin) :
  Necc(ne), Ncirc(nc), deltat_est(fabs(dt_est)), 
  sign((int)(dt_est/fabs(dt_est))), M(Mass), Movermu(Massratio), a(spin)
{
  errencount=false;
  //
  ios::sync_with_stdio();
  //
  x = Realvector(1, 4);
  //
  int i;
  //
  // Take the trajectory data and store it in different structures;
  // this adds a bit of initial overhead but makes later
  // interpolations much easier.
  //
  if (Necc > 0) {
    t_ecc = Realvector(1, Necc);
    p_ecc = Realvector(1, Necc);
    e_ecc = Realvector(1, Necc);
    cosi_ecc = Realvector(1, Necc);
    E_ecc = Realvector(1, Necc);
    Lz_ecc = Realvector(1, Necc);
    Q_ecc = Realvector(1, Necc);
    p3_ecc = Realvector(1, Necc);
    p4_ecc = Realvector(1, Necc);
    beta_wilkins_ecc = Realvector(1, Necc);
    zedminus_ecc = Realvector(1, Necc);
    betazedplus_ecc = Realvector(1, Necc);
    for (i = 1; i <= Necc; i++) {
      t_ecc[i] = sign*EccTraj[i].t;
      p_ecc[i] = EccTraj[i].p;
      e_ecc[i] = EccTraj[i].ecc;
      cosi_ecc[i] = EccTraj[i].cosiota;
      E_ecc[i] = EccTraj[i].E;
      Lz_ecc[i] = EccTraj[i].Lz;
      Q_ecc[i] = EccTraj[i].Q;
      p3_ecc[i] = EccTraj[i].p3;
      p4_ecc[i] = EccTraj[i].p4;
      beta_wilkins_ecc[i] = EccTraj[i].beta_wilkins;
      zedminus_ecc[i] = EccTraj[i].zedminus;
      betazedplus_ecc[i] = EccTraj[i].betazedplus;
    }
  }
  if (Ncirc > 0) {
    t_circ = Realvector(1, Ncirc);
    p_circ = Realvector(1, Ncirc);
    cosi_circ = Realvector(1, Ncirc);
    E_circ = Realvector(1, Ncirc);
    Lz_circ = Realvector(1, Ncirc);
    Q_circ = Realvector(1, Ncirc);
    p3_circ = Realvector(1, Ncirc);
    p4_circ = Realvector(1, Ncirc);
    beta_wilkins_circ = Realvector(1, Ncirc);
    zedminus_circ = Realvector(1, Ncirc);
    betazedplus_circ = Realvector(1, Ncirc);
    for (i = 1; i <= Ncirc; i++) {
      t_circ[i] = sign*CircTraj[i].t;
      p_circ[i] = CircTraj[i].p;
      cosi_circ[i] = CircTraj[i].cosiota;
      E_circ[i] = CircTraj[i].E;
      Lz_circ[i] = CircTraj[i].Lz;
      Q_circ[i] = CircTraj[i].Q;
      p3_circ[i] = CircTraj[i].p3;
      p4_circ[i] = CircTraj[i].p4;
      beta_wilkins_circ[i] = CircTraj[i].beta_wilkins;
      zedminus_circ[i] = CircTraj[i].zedminus;
      betazedplus_circ[i] = CircTraj[i].betazedplus;
    }
  }
  //
  // OK, now that we've read in the data file, we want to rescale the
  // time series.  t_ecc and t_circ are read in radiation reaction
  // time units; we hit with M^2/\mu to rescale it up.
  //
  for (i = 1; i <= Necc; i++)
    t_ecc[i] *= M*Movermu;
  for (i = 1; i <= Ncirc; i++)
    t_circ[i] *= M*Movermu;
  //
  // Figure out the times corresponding to pinit and pfin.  Find the
  // indices that bracket those p values in the p arrays first.
  //
  Real frac;
  if (Necc != 0) { // Begins eccentric
    if (sign*pinit >= sign*p_ecc[1]) {
      pinit = p_ecc[1];
      tecc_init = t_ecc[1];
    } else { // Bracket
      i = 1;
      while (sign*pinit < sign*p_ecc[i]) i++;
      frac = (p_ecc[i] - pinit)/(p_ecc[i] - p_ecc[i - 1]);
      tecc_init = frac*t_ecc[i - 1] + (1. - frac)*t_ecc[i];
    }
  } else { // Begins circular
    if (sign*pinit >= sign*p_circ[1]) {
      pinit = p_circ[1];
      tcirc_init = t_circ[1];
    } else { // Bracket
      i = 1;
      while (sign*pinit < sign*p_circ[i]) i++;
      frac = (p_circ[i] - pinit)/(p_circ[i] - p_circ[i - 1]);
      tcirc_init = frac*t_circ[i - 1] + (1. - frac)*t_circ[i];
    }
  }
 
  if (Ncirc == 0) { // Ends eccentric
    if (sign*pfin <= sign*p_ecc[Necc]) {
      pfin = p_ecc[Necc];
      tecc_fin = t_ecc[Necc]; 
    } else { // Bracket
      i = Necc;
      while (sign*pfin > sign*p_ecc[i]) i--; // This overshoots a touch, so correct.
      i++;
      frac = (p_ecc[i] - pfin)/(p_ecc[i] - p_ecc[i - 1]);
      tecc_fin = frac*t_ecc[i - 1] + (1. - frac)*t_ecc[i];
    }
  } else { // Ends circular
    if (sign*pfin <= sign*p_circ[Ncirc]) {
      pfin = p_circ[Ncirc];
      tcirc_fin = t_circ[Ncirc];
    } else { // Bracket
      i = Ncirc;
      while (sign*pfin > sign*p_circ[i]) i--; // This overshoots a touch, so correct.
      i++;
      frac = (p_circ[i] - pfin)/(p_circ[i] - p_circ[i - 1]);
      tcirc_fin = frac*t_circ[i - 1] + (1. - frac)*t_circ[i];
    }
  }
  //
  // Be careful if we have both circular and eccentric inspirals.
  //
  if (Necc != 0 && Ncirc != 0) {
    tecc_fin = t_ecc[Necc];
    tcirc_init = t_circ[1];
  }
}

GKInsp::~GKInsp()
{
  if (Necc > 0) {
    free_Realvector(t_ecc, 1, Necc);
    free_Realvector(p_ecc, 1, Necc);
    free_Realvector(e_ecc, 1, Necc);
    free_Realvector(cosi_ecc, 1, Necc);
    free_Realvector(E_ecc, 1, Necc);
    free_Realvector(Lz_ecc, 1, Necc);
    free_Realvector(Q_ecc, 1, Necc);
    free_Realvector(p3_ecc, 1, Necc);
    free_Realvector(p4_ecc, 1, Necc);
    free_Realvector(beta_wilkins_ecc, 1, Necc);
    free_Realvector(zedminus_ecc, 1, Necc);
    free_Realvector(betazedplus_ecc, 1, Necc);
  }
  if (Ncirc > 0) {
    free_Realvector(t_circ, 1, Ncirc);
    free_Realvector(p_circ, 1, Ncirc);
    free_Realvector(cosi_circ, 1, Ncirc);
    free_Realvector(E_circ, 1, Ncirc);
    free_Realvector(Lz_circ, 1, Ncirc);
    free_Realvector(Q_circ, 1, Ncirc);
    free_Realvector(p3_circ, 1, Ncirc);
    free_Realvector(p4_circ, 1, Ncirc);
    free_Realvector(beta_wilkins_circ, 1, Ncirc);
    free_Realvector(zedminus_circ, 1, Ncirc);
    free_Realvector(betazedplus_circ, 1, Ncirc);
  }
}

Real GKInsp::LinInterp(const Real t, Real timearr[], Real dataarr[],
		       const int N)
{
  //
  // Sometimes need to advance ihi ...
  //
  while (timearr[ihi] < t && ihi < N)
    ihi++;
  //
  // ... and sometimes need to dial it back.
  //
  while (timearr[ihi - 1] > t && ihi > 1)
    ihi--;
  //
  const Real frac = (timearr[ihi] - t)/
    (timearr[ihi] - timearr[ihi - 1]);
  return(frac*dataarr[ihi - 1] + (1. - frac)*dataarr[ihi]);
}

void GKInsp::getorbitalelements(const Real t, Real *elts)
{
  if (DoingEcc) {
    elts[0]=LinInterp(sign*t, t_ecc, p_ecc, Necc);
    elts[1]=LinInterp(sign*t, t_ecc, e_ecc, Necc);
    elts[2]=LinInterp(sign*t, t_ecc, cosi_ecc, Necc);
  } else {
    elts[0]=LinInterp(sign*t, t_circ, p_circ, Ncirc);
    elts[1]=0.;
    elts[2]=LinInterp(sign*t, t_circ, cosi_circ, Ncirc);
  }
  return;
}

void GKInsp::getorbitalspeeds(const Real t, Real *speeds, Real psi, Real chi)
{
  Real p, ecc, E, Lz, p3, p4, beta_wilkins, zedminus, betazedplus;
  if (DoingEcc) {
    p = LinInterp(sign*t, t_ecc, p_ecc, Necc);
    ecc = LinInterp(sign*t, t_ecc, e_ecc, Necc);
    E = LinInterp(sign*t, t_ecc, E_ecc, Necc);
    Lz = LinInterp(sign*t, t_ecc, Lz_ecc, Necc);
    p3 = LinInterp(sign*t, t_ecc, p3_ecc, Necc);
    p4 = LinInterp(sign*t, t_ecc, p4_ecc, Necc);
    beta_wilkins = LinInterp(sign*t, t_ecc, beta_wilkins_ecc, Necc);
    zedminus = LinInterp(sign*t, t_ecc, zedminus_ecc, Necc);
    betazedplus = LinInterp(sign*t, t_ecc, betazedplus_ecc, Necc);
  } else {
    p = LinInterp(sign*t, t_circ, p_circ, Ncirc);
    ecc = 0.0;
    E = LinInterp(sign*t, t_circ, E_circ, Ncirc);
    Lz = LinInterp(sign*t, t_circ, Lz_circ, Ncirc);
    p3 = LinInterp(sign*t, t_circ, p3_circ, Ncirc);
    p4 = LinInterp(sign*t, t_circ, p4_circ, Ncirc);
    beta_wilkins = LinInterp(sign*t, t_circ, beta_wilkins_circ, Ncirc);
    zedminus = LinInterp(sign*t, t_circ, zedminus_circ, Ncirc);
    betazedplus = LinInterp(sign*t, t_circ, betazedplus_circ, Ncirc);
  }
  //
  const Real cpsi = cos(psi);
  const Real spsi = sin(psi);
  const Real r = p/(1. + ecc*cpsi);
  //
  const Real cchi = cos(chi);
  const Real schi = sin(chi);
  const Real z = zedminus*cchi*cchi;
  //
  const Real delta = Kerr::Delta(r, a);
  const Real r2pa2 = r*r + a*a;
  const Real gamma = E*(r2pa2*r2pa2/delta - a*a) - 2.*r*a*Lz/delta;
  const Real gampa2Ez = gamma + a*a*E*z;
  //
  const Real psitmp1 = (p - p3) - ecc*(p + p3*cpsi);
  const Real psitmp2 = (p - p4) + ecc*(p - p4*cpsi);
  //
  const Real dpsidt = sqrt((1. - E)*(1. + E)*psitmp1*psitmp2)/
    (gampa2Ez*(1. - ecc)*(1. + ecc)*M);
  //
  const Real dchidt = sqrt(betazedplus - beta_wilkins*z)/(gampa2Ez*M);
  //
  speeds[0] = r*r*ecc*spsi*dpsidt/p;
  speeds[1] = sqrt(zedminus/(1.-z))*schi*dchidt;
  speeds[2] = (Lz/(1. - z) + 2.*a*E*r/delta - a*a*Lz/delta)/(gampa2Ez*M);

  return;
}

void GKInsp::getorbitalconsts(const Real t, Real *elts)
{
  if (DoingEcc) {
    elts[0]=LinInterp(sign*t, t_ecc, E_ecc, Necc);
    elts[1]=LinInterp(sign*t, t_ecc, Lz_ecc, Necc);
    elts[2]=LinInterp(sign*t, t_ecc, Q_ecc, Necc);
  } else {
    elts[0]=LinInterp(sign*t, t_circ, E_circ, Ncirc);
    elts[1]=LinInterp(sign*t, t_circ, Lz_circ, Ncirc);
    elts[2]=LinInterp(sign*t, t_circ, Q_circ, Ncirc);
  }
  return;
}

Real GKInsp::rfunc(const Real t, const Real psi)
{
  if (DoingEcc) {
    return(LinInterp(sign*t, t_ecc, p_ecc, Necc)/
	   (1. + LinInterp(sign*t, t_ecc, e_ecc, Necc)*cos(psi)));
  } else {
    return(LinInterp(sign*t, t_circ, p_circ, Ncirc));
  }
}

Real GKInsp::rstarfunc(const Real t, const Real psi)
{
  double rval,rp,rm,afact;
  if (DoingEcc) {
    rval=(LinInterp(sign*t, t_ecc, p_ecc, Necc)/
	   (1. + LinInterp(sign*t, t_ecc, e_ecc, Necc)*cos(psi)));
  } else {
    rval=(LinInterp(sign*t, t_circ, p_circ, Ncirc));
  }
  rp=Kerr::rplus(a);
  rm=Kerr::rminus(a);
  afact=1./sqrt(1.-a*a);
  return(rval+(rp*log(fabs(rval/rp-1.))-rm*log(fabs(rval/rm-1.)))*afact);
}

Real GKInsp::psifunc(const Real t, const Real psi, const Real phi)
{
  double rval,rp,rm,afact;
  if (DoingEcc) {
    rval=(LinInterp(sign*t, t_ecc, p_ecc, Necc)/
	   (1. + LinInterp(sign*t, t_ecc, e_ecc, Necc)*cos(psi)));
  } else {
    rval=(LinInterp(sign*t, t_circ, p_circ, Ncirc));
  }
  rp=Kerr::rplus(a);
  rm=Kerr::rminus(a);
  afact=1./sqrt(1.-a*a);
  return(phi+a*log(fabs((rval-rp)/(rval-rm)))*0.5*afact);
}

Real GKInsp::costhetafunc(const Real t, const Real chi)
{
  Real chi_here = chi;
  
  if (chi_here > 0.) {
    while (chi_here > 2.*M_PI)
      chi_here -= 2.*M_PI;
  }
  else {
    while (chi_here < 0.)
      chi_here += 2.*M_PI;
  }

  Real sgn;
  /* Use the convention that when cos chi is negative, cos theta is also negative. */
  if (chi_here < 0.5*M_PI || chi_here > 1.5*M_PI)
    sgn = 1.;
  else
    sgn = -1.;

  const Real cchi = cos(chi_here);
  Real E, Lz, Q, zedminus;
  if (DoingEcc) {
    E = LinInterp(sign*t, t_ecc, E_ecc, Necc);
    Lz = LinInterp(sign*t, t_ecc, Lz_ecc, Necc);
    Q = LinInterp(sign*t, t_ecc, Q_ecc, Necc);
    zedminus = LinInterp(sign*t, t_ecc, zedminus_ecc, Necc);
  } else {
    E = LinInterp(sign*t, t_circ, E_circ, Ncirc);
    Lz = LinInterp(sign*t, t_circ, Lz_circ, Ncirc);
    Q = LinInterp(sign*t, t_circ, Q_circ, Ncirc);
    zedminus = LinInterp(sign*t, t_circ, zedminus_circ, Ncirc);
  }
  //cout << "GKInsp z- = " << zedminus << endl; 
  if (zedminus*cchi*cchi > 2.) zedminus=0.;
  return(sgn*sqrt(fabs(zedminus*cchi*cchi)));
}

void GKInsp::TakeAStep(Real x[], const Real t1, const Real t2)
{
  const Real EPS = EPSILON*1000.;
  const Real t1call=sign*t1;
  const Real t2call=sign*t2;
  odeint(x, 4, t1call, t2call, EPS, t2call - t1call, EPSILON);
}

void GKInsp::odeint(Real ystart[], const int nvar, Real x1,
		    Real x2, const Real eps, Real h1,
		    const Real hmin)
{
  int nstp,i;
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
    GKInsp_derivs(x, y, dydx);
    for (i = 1; i <= nvar; i++) 
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
    if ((x + h - x2)*(x + h - x1) > 0.0) h = x2 - x;
    bsstep(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext);
    if (errencount)
      return; 
    if ((x - x2)*(x2 - x1) >= 0.0) {
      for (i = 1; i <= nvar; i++) ystart[i] = y[i];
      free_Realvector(dydx, 1, nvar);
      free_Realvector(y, 1, nvar);
      free_Realvector(yscal, 1, nvar);
      return;
    }
    if (fabs(hnext) <= hmin) {
      //Die("Step size too small in odeint");
      errencount=true;
      return;
    }
    h = hnext;
  }
  //Die("Too many steps in routine odeint");
  errencount=true;
  return;
}

void GKInsp::bsstep( Real y[], Real dydx[], int nv, Real *xx,
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
    for (k = 1; k <= kmax; k++){
      xnew = (*xx) + h;
      if (xnew == (*xx)) {
	//Die("step size underflow in bsstep in GKInsp.cc");
	errencount=true;
	return;
      }
      mmid(ysav, dydx, nv, *xx, h, nseq[k], yseq);
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
}

//
// The geodesic equations are encoded here.  Note the
// notation: x[1] = psi; x[2] = chi; x[3] = phi; x[4] = tau
//
void GKInsp::GKInsp_derivs(Real t, Real *x, Real *dxdt)
{
  Real p, ecc, E, Lz, Q, p3, p4, beta_wilkins, zedminus, betazedplus;
  if (DoingEcc) {
    p = LinInterp(t, t_ecc, p_ecc, Necc);
    ecc = LinInterp(t, t_ecc, e_ecc, Necc);
    E = LinInterp(t, t_ecc, E_ecc, Necc);
    Lz = LinInterp(t, t_ecc, Lz_ecc, Necc);
    Q = LinInterp(t, t_ecc, Q_ecc, Necc);
    p3 = LinInterp(t, t_ecc, p3_ecc, Necc);
    p4 = LinInterp(t, t_ecc, p4_ecc, Necc);
    beta_wilkins = LinInterp(t, t_ecc, beta_wilkins_ecc, Necc);
    zedminus = LinInterp(t, t_ecc, zedminus_ecc, Necc);
    betazedplus = LinInterp(t, t_ecc, betazedplus_ecc, Necc);
  } else {
    p = LinInterp(t, t_circ, p_circ, Ncirc);
    ecc = 0.0;
    E = LinInterp(t, t_circ, E_circ, Ncirc);
    Lz = LinInterp(t, t_circ, Lz_circ, Ncirc);
    Q = LinInterp(t, t_circ, Q_circ, Ncirc);
    p3 = LinInterp(t, t_circ, p3_circ, Ncirc);
    p4 = LinInterp(t, t_circ, p4_circ, Ncirc);
    beta_wilkins = LinInterp(t, t_circ, beta_wilkins_circ, Ncirc);
    zedminus = LinInterp(t, t_circ, zedminus_circ, Ncirc);
    betazedplus = LinInterp(t, t_circ, betazedplus_circ, Ncirc);
  }
  //  
  const Real psi = x[1];
  const Real cpsi = cos(psi);
  const Real r = p/(1. + ecc*cpsi);
  //
  const Real chi = x[2];
  const Real cchi = cos(chi);
  const Real z = zedminus*cchi*cchi;
  //
  const Real sigma = Kerr::Sigma(r, a, z);
  const Real delta = Kerr::Delta(r, a);
  const Real r2pa2 = r*r + a*a;
  const Real gamma = E*(r2pa2*r2pa2/delta - a*a) - 2.*r*a*Lz/delta;
  const Real gampa2Ez = gamma + a*a*E*z;
  //
  const Real psitmp1 = (p - p3) - ecc*(p + p3*cpsi);
  const Real psitmp2 = (p - p4) + ecc*(p - p4*cpsi);
  //
  // NOTE: these are all derivatives wrt time and therefore we need
  // to worry about mass normalization!
  //
  dxdt[1] = sign*sqrt((1. - E)*(1. + E)*psitmp1*psitmp2)/
    (gampa2Ez*(1. - ecc)*(1. + ecc)*M);
  //
  dxdt[2] = sign*sqrt(betazedplus - beta_wilkins*z)/(gampa2Ez*M);
  //
  dxdt[3] = sign*(Lz/(1. - z) + 2.*a*E*r/delta - a*a*Lz/delta)/(gampa2Ez*M);
  //
  dxdt[4] = sign*sigma/(gampa2Ez*M);
  //
}

void GKInsp::GKInsp_QuadrupoleMoments(Real t, Real *x, Real *Quad)
{
  Real p, ecc, E, Lz, Q, p3, p4, beta_wilkins, zedminus, betazedplus;
  if (DoingEcc) {
    p = LinInterp(sign*t, t_ecc, p_ecc, Necc);
    ecc = LinInterp(sign*t, t_ecc, e_ecc, Necc);
    E = LinInterp(sign*t, t_ecc, E_ecc, Necc);
    Lz = LinInterp(sign*t, t_ecc, Lz_ecc, Necc);
    Q = LinInterp(sign*t, t_ecc, Q_ecc, Necc);
    p3 = LinInterp(sign*t, t_ecc, p3_ecc, Necc);
    p4 = LinInterp(sign*t, t_ecc, p4_ecc, Necc);
    beta_wilkins = LinInterp(sign*t, t_ecc, beta_wilkins_ecc, Necc);
    zedminus = LinInterp(sign*t, t_ecc, zedminus_ecc, Necc);
    betazedplus = LinInterp(sign*t, t_ecc, betazedplus_ecc, Necc);
  } else {
    p = LinInterp(sign*t, t_circ, p_circ, Ncirc);
    ecc = 0.0;
    E = LinInterp(sign*t, t_circ, E_circ, Ncirc);
    Lz = LinInterp(sign*t, t_circ, Lz_circ, Ncirc);
    Q = LinInterp(sign*t, t_circ, Q_circ, Ncirc);
    p3 = LinInterp(sign*t, t_circ, p3_circ, Ncirc);
    p4 = LinInterp(sign*t, t_circ, p4_circ, Ncirc);
    beta_wilkins = LinInterp(sign*t, t_circ, beta_wilkins_circ, Ncirc);
    zedminus = LinInterp(sign*t, t_circ, zedminus_circ, Ncirc);
    betazedplus = LinInterp(sign*t, t_circ, betazedplus_circ, Ncirc);
  }
  //
  const Real psi = x[1];
  const Real cpsi = cos(psi);
  const Real spsi = sin(psi);
  const Real r = p/(1. + ecc*cpsi);
  //
  const Real chi = x[2];
  const Real cchi = cos(chi);
  const Real schi = sin(chi);
  const Real z = zedminus*cchi*cchi;
  //
  const Real cphi=cos(x[3]);
  const Real sphi=sin(x[3]);
  const Real cth = costhetafunc(t,chi);
  const Real sth=sqrt(1.-cth*cth);
  //
  const Real sigma = Kerr::Sigma(r, a, z);
  const Real delta = Kerr::Delta(r, a);
  const Real r2pa2 = r*r + a*a;
  const Real gamma = E*(r2pa2*r2pa2/delta - a*a) - 2.*r*a*Lz/delta;
  const Real gampa2Ez = gamma + a*a*E*z;
  //
  const Real psitmp1 = (p - p3) - ecc*(p + p3*cpsi);
  const Real psitmp2 = (p - p4) + ecc*(p - p4*cpsi);
  //
  // NOTE: these are all derivatives wrt time and therefore we need
  // to worry about mass normalization!
  //
  const Real dpsidt = sqrt((1. - E)*(1. + E)*psitmp1*psitmp2)/
    (gampa2Ez*(1. - ecc)*(1. + ecc)*M);
  //
  const Real dchidt = sqrt(betazedplus - beta_wilkins*z)/(gampa2Ez*M);
  //
  const Real dphidt = (Lz/(1. - z) + 2.*a*E*r/delta - a*a*Lz/delta)/(gampa2Ez*M);
  //
  // Work out second derivatives, since these are what we need to compute the quadrupole moments for the waveforms.
  // Note that we include the correct mass normalisation, so these are also second derivatives with respect to the 
  // actual time.
  // 
  const Real drdt = r*r*ecc*spsi*dpsidt/p;
  // I think this is correct sign for theta derivative, but can check.
  const Real dthdt = sqrt(zedminus/(1.-z))*schi*dchidt;
  //
  const Real dgampa2Ezdt = 2.*(2.*E*r*r2pa2-a*Lz-(r-1.)*(E*r2pa2*r2pa2-2.*a*r*Lz)/delta)*drdt/delta
    -2.*a*a*E*cchi*schi*zedminus*dchidt;
  //
  const Real d2psidt2=(ecc*spsi*dpsidt*(p3/psitmp1+p4/psitmp2)/2.-dgampa2Ezdt/gampa2Ez)*dpsidt;
  //
  const Real d2chidt2=(beta_wilkins*zedminus*cchi*schi*dchidt/(betazedplus - beta_wilkins*z)-dgampa2Ezdt/gampa2Ez)*dchidt;
  //
  const Real d2phidt2=(-2.*dchidt*cchi*schi*zedminus*Lz/((1.-z)*(1.-z))+2.*(a*E/delta+(r-1.)*(a*a*Lz-2.*a*E*r)/(delta*delta))*
	     drdt)/(gampa2Ez*M)-dgampa2Ezdt*dphidt/gampa2Ez;
  //
  const Real d2rdt2 = r*r*ecc*(spsi*d2psidt2+r*dpsidt*dpsidt*(cpsi+ecc*(1.+spsi*spsi))/p)/p;
  const Real d2thdt2 = sqrt(zedminus)*(cchi*(1.-zedminus*schi*schi/(1.-z))*dchidt*dchidt+schi*d2chidt2)/sqrt(1.-z);
  // 
  const Real s2th=2.*sth*cth;
  const Real c2th=cth*cth-sth*sth;
  const Real s2ph=2.*sphi*cphi;
  const Real c2ph=cphi*cphi-sphi*sphi;
  const Real d2r2=(r*d2rdt2+drdt*drdt);
  const Real d2rphi=drdt*dphidt;
  const Real d2rth=drdt*dthdt;
  const Real d2phith=dphidt*dthdt;
  const Real d2phi2=dphidt*dphidt;
  const Real d2th2=dthdt*dthdt;
  //
  Quad[1]=2.*d2r2*(sth*sth*cphi*cphi-1./3.)+4.*r*(s2th*d2rth*cphi*cphi-sth*sth*s2ph*d2rphi)
    +r*r*(-2.*s2th*s2ph*d2phith-sth*sth*(2.*c2ph*d2phi2+s2ph*d2phidt2)+cphi*cphi*(2.*c2th*d2th2+d2thdt2*s2th));
  //
  Quad[2]=2.*d2r2*(sth*sth*sphi*sphi-1./3.)+4.*r*(s2th*d2rth*sphi*sphi+sth*sth*s2ph*d2rphi)
    +r*r*(2.*s2th*s2ph*d2phith+sth*sth*(2.*c2ph*d2phi2+s2ph*d2phidt2)+sphi*sphi*(2.*c2th*d2th2+d2thdt2*s2th));
  //
  Quad[3]=d2r2*s2ph*sth*sth+2.*r*(s2th*s2ph*d2rth+2.*c2ph*sth*sth*d2rphi)
    +r*r*(2.*s2th*c2ph*d2phith+sth*sth*(c2ph*d2phidt2-2.*s2ph*d2phi2)+s2ph*(c2th*d2th2+s2th*d2thdt2/2.));
  //
  Quad[4]=d2r2*s2th*cphi+2.*r*(2.*c2th*cphi*d2rth-s2th*sphi*d2rphi)
    +r*r*(cphi*(c2th*d2thdt2-2.*s2th*d2th2)-2.*c2th*sphi*d2phith-s2th*(cphi*d2phi2+sphi*d2phidt2)/2.);
  //
  Quad[5]=d2r2*s2th*sphi+2.*r*(2.*c2th*sphi*d2rth+s2th*cphi*d2rphi)
    +r*r*(sphi*(c2th*d2thdt2-2.*s2th*d2th2)+2.*c2th*cphi*d2phith-s2th*(sphi*d2phi2-cphi*d2phidt2)/2.);
  //
}

void GKInsp::mmid(Real y[], Real dydx[], int nvar, Real xs,
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
  GKInsp_derivs(x, yn, yout);

  h2 = 2.0*h;
  for (n = 2; n <= nstep; n++) {
    for (i = 1; i <= nvar; i++) {
      swap = ym[i] + h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    GKInsp_derivs(x, yn, yout);
  }
  for (i = 1; i<= nvar; i++)
    yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);
  free_Realvector(yn, 1, nvar);
  free_Realvector(ym, 1, nvar);
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
