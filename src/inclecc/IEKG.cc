#include <math.h>
#include "Globals.h"
#include "IEKG.h"
#include "GKG.h"
#include "NRUtil.h"
#include "NRCKG.h"

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

IEKG::IEKG(const Real semilat, const Real eee, const Real costilt,
	   const Real spin) :
  p(semilat), ecc(eee), cosiota(costilt), a(spin)
{
  //
  // Test whether ecc = 0 has been passed in, which doesn't work well:
  // in that case the two equations [R(peri) = 0, R(ap) = 0] don't
  // provide separate conditions.
  //
  if (ecc == 0)
    Die("Zero eccentricity?  Use the circular code, punky.");
  //
  // Compute peribothron, apbothron
  //
  peri = p/(1. + ecc);
  ap = p/(1. - ecc);
  tansqiota=(1. - cosiota)*(1. + cosiota)/(cosiota*cosiota);
  taniota=sqrt(tansqiota);
  //
  // Need a GKG object for some things.
  //
  gkg = new GKG;
  //
  // Compute E, Lz, Q
  //
  Constants();
  //
  // Quantities used to determine the range in theta
  //
  alpha_wilkins = Q + Lz*Lz;
  beta_wilkins = a*a*(1. - E*E);
  const Real tmp = alpha_wilkins + beta_wilkins;
  if (beta_wilkins > EPSILON) {
    zedplus = (tmp + sqrt(tmp*tmp - 4.*Q*beta_wilkins));
    zedplus /= 2.*beta_wilkins;
    betazedplus = (tmp + sqrt(tmp*tmp - 4.*Q*beta_wilkins))/2.;
    zedminus = (tmp - sqrt(tmp*tmp - 4.*Q*beta_wilkins));
    zedminus /= 2.*beta_wilkins;
    if (fabs(zedminus) < 1.e-15) zedminus = 0.0;
  } else {
    zedplus = tmp/beta_wilkins - (Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
    betazedplus = tmp - (Q*beta_wilkins/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
    zedminus = (Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
  }
  thetamin = acos(sqrt(zedminus));
  thetamax = M_PI - thetamin;
  cosiota = Lz/sqrt(alpha_wilkins);
  //
  // Find roots of R(r)
  //
  Rroots(E, Lz, Q, r1, r2, r3, r4, Stable);
  p3 = r3*(1 - ecc);
  p4 = r4*(1 + ecc);
  //
  // Final check: if the orbit wants to go inside the horizon, mark
  // the orbit as unstable.  Note we don't want to 'Die()' on this since
  // the codes which map the orbital parameter space will get botched in
  // that case.
  //
  if (peri < Kerr::rplus(a)) Stable = -1;
}

IEKG::~IEKG()
{
  delete gkg;
}

void IEKG::Constants()
{
  //
  // First, assemble the initial guess for the constants, using the
  // algorithm I came up with on page 8 of my notes.
  //
  const Real semimaj = peri/(1. - ecc);
  Real Epro = Kerr::Eeqpro(semimaj, a);
  Real Eret = Kerr::Eeqret(semimaj, a);
  if (isnan(Epro) || Epro > 1.-1.e-6) Epro = 1.-1.e-6;
  if (isnan(Eret) || Eret > 1.) Eret = 1.;
  const Real Eguess = 0.5*(cosiota + 1.)*Epro + 0.5*(1. - cosiota)*Eret;
  const Real Lzguess = cosiota*sqrt((1. - ecc)*(1. + ecc)/(2. - 2.*Eguess));

  const Real tolx = 1.e-14;
  const Real tolf = 1.e-14;
  Real *x; x = Realvector(1, 2);

  x[1] = Eguess; x[2] = Lzguess;

  mnewt(200, x, 2, tolx, tolf);
  
  E = x[1]; Lz = x[2]; Q = Lz*Lz*tansqiota;
  free_Realvector(x, 1, 2);
}

void IEKG::Rroots(const Real Ehere, const Real Lzhere, const Real Qhere,
		  Real &arr1, Real &arr2, Real &arr3, Real &arr4, int & stable)
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
  // Check stability.  If arr1 agrees with apbothron and arr2 agrees
  // with peribothron, it is stable.
  //
  if (fabs(arr1 - ap) < 1e-7 && fabs(arr2 - peri) < 1e-7) {
    stable = 1;
    //
    // Further, if arr3 and arr2 agree, it is marginally stable.
    //
    if (fabs(arr2 - arr3) < 1e-7)
      stable = 0;
  }
  //
  // If the roots are inconsistent, it was not a stable orbit.
  //
  if (fabs(arr2 - peri) > 1e-7)
    stable = -1;

  free_Complexvector(R_roots, 1, N_R);
}

Real IEKG::rfunc(const Real psi)
{
  return(p/(1. + ecc*cos(psi)));
}

Real IEKG::costhetafunc(const Real chi)
{
  Real chi_here = chi;
  while (chi_here > 2.*M_PI)
    chi_here -= 2.*M_PI;

  Real sign;
  if (chi_here < 0.5*M_PI || chi_here > 1.5*M_PI)
    sign = -1.;
  else
    sign = 1.;

  const Real cchi = cos(chi_here);
  return(sign*sqrt(fabs(zedminus)*cchi*cchi));
}

/* Uses the expressions derived in Schmidt (2002) to compute the three frequencies that characterize the 
   geodesic. Note that these are in units of 1/M, since the IEKG class does not know the mass of the 
   central black hole. */

void IEKG::Frequencies(Real *Omega)
{
  Real thk,rk,rH,rC,Xint,Yint,Zint,Wint,ElK,ElE,ElEth,ElPiR4,ElPiRH,ElPiRC,gamma,lambda,fact,rtfact,art,
    DeltaR4,WfactorI,WfactorII;

  rk=sqrt((ap-peri)*(r3-r4)/((ap-r3)*(peri-r4)));
  thk=sqrt(zedminus/zedplus);
  fact=p/((1.-ecc*ecc)*sqrt(1.-E*E));
  rtfact=sqrt((ap-r3)*(peri-r4));
  art=sqrt(1.-a*a);
  rH=1.+art;
  rC=1.-art;

  ElK=ellf(M_PI/2.,thk);
  ElEth=elle(M_PI/2.,thk);
  ElE=elle(M_PI/2.,rk);
  ElPiR4=(ap-r4)*ellpi(M_PI/2.,-(peri-ap)/(peri-r4),rk)/rtfact;
  ElPiRH=(ap-r4)*ellpi(M_PI/2.,-(ap-peri)*(rH-r4)/((peri-r4)*(ap-rH)),rk)/((ap-rH)*(rH-r4)*rtfact*art);
  ElPiRC=(ap-r4)*ellpi(M_PI/2.,-(ap-peri)*(rC-r4)/((peri-r4)*(ap-rC)),rk)/((ap-rC)*(rC-r4)*rtfact*art);
  DeltaR4=1./(r4*r4-2.*r4+a*a);

  Xint=2.*fact*ellf(M_PI/2.,rk)/rtfact;
  Yint=fact*((ap+peri+r3+r4)*ElPiR4+rtfact*ElE)+(r4*r4-0.5*(ap-r4)*(peri-r4))*Xint;
  Zint=fact*((2.*E*rH-a*Lz)*a*ElPiRH + a*(a*Lz-2.*E*rC)*ElPiRC)+(Lz*r4+2.*a*E-2.*Lz)*r4*Xint*DeltaR4;

  WfactorI=8.*E-2.*a*Lz;
  WfactorII=4.*E*a*a;
  Wint=E*Yint+(2.*E*r4+4.*E+(WfactorI*r4-WfactorII)*DeltaR4)*Xint +fact*(4.*E*ElPiR4+(WfactorI*rH-WfactorII)*ElPiRH+(WfactorII-WfactorI*rC)*ElPiRC);
  
  lambda=(Yint+a*a*zedplus*Xint)*ElK-a*a*zedplus*Xint*elle(M_PI/2.,thk);
  gamma=((Wint+a*a*zedplus*E*Xint)*ElK-a*a*zedplus*E*Xint*ElEth)/lambda;

  /* Radial frequency. We define frequencies as the number of complete cycles per second, which is 1/T_r etc. 
     This is slightly different to the omega_r etc. defined by Schmidt, which are the number of radians per
     second of the angle variables. */
  Omega[0]=p*ElK/(2.*gamma*(1.-ecc*ecc)*lambda);
  /* Theta frequency */
  Omega[1]=sqrt(betazedplus)*Xint/(4.*gamma*lambda);
  /* Phi frequency */
  Omega[2]=((Zint-Lz*Xint)*ElK+Lz*Xint*ellpi(M_PI/2.,-zedminus,thk))/(2.*M_PI*gamma*lambda);
}

void IEKG::TakeAStep(Real x[], const Real t1, const Real t2)
{
  const Real EPS = EPSILON*1000.;
  odeint(x, 4, t1, t2, EPS, t2 - t1, EPSILON);
}

void IEKG::odeint(Real ystart[], const int nvar, const Real x1,
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
    IEKG_derivs(x, y, dydx);
    for (i = 1; i <= nvar; i++)
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
    if ((x + h - x2)*(x + h - x1) > 0.0) h = x2 - x;
    bsstep(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext);
    if ((x - x2)*(x2 - x1) >= 0.0) {
      for (i = 1; i <= nvar; i++) ystart[i] = y[i];
      free_Realvector(dydx, 1, nvar);
      free_Realvector(y, 1, nvar);
      free_Realvector(yscal, 1, nvar);
      return;
    }
    if (fabs(hnext) <= hmin) Die("Step size too small in odeint");
    h = hnext;
  }
  Die("Too many steps in routine odeint");
}

void IEKG::bsstep(Real y[], Real dydx[], int nv, Real *xx,
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
      if (xnew == (*xx)) 
	Die("step size underflow in bsstep in IEKG.cc");
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
void IEKG::IEKG_derivs(Real t, Real *x, Real *dxdt)
{
  const Real psi = x[1];
  const Real cpsi = cos(psi);
  const Real r = p/(1. + ecc*cpsi);

  const Real chi = x[2];
  const Real cchi = cos(chi);
  const Real z = zedminus*cchi*cchi;

  const Real sigma = Kerr::Sigma(r, a, z);
  const Real delta = Kerr::Delta(r, a);
  const Real r2pa2 = r*r + a*a;
  const Real gamma = E*(r2pa2*r2pa2/delta - a*a) - 2.*r*a*Lz/delta;
  const Real gampa2Ez = gamma + a*a*E*z;
  //
  const Real psitmp1 = (p - p3) - ecc*(p + p3*cpsi);
  const Real psitmp2 = (p - p4) + ecc*(p - p4*cpsi);
  dxdt[1] = sqrt((1. - E)*(1. + E)*psitmp1*psitmp2)/
    (gampa2Ez*(1. - ecc)*(1. + ecc));
  //
  dxdt[2] = sqrt(betazedplus - beta_wilkins*z)/gampa2Ez;
  //
  dxdt[3] = (Lz/(1. - z) + 2.*a*E*r/delta - a*a*Lz/delta)/gampa2Ez;
  //
  dxdt[4] = sigma/gampa2Ez;
}

void IEKG::mmid(Real y[], Real dydx[], int nvar, Real xs,
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
  IEKG_derivs(x, yn, yout);

  h2 = 2.0*h;
  for (n = 2; n <= nstep; n++) {
    for (i = 1; i <= nvar; i++) {
      swap = ym[i] + h2*yout[i];
      ym[i] = yn[i];
      yn[i] = swap;
    }
    x += h;
    IEKG_derivs(x, yn, yout);
  }
  for (i = 1; i<= nvar; i++)
    yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);
  free_Realvector(yn, 1, nvar);
  free_Realvector(ym, 1, nvar);
}

void IEKG::func_and_jac(Real *x, int n, Real *fvec, Real **fjac)
{
  const Real e = x[1];
  const Real lz = x[2];

  const Real tisqr = taniota*taniota;
  const Real cisqr = cosiota*cosiota;
  const Real q = lz*lz*tisqr;

  //
  // fvec[1] = [R(r_p) == 0]
  // fvec[2] = [R(r_a) == 0]
  //
  fvec[1] = gkg->RFunc(peri, a, 0, e, lz, q);
  fvec[2] = gkg->RFunc(ap, a, 0, e, lz, q);
  //
  // Now the jacobian: fjac[i][j] = dR(r_i)/dx_j
  // r_1 = r_p, r_2 = r_a
  // x_1 = e, x_2 = lz
  //
  const Real aelz = a*e - lz;
  fjac[1][1] = 2.*peri*(2.*a*aelz + peri*(a*a*e + peri*peri*e));
  fjac[1][2] = -2.*a*a*lz*tisqr + 2.*peri*(2.*(lz*tisqr - aelz)
					   - peri*lz/cisqr);
  fjac[2][1] = 2.*ap*(2.*a*aelz + ap*(a*a*e + ap*ap*e));
  fjac[2][2] = -2.*a*a*lz*tisqr + 2.*ap*(2.*(lz*tisqr - aelz)
					 - ap*lz/cisqr);
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
