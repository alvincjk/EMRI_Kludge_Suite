#include "Globals.h"
#include "NRUtil.h"

#define MAXSTP 10000
#define TINY 1.0e-30
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void odeint(Real ystart[], int nvar, Real x1, Real x2, Real eps,
	    Real h1, Real hmin,
	    void (*derivs)(Real, Real [], Real []),
	    void (*rkqs)(Real [], Real [], int, Real *, Real, Real,
			 Real [], Real *, Real *,
			 void (*)(Real, Real [], Real [])))
{
  int nstp, i;
  Real xsav, x, hnext, hdid, h;
  Real *yscal, *y, *dydx;
  
  yscal = Realvector(1, nvar);
  y = Realvector(1, nvar);
  dydx = Realvector(1, nvar);
  x = x1;
  h = SIGN(h1, x2 - x1);
  for (i = 1; i <= nvar; i++) y[i] = ystart[i];
  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    (*derivs)(x, y, dydx);
    for (i = 1; i <= nvar; i++)
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
    if ((x + h - x2)*(x + h - x1) > 0.0) h = x2 - x;
    (*rkqs)(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if ((x - x2)*(x2 - x1) >=  0.0) {
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

void rkqs(Real y[], Real dydx[], int n, Real *x, Real htry, Real eps,
	  Real yscal[], Real *hdid, Real *hnext,
	  void (*derivs)(Real, Real [], Real []))
{
  void rkck(Real y[], Real dydx[], int n, Real x, Real h,
	    Real yout[], Real yerr[],
	    void (*derivs)(Real, Real [], Real []));
  int i;
  Real errmax, h, xnew, *yerr, *ytemp;
  
  yerr = Realvector(1, n);
  ytemp = Realvector(1, n);
  h = htry;
  for (;;) {
    rkck(y, dydx, n, *x, h, ytemp, yerr, derivs);
    errmax = 0.0;
    for (i = 1; i <= n; i++) {
      errmax = FMAX(errmax, fabs(yerr[i]/yscal[i]));
    }
    errmax /= eps;
    if (errmax > 1.0) {
      h = SAFETY*h*pow(errmax, PSHRNK);
      if (h < 0.1*h) h *= 0.1;
      xnew = (*x) + h;
      if (xnew == *x) Die("stepsize underflow in rkqs");
      continue;
    } else {
      if (errmax > ERRCON) *hnext = SAFETY*h*pow(errmax, PGROW);
      else *hnext = 5.0*h;
      *x += (*hdid = h);
      for (i = 1; i <= n; i++) y[i] = ytemp[i];
      break;
    }
  }
  free_Realvector(ytemp, 1, n);
  free_Realvector(yerr, 1, n);
}

void rkck(Real y[], Real dydx[], int n, Real x, Real h, Real yout[],
	  Real yerr[], void (*derivs)(Real, Real [], Real []))
{
  int i;
  static Real a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
    b21 = 0.2, b31 = 3.0/40.0, b32 = 9.0/40.0, b41 = 0.3, b42 = -0.9,
    b43 = 1.2, b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0,
    b54 = 35.0/27.0, b61 = 1631.0/55296.0, b62 = 175.0/512.0,
    b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
    c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    dc5 = -277.0/14336.0;
  Real dc1 = c1 - 2825.0/27648.0, dc3 = c3 - 18575.0/48384.0,
    dc4 = c4 - 13525.0/55296.0, dc6 = c6 - 0.25;
  Real *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;
  
  ak2 = Realvector(1, n);
  ak3 = Realvector(1, n);
  ak4 = Realvector(1, n);
  ak5 = Realvector(1, n);
  ak6 = Realvector(1, n);
  ytemp = Realvector(1, n);
  for (i = 1; i <= n; i++)
    ytemp[i] = y[i] + b21*h*dydx[i];
  (*derivs)(x + a2*h, ytemp, ak2);
  for (i = 1; i <= n; i++)
    ytemp[i] = y[i] + h*(b31*dydx[i] + b32*ak2[i]);
  (*derivs)(x + a3*h, ytemp, ak3);
  for (i = 1; i <= n; i++)
    ytemp[i] = y[i] + h*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]);
  (*derivs)(x + a4*h, ytemp, ak4);
  for (i = 1; i <= n; i++)
    ytemp[i] = y[i] + h*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);
  (*derivs)(x + a5*h, ytemp, ak5);
  for (i = 1; i <= n; i++)
    ytemp[i] = y[i] + h*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i]
			 + b65*ak5[i]);
  (*derivs)(x + a6*h, ytemp, ak6);
  for (i = 1; i <= n; i++)
    yout[i] = y[i] + h*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]);
  for (i = 1; i <= n; i++) {
    yerr[i] = h*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i]
		 + dc6*ak6[i]);
  }
  free_Realvector(ytemp, 1, n);
  free_Realvector(ak6, 1, n);
  free_Realvector(ak5, 1, n);
  free_Realvector(ak4, 1, n);
  free_Realvector(ak3, 1, n);
  free_Realvector(ak2, 1, n);
}

#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI
#undef MAXSTP
#undef TINY
