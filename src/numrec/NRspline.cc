#include "Globals.h"
#include "NRUtil.h"

void spline(Real x[], Real y[], int n, Real yp1, Real ypn, Real y2[])
{
  int i, k;
  Real p, qn, sig, un, *u;

  u = Realvector(1, n - 1);
  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0.0;
  else {
    y2[1] = -0.5;
    u[1] = (3.0/(x[2] - x[1]))*((y[2] - y[1])/(x[2] - x[1]) - yp1);
  }
  for (i = 2; i <= n - 1; i++) {
    sig = (x[i] - x[i - 1])/(x[i + 1] - x[i - 1]);
    p = sig*y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0)/p;
    u[i] = (y[i + 1] - y[i])/(x[i + 1] - x[i]) -
      (y[i] - y[i - 1])/(x[i] - x[i - 1]);
    u[i] = (6.0*u[i]/(x[i + 1] - x[i - 1]) - sig*u[i - 1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n] - x[n - 1]))*(ypn - (y[n] - y[n - 1])/(x[n] - x[n - 1]));
  }
  y2[n]=(un - qn*u[n - 1])/(qn*y2[n - 1] + 1.0);
  for (k = n - 1; k >= 1; k--)
    y2[k] = y2[k]*y2[k + 1] + u[k];
  free_Realvector(u, 1, n - 1);
}

//
// Spline for complex functions of real arguments
//
void spline_complex(Real x[], Complex y[], int n, Complex yp1,
		    Complex ypn, Complex y2[])
{
  int i, k;
  Real qn, sig;
  Complex p, un, *u;

  u = Complexvector(1, n - 1);
  if (yp1.real() > 0.99e30)
    y2[1] = u[1] = 0.0;
  else {
    y2[1] = -0.5;
    u[1] = (3.0/(x[2] - x[1]))*((y[2] - y[1])/(x[2] - x[1]) - yp1);
  }
  for (i = 2; i <= n - 1; i++) {
    sig = (x[i] - x[i - 1])/(x[i + 1] - x[i - 1]);
    p = sig*y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0)/p;
    u[i] = (y[i + 1] - y[i])/(x[i + 1] - x[i]) -
      (y[i] - y[i - 1])/(x[i] - x[i - 1]);
    u[i] = (6.0*u[i]/(x[i + 1] - x[i - 1]) - sig*u[i - 1])/p;
  }
  if (ypn.real() > 0.99e30) {
    qn = 0.0;
    un = Complex(0.0, 0.0);
  }
  else {
    qn = 0.5;
    un = (3.0/(x[n] - x[n - 1]))*(ypn - (y[n] - y[n - 1])/(x[n] - x[n - 1]));
  }
  y2[n]=(un - qn*u[n - 1])/(qn*y2[n - 1] + 1.0);
  for (k = n - 1; k >= 1; k--)
    y2[k] = y2[k]*y2[k + 1] + u[k];
  free_Complexvector(u, 1, n - 1);
}
