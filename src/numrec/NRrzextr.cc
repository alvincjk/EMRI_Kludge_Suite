#include "Globals.h"
#include "NRUtil.h"

void rzextr(int iest, Real xest, Complex yest[], Complex yz[],
	    Complex dy[], Complex **d, Real *x, int nv)
{
  int k, j;
  Real *fx;
  Complex yy, v, ddy, c, b1, b;
  
  fx = Realvector(1, iest);
  x[iest] = xest;
  if (iest == 1)
    for (j = 1; j <= nv; j++) {
      yz[j] = yest[j];
      d[j][1] = yest[j];
      dy[j] = yest[j];
    } else {
      for (k = 1; k < iest; k++)
	fx[k+1] = x[iest-k]/xest;
      for (j = 1; j <= nv; j++) {
	v = d[j][1];
	d[j][1] = yy = c = yest[j];
	for (k = 2; k <= iest; k++) {
	  b1 = fx[k]*v;
	  b = b1-c;
	  if (b != Complex(0., 0.)) {
	    b = (c-v)/b;
	    ddy = c*b;
	    c = b1*b;
	  } else
	    ddy = v;
	  if (k != iest) v = d[j][k];
	  d[j][k] = ddy;
	  yy += ddy;
	}
	dy[j] = ddy;
	yz[j] = yy;
      }
    }
  free_Realvector(fx, 1, iest);
}

void rzextr(int iest, Real xest, Real yest[], Real yz[], Real dy[],
	    Real **d, Real *x, int nv)
{
  int k, j;
  Real *fx;
  Real yy, v, ddy, c, b1, b;
  
  fx = Realvector(1, iest);
  x[iest] = xest;
  if (iest == 1)
    for (j = 1; j <= nv; j++) {
      yz[j] = yest[j];
      d[j][1] = yest[j];
      dy[j] = yest[j];
    } else {
      for (k = 1; k < iest; k++)
	fx[k+1] = x[iest-k]/xest;
      for (j = 1; j <= nv; j++) {
	v = d[j][1];
	d[j][1] = yy = c = yest[j];
	for (k = 2; k <= iest; k++) {
	  b1 = fx[k]*v;
	  b = b1-c;
	  if (b != 0.) {
	    b = (c-v)/b;
	    ddy = c*b;
	    c = b1*b;
	  } else
	    ddy = v;
	  if (k != iest) v = d[j][k];
	  d[j][k] = ddy;
	  yy += ddy;
	}
	dy[j] = ddy;
	yz[j] = yy;
      }
    }
  free_Realvector(fx, 1, iest);
}
