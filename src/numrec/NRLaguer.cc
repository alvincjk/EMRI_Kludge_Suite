#include "Globals.h"
#include "NRUtil.h"
#define EPSS 1.0e-15
#define MR 15
#define MT 20
#define MAXIT (MT*MR)

void laguer(Complex a[], int m, Complex *x, int *its)
{
  int iter, j;
  Real abx, abp, abm, err;
  Complex dx, x1, b, d, f, g, h, sq, gp, gm, g2;
  static Real frac[MR + 1] =
  {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};

  for (iter = 1; iter <= MAXIT; iter++) {
    *its = iter;
    b = a[m];
    err = abs(b);
    d = f = Complex(0.0, 0.0);
    abx = abs(*x);
    for (j = m - 1; j >= 0; j--) {
      f = (*x)*f + d;
      d = (*x)*d + b;
      b = (*x)*b + a[j];
      err = abs(b) + abx*err;
    }
    err *= EPSS;
    if (abs(b) <= err) return;
    g = d/b;
    g2 = g*g;
    h = g2 - (2.0 * f/b);
    sq = sqrt(((Real)(m - 1))*(((Real)m)*h - g2));
    gp = g + sq;
    gm = g - sq;
    abp = abs(gp);
    abm = abs(gm);
    if (abp < abm) gp = gm;
    dx = (FMAX(abp, abm) > 0.0 ?
	  Complex((float) m, 0.0)/gp :
	  exp(log(1 + abx)) * Complex(cos((float)iter), sin((float)iter)));
    x1 = (*x) - dx;
    if (x->real() == x1.real() && x->imag() == x1.imag())
      return;
    if (iter % MT)
      *x = x1;
    else
      *x = *x - frac[iter/MT]*dx;
  }
  Die("too many iterations in laguer");
  return;
}
#undef EPSS
#undef MR
#undef MT
#undef MAXIT
#undef NRANSI
