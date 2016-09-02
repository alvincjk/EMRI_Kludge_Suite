#include "Globals.h"
#define MAXM 100

void zroots(Complex a[], int m, Complex roots[], int polish, Real EPS)
{
  void laguer(Complex a[], int m, Complex *x, int *its);
  int i, its, j, jj;
  Complex x, b, c, ad[MAXM];
  
  for (j = 0; j <= m; j++)
    ad[j] = a[j];
  for (j = m; j >= 1; j--) {
    x = Complex(0.0, 0.0);
    laguer(ad, j, &x, &its);
    if (fabs(x.imag()) <= 2.0*EPS*fabs(x.real()))
      x = Complex(x.real(), 0.0);
    roots[j] = x;
    b = ad[j];
    for (jj = j - 1; jj >= 0; jj--) {
      c = ad[jj];
      ad[jj] = b;
      b = c + x*b;
    }
  }
  if (polish)
    for (j = 1; j <= m; j++)
      laguer(a, m, &roots[j], &its);
  for (j = 2; j <= m; j++) {
    x = roots[j];
    for (i = j - 1; i >= 1; i--) {
      if (roots[i].real() <= x.real()) break;
      roots[i + 1] = roots[i];
    }
    roots[i + 1] = x;
  }
}
