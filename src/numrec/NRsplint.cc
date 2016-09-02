#include "Globals.h"
#include "NRUtil.h"

void splint(Real xa[], Real ya[], Real y2a[], int n, Real x, Real *y)
{
  int klo, khi, k;
  Real h, b, a;
  
  klo = 1; khi = n;
  while (khi - klo > 1) {
    k = (khi + klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0) Die("Bad xa input to routine splint");
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;
  *y = a*ya[klo] + b*ya[khi] +
    ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}

//
// Splint for complex functions of real arguments
//
void splint_complex(Real xa[], Complex ya[], Complex y2a[], int n,
		    Real x, Complex *y)
{
  int klo, khi, k;
  Real h, b, a;
  
  klo = 1; khi = n;
  while (khi - klo > 1) {
    k = (khi + klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0) Die("Bad xa input to routine splint_complex");
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;
  *y = a*ya[klo] + b*ya[khi] +
    ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}
