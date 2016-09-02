#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

Real factrl(const int n)
{
  static int ntop=4;
  Real gammln(const Real xx);
  static Real a[33]={1.0, 1.0, 2.0, 6.0, 24.0};
  int j;
  
  if (n < 0) Die("Negative factorial in routine factrl");
  if (n > 32) return exp(gammln(n+1.0));
  while (ntop<n) {
    j=ntop++;
    a[ntop]=a[j]*ntop;
  }
  return a[n];
}
