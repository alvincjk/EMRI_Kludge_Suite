// NK: Integral

#ifndef _INTEGRAL_H
#define _INTEGRAL_H

#include "Globals.h"

//
// Prototype for ratint()
//
void ratint(const Real xa[], const Complex ya[], const int n,
	    const Real x, Complex *y, Complex *dy);
//
// The Integral class.  This will be inherited by all classes that
// use the integrator.
//
class Integral {
public:
  Integral(const Real epsilon) : EPSILON(epsilon) {};
  virtual Complex Integrand(const Real x)
    {return(Complex(0.,0.));};
  Complex qromb(const Real a, const Real b);

private:
  Real EPSILON;
  Complex trapzd(const Real a, const Real b, const int n);
  inline Complex FUNC(const Real x)
    {
      return(Integrand(x));
    };
};

#endif
