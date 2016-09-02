// NK: Globals

#ifndef _Globals_H
#define _Globals_H

#include <iostream>
#include <math.h>
#include <complex>

using namespace std;

// For convenience, can change the floating-point type with a
// single preprocessor symbol
#define USE_FLOAT       0
#define USE_DOUBLE      1
#define USE_LONG_DOUBLE 2

#if !defined(REAL_TYPE) || REAL_TYPE == USE_DOUBLE
typedef double Real;
#elif REAL_TYPE == USE_FLOAT
typedef float Real;
#elif REAL_TYPE == USE_LONG_DOUBLE
typedef long double Real;
#else
Illegal value!
#endif

#define MSUN_CM (147661.2814476609)
#define C (299792458.)
#define MSUN_SEC (MSUN_CM/(100.*C))
#define G (6.6726e-11)
#define YEAR (31536000.) /* One year in seconds. */
#define Msun (1.9889e30)
#define SOLARMASSINSEC (G*Msun/(C*C*C))
#define AU (1.4959787066e11)
#define pc (3.08567818585e16)
#define DIM_BASIS (5) /* Define the dimension of the basis of intrinsic waveforms - 5 for pure quadrupole. */

typedef complex<Real> Complex;
const Complex II = Complex(0.,1.);

//
// Container for some useful global functions.
//
class Kerr {
 public:
  //
  // r_+ = M + \sqrt{M^2 + a^2}; M \equiv 1
  //
  static Real rplus(const Real a)
    {
      return(1. + sqrt((1. - a)*(1. + a)));
    };
  //
  // r_- = M - \sqrt{M^2 - a^2}
  //
  static Real rminus(const Real a)
    {
      return(1. - sqrt((1. - a)*(1. + a)));
    };
  //
  // r^* = r + 2 r_+/(r_+ - r_-) \ln{(r - r_+)/2M}
  //         - 2 r_-/(r_+ - r_-) \ln{(r - r_-)/2M}
  //
  static Real rstar(const Real r, const Real a)
    {
      const Real rm = rminus(a);
      const Real rp = rplus(a);
      
      return(r + ((2.*rp)/(rp - rm))*log((r - rp)/2.) -
	     ((2.*rm)/(rp - rm))*log((r - rm)/2.));
    };
  //
  // \Delta  = r^2 - 2 M r + a^2
  //
  static Real Delta(const Real r, const Real a)
    {
      return(r*r - 2.*r + a*a);
    }
  //
  // d\Delta/dr  = 2 r - 2 M
  //
  static Real dr_Delta(const Real r)
    {
      return(2.*(r - 1.));
    }
  //
  // All right, so this one's kind of silly ...
  //
  // d^2\Delta/dr^2  = 2
  //
  static Real ddr_Delta()
    {
      return(2.);
    }
  //
  // \Sigma = r^2 + a^2 \cos^2\theta
  //
  static Real Sigma(const Real r, const Real a, const Real z)
    {
      return(r*r + a*a*z);
    };
  //
  // d\Sigma/dr = 2 r
  //
  static Real dr_Sigma(const Real r)
    {
      return(2.*r);
    };
  //
  // All right, so this one's kind of silly too ...
  //
  // d^2\Sigma/dr^2 = 2
  //
  static Real ddr_Sigma()
    {
      return(2.);
    };

  static Real Eeqpro(const Real r, const Real a)
    {
      const Real v = 1./sqrt(r);

      const Real numer = 1. - v*v*(2. - a*v);
      const Real denom = 1. - v*v*(3. - 2.*a*v);

      return(numer/sqrt(denom));
    };

  static Real Eeqret(const Real r, const Real a)
    {
      const Real v = 1./sqrt(r);
      
      const Real numer = 1. - v*v*(2. + a*v);
      const Real denom = 1. - v*v*(3. + 2.*a*v);
      
      return(numer/sqrt(denom));
    };
  
  static Real Lzeqpro(const Real r, const Real a)
    {
      const Real v = 1./sqrt(r);
      
      const Real numer = 1. - a*v*v*v*(2. - a*v);
      const Real denom = 1. - v*v*(3. - 2.*a*v);
      
      return(r*v*numer/sqrt(denom));
    };

  static Real Lzeqret(const Real r, const Real a)
    {
      const Real v = 1./sqrt(r);
      
      const Real numer = 1. + a*v*v*v*(2. + a*v);
      const Real denom = 1. - v*v*(3. + 2.*a*v);
      
      return(-r*v*numer/sqrt(denom));
    };

  static Real Omega_phi_eqpro(const Real r, const Real a)
    {
      return(1./(sqrt(r*r*r) + a));
    };

  static Real Omega_phi_eqret(const Real r, const Real a)
    {
      return(-1./(sqrt(r*r*r) - a));
    };

  static Real isco_pro(const Real a)
    {
      const Real Z1 = 1. + (pow(1. + a, 1./3.) + pow(1. - a, 1./3.))*
	pow((1. + a)*(1. - a), 1./3.);
      const Real Z2 = sqrt(3.*a*a + Z1*Z1);

      return(3. + Z2 - sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)));
    };

  static Real isco_ret(const Real a)
    {
      const Real Z1 = 1. + (pow(1. + a, 1./3.) + pow(1. - a, 1./3.))*
	pow((1. + a)*(1. - a), 1./3.);
      const Real Z2 = sqrt(3.*a*a + Z1*Z1);

      return(3. + Z2 + sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)));
    };
};
//
// Holder for the trajectory data.  Each element of such a structure
// holds one point along the parameter space trajectory.
//
struct TrajData {
  Real t; // This is in radiation reaction time!
  Real p, ecc, cosiota;
  Real E, Lz, Q;
  Real p3, p4;
  Real beta_wilkins, zedminus, betazedplus;
};
//
// Index range
//
#define L_MAX 100
#define K_MAX 100

#endif
