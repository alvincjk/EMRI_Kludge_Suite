// NK: Fluxes

#ifndef _RRGW_H
#define _RRGW_H

#include <math.h>
#include "Globals.h"
#include "SWSH.h"

class RRGW {
public:
  void Flux_Infinity(const Real a, const int m, const Real lamb,
		     const Real w, const Real p, const Complex ZedH,
		     Real & Edot, Real & Lzdot);

  void Flux_Horizon(const Real a, const int m, const Real lamb,
		    const Real w, const Real p, const Complex ZedInf,
		    Real & Edot, Real & Lzdot);

  void Qdotrdot(const Real r, const Real a, const Real Q, const Real E,
		const Real Lz, const Real Edot, const Real Lzdot,
		Real & Qdot, Real & rdot);

  // Quasi-general purpose, but does not work well for inspirals!
  void Wave(const int m, const Real t_ret, const Real phi,
	    const Real S, const Real w, const Complex ZedH,
	    Real & hplus, Real & hcross);

  // Good for restricted inspirals, with indices m & k.
  void Wave(const int m, const int k, const Real N_m, const Real N_k,
	    const Real phi, const Real S, const Real w,
	    const Complex ZedH, Real & hplus, Real & hcross);

private:
  Real alpha_func(const Real a, const int m, const Real lamb,
		  const Real w, const Real p);
};
#endif
