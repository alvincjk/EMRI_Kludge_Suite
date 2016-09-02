// NK: Generic Kerr geodesics

#ifndef _GKG_H
#define _GKG_H

#include "Globals.h"

class GKG {
public:
  //
  // Sigma^2 * rdot^2 && its second derivative.
  //
  Real RFunc(const Real r, const Real a, const Real z,
	     const Real E, const Real Lz, const Real Q);
  Real dr_RFunc(const Real r, const Real a, const Real z,
		const Real E, const Real Lz, const Real Q);
  Real ddr_RFunc(const Real r, const Real a, const Real z,
		 const Real E, const Real Lz, const Real Q);
  //
  // Sigma^2 * thdot^2
  //
  Real ThetaFunc(const Real r, const Real a, const Real z,
		 const Real E, const Real Lz, const Real Q);
  //
  // Sigma phdot
  //
  Real PhiFunc(const Real r, const Real a, const Real z,
	       const Real E, const Real Lz, const Real Q);
  //
  // Sigma tdot
  //
  Real TFunc(const Real r, const Real a, const Real z,
	     const Real E, const Real Lz, const Real Q);
};

#endif
