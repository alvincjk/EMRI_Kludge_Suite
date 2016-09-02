// NK: Lattice Boltzmann

#ifndef _LB_H
#define _LB_H

#include "Globals.h"
#include "GKG.h"

class LB {
public:
  LB(const Real arr, const Real ayy) : r(arr), a(ayy) {
    gkg = new GKG;
  };
  ~LB() {
    delete gkg;
  };
  void mnewt(const int ntrial, Real x[], const int n, const Real tolx,
	     const Real tolf);
  GKG *gkg;

private:
  void func_and_jac(Real *x, int n, Real *fvec, Real **fjac);

  Real r, a;
};

#endif
