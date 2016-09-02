// NK: White dwarf binaries

#ifndef _WDBINARY_H
#define _WDBINARY_H

#include <math.h>
#include "Globals.h"


class WDBinary {
 public:
  Real f,inc,M,eta,csthl,pL,csthsky,pS;
  Real a,ainM,dphys,forb;
  Real cphis,sphis,sinthl,sinthsky,d,phi0,cphil,sphil;
  WDBinary(const Real freq, const Real totmass, const Real massrat, const Real csthetaL, const Real phiL, const Real csthetaS, const Real phiS, 
		     const Real dist, const Real phase);
  ~WDBinary();
  void hpluscross(Real *hp, Real *hx, Real dt, int npts);
  void hLISA(Real *hI, Real *hII, Real dt, int npts, Real t0);

 private:
  Real amp,cosi;
};

#endif
