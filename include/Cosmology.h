// NK: Cosmology

#ifndef _COSMOLOGY_H
#define _COSMOLOGY_H

#include <math.h>
#include "Globals.h"
#include "NRUtil.h"

class Cosmology {
 public:
  Cosmology(Real hubble, Real OmegaM, Real OmegaL);
  ~Cosmology();

  Real H0,Om,OL;
  Real Ok,q0,sigma0;

  int curv;
  Real a0,rhocrit;

  Real ProperD(Real z);
  Real LumD(Real z);
  Real AngD(Real z);
  Real zofDp(Real d);
  Real zofDL(Real d);
  Real zofDA(Real d, int rt);

 private:
  Real uc,p,k,zeropt;
  Real zofMaxDA();
  void SolveCubic(const Real p0, const Real p1, const Real p2, Real *roots, bool &cmplex);
  void CubicSolver(const Real p0, const Real p1, const Real p2, Real *roots, bool &cmplex);
  void SolveQuartic(const Real p0, const Real p1, const Real p2, const Real p3, Real *roots, int &cmplex);
  Real raise(Real x, int n);
};

#endif
