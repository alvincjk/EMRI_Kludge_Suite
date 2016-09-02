// NK: Inspiral

#ifndef _INSPIRAL_H
#define _INSPIRAL_H

#include <math.h>
#include "Globals.h"
#include "CKG.h"
#include "IEKG.h"
#include "NRUtil.h"
#include "GKG.h"
#include "NRLB.h"
#include "SWSH.h"
#include "GKTraj.h"
#include "GKInsp.h"

class Inspiral {
 public:
  Inspiral(const Real *params, double deltat, const int Nfrst, const int Nlast, bool qd, const bool datout, char *filename);
  ~Inspiral();

  const Real *p0;
  double dt;
  const int Nzero;
  int Nmx;
  int length;
  enum { ORBIT_P,    /* geodesic semilatus rectum / M */
	 LN_ORBIT_ETA,  /* logarithm of reduced mass ratio for nongeodesic orbits */
	 LN_SYSTEM_M,   /* log of central BH mass M (in solar mass units) */
	 ORBIT_E,    /* geodesic eccentricity */
	 ORBIT_PSI0, /* geodesic initial radial cyclic parameter (radians) */
	 ORBIT_I,    /* geodesic inclination */
	 ORBIT_CHI0, /* geodesic initial vertical cyclic parameter (radians) */
	 SYSTEM_A,      /* central BH spin / M */
	 NUM_PARAMS  /* number of search parameters */ };

  Real *x,*y,*z,*Quad[5];
  const bool genquad; /* Flag to tell code if it should generate the quadrupole moments of the inspiral as it integrates. This defaults to true. */ 
  const bool outputdata; /* Flag to tell code to output trajectory information. */
  char *filenm;
  Real p_end,e_end,i_end,zmin_end,pdot_end,edot_end,idot_end,zmindot_end;
  Real p_st,e_st,i_st,zmin_st,pdot_st,edot_st,idot_st,zmindot_st;

 private:
  Real *xx,*yy,*zz,**Qdmom;
};

#endif
