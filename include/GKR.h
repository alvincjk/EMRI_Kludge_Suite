// NK: Generic Kerr radiation reaction

#ifndef _GKR_H
#define _GKR_H

#include <math.h>
#include "Globals.h"
#include "RRGW.h"

class GKR {
public:
  GKR(const Real semilatus, const Real ecc, const Real energy,
      const Real angmom, const Real carter, const Real spin);
  Real p_LSO(Real csi);

  Real Edot, Lzdot, Qdot;
  Real radot, rpdot;
  Real pdot, edot, idot,cosidot,Ldotcirc,Qdotcirc;

private:
  Real Lzdot_GG();
  Real Qdot_GG();
  Real Edot_GG();
 
  Real Denom(const Real r);
  Real Numer1(const Real r);
  Real Numer2(const Real r);
  Real Numer3(const Real r);

  RRGW *rrgw;

  Real p, e, E, Lz, Q, kerr_a, Ecirc, Lzcirc, Qcirc, Lovercosicirc;
  Real cosiota, siniotasqr, sini,ra, rp, Lovercosi, a_sm, ome, ope, ome2;  
};
#endif
