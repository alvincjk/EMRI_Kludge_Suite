// NK: Circular Kerr radiation

#ifndef _CKR_H
#define _CKR_H

#include "Globals.h"
#include "Integral.h"
#include "SWSH.h"
#include "GKG.h"
#include "CKG.h"
#include "SNT.h"

class CKR : public Integral {
public:
  CKR(SWSH *swsh_in, SNT *snt_in, CKG *ckg_in, const Real epsilon);
  CKR(SWSH *swsh_in, SNT *snt_in, CKG *ckg_in, const Real epsilon,
      const int usehorz);

  Complex dZed(const int seg, const int n, const Real chi);

  Complex ZedH;
  Complex ZedInf;

private:
  Complex rho_func(const int n, const Real z);
  Complex Cnn(const int seg, const int n, const Real z);
  Complex Cnmbar(const int seg, const int n, const Real z);
  Complex Cmbarmbar(const int seg, const int n, const Real z);
  Complex Ann0(const int seg, const int n, const Real z);
  Complex Anmbar0(const int seg, const int n, const Real z);
  Complex Ambarmbar0(const int seg, const int n, const Real z);
  Complex Anmbar1(const int seg, const int n, const Real z);
  Complex Ambarmbar1(const int seg, const int n, const Real z);
  Complex Ambarmbar2(const int seg, const int n, const Real z);

  Real Psi_mk(const Real chi);

  // The above methods are used to put together ...

  // A wrapper for dZed, for the integrator.
  Complex Integrand(const Real chi);

  // What you get after dZed is integrated and summed.
  Complex Zed(const Complex R,
	      const Complex dr_R,
	      const Complex ddr_R);

  // Other stuff we need in lots of places.
  int l, m, k;
  SWSH *swsh;
  SNT *snt;
  GKG gkg;
  CKG *ckg;
  Real r, a, E, Lz, Q;
  Real wmk, pmk;
  Real Delta, dr_Delta;
  Real SN_Kay, dr_SN_Kay;

  // Used by overloaded version to possibly skip calc of the horizon
  // flux.
  int UseHoriz;

  // Temporary storage of Teukolsky functions.
  Complex tmp_TeukR, tmp_dr_TeukR, tmp_ddr_TeukR;
};

#endif
