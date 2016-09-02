#include <math.h>
#include "Globals.h"
#include "RRGW.h"
#include "SWSH.h"

// Calculates Edot, Lzdot at infinity.
void RRGW::Flux_Infinity(const Real a, const int m, const Real lamb,
			 const Real w, const Real p,
			 const Complex ZedH,
			 Real & Edot, Real & Lzdot)
{
  const Real zed = abs(ZedH);

  Edot = zed*zed/(4.*M_PI*w*w);
  Lzdot = m*zed*zed/(4.*M_PI*w*w*w);
}

// Calculates Edot, Lzdot down the horizon.
void RRGW::Flux_Horizon(const Real a, const int m, const Real lamb,
			const Real w, const Real p,
			const Complex ZedInf,
			Real & Edot, Real & Lzdot)
{
  const Real zed = abs(ZedInf);
  const Real alpha = alpha_func(a, m, lamb, w, p);

  Edot = zed*zed*alpha/(4.*M_PI*w*w);
  Lzdot = m*zed*zed*alpha/(4.*M_PI*w*w*w);
}

// Calculates Qdot and rdot from Edot, Lzdot.
// Kind of messy ... check this carefully.
void RRGW::Qdotrdot(const Real r, const Real a, const Real Q,
		    const Real E, const Real Lz,
		    const Real Edot, const Real Lzdot,
		    Real & Qdot, Real & rdot)
{
  //
  // Call this routine with Edot, Lzdot in radiation.
  //
  const Real Edot_particle = -Edot;
  const Real Lzdot_particle = -Lzdot;

  const Real lz2pq = Lz*Lz + Q;
  const Real ome2 = 1. - E*E;
  const Real lzmae = Lz - a*E;
  const Real qplzmae2 = Q + lzmae*lzmae;

  const Real c11_6 = -4.*E*ome2;
  const Real c11_5 = 12.*E;
  const Real c11_4 = -2.*E*(a*a*ome2 + 3.*lz2pq);
  const Real c11_3 = 8.*(a*a*E*(1. + ome2) + E*lz2pq - 2.*a*Lz);
  const Real c11_2 = -2.*a*(a*E*(6. + lz2pq + a*a*ome2) - 6.*Lz);
  const Real c11_1 = 4.*a*a*E*qplzmae2;
  const Real c11_0 = -4.*a*lzmae*qplzmae2;
  const Real c11 = c11_0 + r*(c11_1 +
			      r*(c11_2 +
				 r*(c11_3 +
				    r*(c11_4 +
				       r*(c11_5 +
					  r*c11_6)))));
  
  const Real c12_4 = -4.*Lz*ome2;
  const Real c12_3 = 16.*ome2*lzmae;
  const Real c12_2 = 2.*(Lz*(a*a*ome2 + lz2pq) - 6.*lzmae);
  const Real c12_1 = -4.*Lz*qplzmae2;
  const Real c12_0 = 4.*lzmae*qplzmae2;
  const Real c12 = c12_0 + r*(c12_1 +
			      r*(c12_2 +
				 r*(c12_3 +
				    r*c12_4)));

  const Real c21_5 = 2.*E;
  const Real c21_4 = -6.*E;
  const Real c21_3 = 4.*a*a*E;
  const Real c21_2 = 2.*a*(lzmae - a*E);
  const Real c21_1 = 2.*a*a*a*a*E;
  const Real c21_0 = -2.*a*a*a*lzmae;
  const Real c21 = c21_0 + r*(c21_1 +
			      r*(c21_2 +
				 r*(c21_3 +
				    r*(c21_4 + 
				       r*c21_5))));
  
  const Real c22 = 2.*a*a*lzmae + r*(-2.*a*a*Lz +
				     r*2.*a*E);

  const Real d_4 = -2*ome2;
  const Real d_3 = 8.*ome2;
  const Real d_2 = lz2pq - 6. - 5.*a*a*ome2;
  const Real d_1 = 2*(a*a*(2. + ome2) + 2.*a*E*Lz - lz2pq);
  const Real d_0 = 2.*lz2pq +
    a*(-4*E*Lz + a*(2.*E*E - lz2pq - a*a*ome2));
  const Real d = d_0 + r*(d_1 +
			  r*(d_2 +
			     r*(d_3 +
				r*d_4)));
  
  Qdot = -(c11*Edot_particle + c12*Lzdot_particle)/d;
  rdot = -(c21*Edot_particle + c22*Lzdot_particle)/d;
}

// Calculates (r times) hplus, hcross.  To speed things up,
// calculate S(ct) in the calling routine.
// This is general, but doesn't really do inspirals well.
void RRGW::Wave(const int m, const Real t_ret,
		const Real phi, const Real S,
		const Real w, const Complex ZedH,
		Real & hplus, Real & hcross)
{
  const Complex h = -2*S*ZedH*exp(II*(m*phi - w*t_ret))/(w*w);

  hplus  = h.real();  hcross = -h.imag();
}

// Calculates (r times) hplus, hcross.  To speed things up,
// calculate S(ct) in the calling routine.
// This does restricted inspirals.
void RRGW::Wave(const int m, const int k,
		const Real N_m, const Real N_k,
		const Real phi, const Real S,
		const Real w, const Complex ZedH,
		Real & hplus, Real & hcross)
{
  const Complex h = -2*S*ZedH*exp(II*(m*phi - 2.*M_PI*(m*N_m + k*N_k)))/(w*w);

  hplus  = h.real();  hcross = -h.imag();
}

// Factor that appears in the fluxes on the horizon; arises from
// the conversion of the Newman-Penrose tetrad to a Hartle-
// Hawking tetrad (used to describe how the horizon grows due to
// stuff falling in).  See Teukolsky and Press 1974.
Real RRGW::alpha_func(const Real a, const int m, const Real lamb,
		      const Real w, const Real p)
{
  const Real r_plus = Kerr::rplus(a);
  const Real eps = sqrt(1. - a*a)/(4.*r_plus);
  const Real lp2 = lamb + 2.;
  const Real aw_m_m = a*w - ((Real)m);

  const Real tmp1 = lp2*lp2 - 4.*a*w*aw_m_m;
  const Real tmp2 = lamb*lamb - 36.*a*w*aw_m_m;
  const Real tmp3 = 48.*a*w*(2.*lamb + 3.)*(2.*a*w - ((Real)m));
  const Real tmp4 = 256.*pow(2.*r_plus, 5.)*p*(p*p + 4.*eps*eps)*
    (p*p + 16.*eps*eps)*w*w*w;

  const Real abssqrClm = tmp1*tmp2 + tmp3 + 144.*w*w*(1. - a*a);

  return(tmp4/abssqrClm);
}

