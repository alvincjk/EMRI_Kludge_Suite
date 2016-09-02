#include <math.h>
#include "Globals.h"
#include "CKR.h"
#include "Integral.h"
#include "SWSH.h"
#include "CKG.h"
#include "SNT.h"

// Constructor.
CKR::CKR(SWSH *swsh_in, SNT *snt_in, CKG *ckg_in, const Real epsilon) :
  Integral(epsilon), swsh(swsh_in), snt(snt_in), ckg(ckg_in)
{
  //
  // Local storage of orbital characteristics.
  //
  r = ckg->r;
  a = ckg->a;
  E = ckg->E;
  Lz = ckg->Lz;
  Q = ckg->Q;
  //
  // Local storage of mode characteristics.
  //
  l = swsh->Gl;
  m = ckg->m;
  k = ckg->k;
  wmk = ckg->wmk;
  pmk = ckg->pmk;
  //
  // Local storage of useful quantities that appear all over the
  // place.
  //
  Delta = Kerr::Delta(r, a);
  dr_Delta = Kerr::dr_Delta(r);
  SN_Kay = snt->K(r);
  dr_SN_Kay = snt->dr_K(r);
  //
  // Compute amplitudes ZedH, ZedInf.
  //
  ZedH = Zed(snt->TeukRH, snt->dr_TeukRH, snt->ddr_TeukRH);
  ZedInf = Zed(snt->TeukRInf, snt->dr_TeukRInf, snt->ddr_TeukRInf);
}

// Overloaded version of constructor: gives option of skipping calc of
// flux down horizon.
CKR::CKR(SWSH *swsh_in, SNT *snt_in, CKG *ckg_in, const Real epsilon,
	 const int usehorz) :
  Integral(epsilon), swsh(swsh_in), snt(snt_in), ckg(ckg_in),
  UseHoriz(usehorz)
{
  //
  // Local storage of orbital characteristics.
  //
  r = ckg->r;
  a = ckg->a;
  E = ckg->E;
  Lz = ckg->Lz;
  Q = ckg->Q;
  //
  // Local storage of mode characteristics.
  //
  l = swsh->Gl;
  m = ckg->m;
  k = ckg->k;
  wmk = ckg->wmk;
  pmk = ckg->pmk;
  //
  // Local storage of useful quantities that appear all over the
  // place.
  //
  Delta = Kerr::Delta(r, a);
  dr_Delta = Kerr::dr_Delta(r);
  SN_Kay = snt->K(r);
  dr_SN_Kay = snt->dr_K(r);
  //
  // Compute amplitudes ZedH, ZedInf.
  //
  ZedH = Zed(snt->TeukRH, snt->dr_TeukRH, snt->ddr_TeukRH);
  if (UseHoriz)
    ZedInf = Zed(snt->TeukRInf, snt->dr_TeukRInf, snt->ddr_TeukRInf);
  else
    ZedInf = 0.0;
}

// The complex amplitudes ZH, ZInf that are used to construct
// Edot, Lzdot, and the waveform.  To get ZH, pass in RH and
// its derivatives; likewise for ZInf.
Complex CKR::Zed(const Complex R,
		   const Complex dr_R,
		   const Complex ddr_R)
{
  tmp_TeukR = R;
  tmp_dr_TeukR = dr_R;
  tmp_ddr_TeukR = ddr_R;

  const Complex I = qromb(0.0, M_PI);
  const Complex denom = II*wmk*snt->Bin*ckg->T_th;

  return(M_PI*I/denom);
}

// The phase.
Real CKR::Psi_mk(const Real chi)
{
  return(wmk*ckg->t(chi) - m*ckg->phi(chi));
}

// The Newman-Penrose quantity rho (with sign as in
// Teukolsky's papers).
Complex CKR::rho_func(const int n, const Real z)
{
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Complex tmp = Complex(r, -a*ct);

  return(-1./tmp);
}

// The C_{ab} functions are projections of the orbiting
// particle's stress-energy tensor onto the Newman-Penrose
// tetrad (modulo some factors).
Complex CKR::Cnn(const int seg, const int n, const Real z)
{
  const Real tmp = E*(r*r + a*a) - a*Lz;
  const Real sig = Kerr::Sigma(r, a, z);

  return(0.25*tmp*tmp/(sig*sig*gkg.TFunc(r, a, z, E, Lz, Q)));
}

Complex CKR::Cnmbar(const int seg, const int n, const Real z)
{
  const Real st = sqrt(1.-z);
  const Complex rho = rho_func(n, z);
  const Real sig = Kerr::Sigma(r, a, z);

  const Real tmp1 = E*(r*r + a*a) - a*Lz;
  Complex tmp2;

  if(seg == 0)
    tmp2 = Complex(sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), st*a*E - Lz/st);
  else
    tmp2 = Complex(-sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), st*a*E - Lz/st);

  return(rho*tmp1*tmp2/(2.*sqrt(2.)*sig*gkg.TFunc(r, a, z, E, Lz, Q)));
}

Complex CKR::Cmbarmbar(const int seg, const int n, const Real z)
{
  const Real st = sqrt(1. - z);
  const Complex rho = rho_func(n, z);

  Complex tmp;

  if (seg == 0)
    tmp = Complex(sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), a*E*st - Lz/st);
  else
    tmp = Complex(-sqrt(gkg.ThetaFunc(r, a, z, E, Lz, Q)), a*E*st - Lz/st);

  return(rho*rho*tmp*tmp/(2.*gkg.TFunc(r, a, z, E, Lz, Q)));
}

// The A_{abj}'s are functions that appear when the source is
// integrated over the Green's function.  See notes for details
// --- it's a touch involved.
Complex CKR::Ann0(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real st = sqrt(1. - z);

  const Complex tmp1 = 2.*II*a*rho*st*swsh->l2dagspheroid(ct);
  const Real tmp2 = swsh->l1dagl2dagspheroid(ct);

  const Complex pref = -2.*Cnn(seg, n, z)/(Delta*Delta*rho*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CKR::Anmbar0(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real st = sqrt(1. - z);
  const Real S = swsh->spheroid(ct);
  const Real L2dagS = swsh->l2dagspheroid(ct);

  const Complex tmp1 = (II*SN_Kay/Delta - rho - rhob)*L2dagS;
  const Complex tmp2 = (II*SN_Kay/Delta + rho +
			rhob)*II*a*st*S*(rho - rhob);

  const Complex pref = -2.*sqrt(2.)*Cnmbar(seg, n, z)/(Delta*rho3);

  return(pref*(tmp1 + tmp2));
}

Complex CKR::Anmbar1(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real st = sqrt(1. - z);
  const Real S = swsh->spheroid(ct);
  const Real L2dagS = swsh->l2dagspheroid(ct);

  const Complex tmp = L2dagS + II*a*st*(rho - rhob)*S;

  const Complex pref = -2.*sqrt(2.)*Cnmbar(seg, n, z)/(Delta*rho3);
  return(pref*tmp);
}

Complex CKR::Ambarmbar0(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real S = swsh->spheroid(ct);

  const Complex tmp1 = (SN_Kay/Delta)*(SN_Kay/Delta);
  const Complex tmp2 = 2.*rho*II*SN_Kay/Delta;
  const Complex tmp3 = II*(dr_SN_Kay - SN_Kay*dr_Delta/Delta)/Delta;

  const Complex pref = S*Cmbarmbar(seg, n, z)*rhob/rho3;

  return(pref*(tmp1 + tmp2 + tmp3));
}

Complex CKR::Ambarmbar1(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real S = swsh->spheroid(ct);

  const Complex tmp = -II*SN_Kay/Delta;

  const Complex pref = 2.*S*rhob*Cmbarmbar(seg, n, z)/rho3;

  return(pref*(rho + tmp));
}

Complex CKR::Ambarmbar2(const int seg, const int n, const Real z)
{
  const Complex rho = rho_func(n, z);
  const Complex rho3 = rho*rho*rho;
  const Complex rhob = conj(rho);
  const Real ct = (n == 0 ? sqrt(z) : -sqrt(z));
  const Real S = swsh->spheroid(ct);

  return(-Cmbarmbar(seg, n, z)*rhob*S/rho3);
}

// "dZed" is what you get when you do the radial integral of
// source and Green's function.  Integrate it and you get
// Zed [modulo factors which are handled in the function
// Zed().]
Complex CKR::dZed(const int seg, const int n, const Real chi)
{
  const Real cchi = cos(chi);
  const Real z = ckg->zedminus*cchi*cchi;
  const Real sign = (seg == 0) ? 1. : -1.;

  const Complex term0 = tmp_TeukR*(Ann0(seg, n, z) +
				   Anmbar0(seg, n, z) +
				   Ambarmbar0(seg, n, z));
  const Complex term1 = -tmp_dr_TeukR*(Anmbar1(seg, n, z) +
				       Ambarmbar1(seg, n, z));
  const Complex term2 = tmp_ddr_TeukR*Ambarmbar2(seg, n, z);

  const Complex pref = ckg->dtdchi(chi)*exp(sign*II*Psi_mk(chi));

  return(pref*(term0 + term1 + term2));
}

// A wrapper so that the integrator easily integrates dZed
// over theta to get Z.
Complex CKR::Integrand(const Real chi)
{
  const int n = (chi < 0.5*M_PI ? 0 : 1);

  return(dZed(0, n, chi) + dZed(1, n, chi));
}
