#include <math.h>
#include "Globals.h"
#include "CKG.h"
#include "NRCKG.h"
#include "NRUtil.h"

#define EPSILON 1.e-9

CKG::CKG(const int emm, const int kay, const Real Lz_or_E_or_cosiota,
	 int Lz_or_E_or_cosiota_flag, const Real rad, const Real spin) :
  m(emm), k(kay), r(rad), a(spin)
{
  // This is needed right away ...
  Delta = Kerr::Delta(r, a);

  //
  // Note that Efunc() is correct only for outside the region of
  // horizon skimming orbits.
  //
  if (Lz_or_E_or_cosiota_flag == USING_Lz) {
    Lz = Lz_or_E_or_cosiota;
    E = Efunc();
    Q = Qfunc();
    Lovercosi=sqrt(Lz*Lz + Q);
    cosiota = Lz/Lovercosi;
  } else if (Lz_or_E_or_cosiota_flag == USING_E) {
    E = Lz_or_E_or_cosiota;
    Lz = Lzfunc();
    Q = Qfunc();
    Lovercosi=sqrt(Lz*Lz + Q);
    cosiota = Lz/Lovercosi;
  } else {
    cosiota = Lz_or_E_or_cosiota;
    //
    // Now, need to bisect this puppy until I get Lz and Q correct.
    //
    if (cosiota ==1.)
      E=Kerr::Eeqpro(r, a);
    else if (cosiota==-1.)
      E=Kerr::Eeqret(r, a);
    else if ((cosiota<EPSILON) && (cosiota > -1.*EPSILON)) {
      cosiota=0.;
      Lz=0.;
      E=Efunc();
      Q=Qfunc();
      //cout << cosiota << " " << Q << endl;
      Lovercosi=sqrt(Q);
    }
    else {
      Real Ehi, Elo, cosiota_guess;
      Elo = Kerr::Eeqpro(r, a);
      if (r >= Kerr::isco_ret(a))
	Ehi = Kerr::Eeqret(r, a);
      else
	Ehi = 1.;
      int BRACKETED = 0;
      while (!BRACKETED) {
	E = 0.5*(Ehi + Elo);
	Lz = Lzfunc();
	Q = Qfunc();
	cosiota_guess = Lz/sqrt(Lz*Lz + Q);
	if (cosiota_guess > cosiota) // Not energetic enough
	  Elo = E;
	else // too energetic
	  Ehi = E;
	if (fabs(cosiota_guess - cosiota) < EPSILON) BRACKETED = 1;
      }
    }
    Lz = Lzfunc();
    if (cosiota) {
      Q = Lz*Lz*(1. - cosiota*cosiota)/(cosiota*cosiota);
      Lovercosi=sqrt(Lz*Lz + Q);
    }
  }
  
  if (cosiota) {
    alpha_wilkins = Q + Lz*Lz;
    beta_wilkins = a*a*(1. - E)*(1. + E);
    const Real tmp = alpha_wilkins + beta_wilkins;
    
    if (beta_wilkins > EPSILON) {
      zedplus = (tmp + sqrt(tmp*tmp - 4.*Q*beta_wilkins));
      zedplus /= 2.*beta_wilkins;
      betazedplus = (tmp + sqrt(tmp*tmp - 4.*Q*beta_wilkins))/2.;
      zedminus = (tmp - sqrt(tmp*tmp - 4.*Q*beta_wilkins));
      zedminus /= 2.*beta_wilkins;
      if (fabs(zedminus) < 1.e-14) zedminus = 0.0;
      k_zmzp = sqrt(zedminus/zedplus);
    } else {
      zedplus = tmp/beta_wilkins - (Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
      betazedplus = tmp - (Q*beta_wilkins/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
      zedminus = (Q/tmp)*(1. - 2.*beta_wilkins/(tmp*tmp));
      k_zmzp = 0.;
    }
    thetamin = acos(sqrt(zedminus));
    thetamax = M_PI - thetamin;
    
    const Real r2pa2 = r*r + a*a;

    gamma = E*(r2pa2*r2pa2/Delta - a*a) + a*Lz*(1. - r2pa2/Delta);
    delta = a*E*(r2pa2/Delta - 1.) - a*a*Lz/Delta;

    phi_pi2 = phi_0(0.5*M_PI);
    t_pi2 = t_0(0.5*M_PI);

    //T_th = Period_Theta();
    //Phi_range = Phi();

    Om_th = Omega_Theta();
    Om_phi = Omega_Phi();
    wmk = omegamk();
    pmk = wmk - m*a/(2.*Kerr::rplus(a));
  }
}

Real CKG::Efunc()
{
  const Real rm1 = r - 1.;
  const Real numer = a*a*Lz*Lz*rm1 + r*Delta*Delta;

  //
  // This guy's a bit of a mess ...
  //
  Real denom = r*r*r*r*(r-3.);
  denom += a*a*r*(Lz*Lz + 2.*r*rm1);
  denom += a*a*a*a*(1.+r);
  denom  = Delta*sqrt(r*denom);
  denom += a*Lz*(r*r-a*a);

  return(numer/denom);
}

Real CKG::Lzfunc()
{
  const Real rm1 = r - 1.;
  const Real tmp1 = r*(1. + r*(E*E - 1.));
  const Real numer = E*(r*r - a*a) - Delta*sqrt(tmp1);
  Real denom = a*rm1;

  return(numer/denom);
}

Real CKG::Qfunc()
{
  const Real Q1 = SQR(E*(a*a + r*r) - a*Lz)/Delta;
  const Real Q2 = r*r + a*E*(a*E - 2.*Lz) + Lz*Lz;
  
  if (fabs(Q1 - Q2) < EPSILON) return(0.);
  else return(Q1 - Q2);
}

Real CKG::Omega_Theta()
{
  return(2.*M_PI/Period_Theta());
}

Real CKG::Phi()
{
  return(4.*phi_pi2);
}

Real CKG::Omega_Phi()
{
  return(Phi()/Period_Theta());
}

Real CKG::Period_Theta()
{
  return(4.*t_pi2);
}

//
// This function should never be used at this point.
// Retained for historical correctness.
//
// Real CKG::dtdz(const Real z)
// {
//   const Real tmp1 = gamma + a*a*E*z;
//   const Real tmp2 = z*beta_wilkins*(z - zedplus)*(z - zedminus);
//
//   return(tmp1/(2.*sqrt(tmp2)));
// }

Real CKG::dtdchi(const Real chi)
{
  const Real cchi = cos(chi);
  const Real z = zedminus*cchi*cchi;
  const Real tmp1 = gamma + a*a*E*z;
  const Real tmp2 = sqrt(betazedplus - beta_wilkins*z);

  return(tmp1/tmp2);
}

Real CKG::t(const Real chi)
{
  if (chi < 0.5*M_PI)
    return(t_pi2 - t_0(0.5*M_PI - chi));
  else
    return(t_pi2 + t_0(chi - 0.5*M_PI));
}

Real CKG::phi(const Real chi)
{
  if (chi < 0.5*M_PI)
    return(phi_pi2 - phi_0(0.5*M_PI - chi));
  else
    return(phi_pi2 + phi_0(chi - 0.5*M_PI));
}

Real CKG::t_0(const Real chi)
{
  const Real tmp1 = gamma/sqrt(betazedplus)*ellf(chi, k_zmzp);
  const Real tmp2 = (E/((1. - E)*(1. + E)))*sqrt(betazedplus)*
    (ellf(chi, k_zmzp) - elle(chi, k_zmzp));
  return(tmp1 + tmp2);
}

Real CKG::phi_0(const Real chi)
{
  const Real tmp = Lz*ellpi(chi, -zedminus, k_zmzp) +
    delta*ellf(chi, k_zmzp);

  return(tmp/sqrt(betazedplus));
}

Real CKG::omegamk()
{
  return(m*Om_phi + k*Om_th);
}
