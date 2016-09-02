// Methods for generic Kerr radiation reaction.
//
// Jonathan Gair, 10 August 2005
//
// The trajectory integration scheme used here is described in Gair & Glampedakis 2005 (gr-qc/0510129, hereafter [1]). 
// We prescribe dLzdt and dQdt, plus the eccentricity pieces of dEdt. The circular part of dEdt is
// determined by the "circular to circular" condition described in [1]. The circular part of dLzdt is
// determined from a fit to Teukolsky based computations supplied by Scott Hughes. The circular part of 
// dQdt is determined from a fit to both dLzdt and d(iota)dt using Scott's data. This was found to perform
// better than fitting dQdt directly. The eccentricity terms of dLzdt and dEdt are taken from 2PN results,
// as described in [1]. The eccentricity parts of dQdt are determined by the "spherical" approximation
// dQdt = 2(Q/Lz) dLzdt, with all divergent terms (proportional to 1/cos(iota)) removed (physically, these 
// will be cancelled by eccentricity dependent pieces of d(iota)dt that are also divergent at the pole).
//
// This "new GHK" kludge gives extremely good results for circular orbits, and produce physically reasonable
// evolutions throughout phase space - the pathologies present in the original version of the kludge have 
// been removed.
//

#include <math.h>
#include "Globals.h"
#include "GKR.h"
#include "RRGW.h"
#include "CKG.h"
#include "IEKG.h"
#include "LB.h"

GKR::GKR(const Real semilatus, const Real ecc, const Real energy,
	 const Real angmom, const Real carter, const Real spin) :
  p(semilatus), e(ecc), E(energy), Lz(angmom), Q(carter), kerr_a(spin)
{
  ome = 1. - e;
  ope = 1. + e;
  ome2 = ome*ope;
  a_sm = p/ome2;
  ra = p/ome;
  rp = p/ope;
  Lovercosi=sqrt(Lz*Lz+Q);
  cosiota=Lz/Lovercosi;
  siniotasqr = Q/(Lz*Lz + Q);
  sini=sqrt(siniotasqr);

  CKG *ckg;
  if (kerr_a) {
    ckg=new CKG(0,0,cosiota,USING_cosiota,p,kerr_a);
    Qcirc=ckg->Q;
    Ecirc=ckg->E;
    Lzcirc=ckg->Lz;
    Lovercosicirc=ckg->Lovercosi;
    delete ckg;
  }
  else {
    Ecirc=sqrt(1.-(p-4.)/(p*(p-3.)));
    Lovercosicirc=p/sqrt(p-3.);
    Qcirc=Lovercosicirc*Lovercosicirc*siniotasqr;
    Lzcirc=Lovercosicirc*cosiota;
  }
  Lzdot=Lzdot_GG();
  Qdot=Qdot_GG();
  Edot = Edot_GG();

  if (Q)
    idot = 0.5*sini*cosiota*Qdot/Q-siniotasqr*Lzdot/sqrt(Q);
  else
    idot = 0.;

  cosidot=-1.*sini*idot;
  
  if (e == 0.0) {
    rrgw = new RRGW;
    Real junk;
    //
    // rrgw->Qdotrdot expects the radiated flux of E and Lz, so
    // we need to switch signs.
    //
    rrgw->Qdotrdot(p, kerr_a, Q, E, Lz, -Edot, -Lzdot, junk, pdot);
    edot = 0.0;
    delete rrgw;
  } else {
    //
    // Calculate radot, rpdot:
    //
    radot = -(Numer1(ra)*Edot/Denom(ra) + Numer2(ra)*Lzdot/(Lz*Denom(ra)) + Numer3(ra)*idot/Denom(ra));
    rpdot = -(Numer1(rp)*Edot/Denom(rp) + Numer2(rp)*Lzdot/(Lz*Denom(rp)) + Numer3(rp)*idot/Denom(rp));
    if (rpdot > 0.) {
      cout << "rp growing!" << endl;
    }
    if (radot > 0.) {
      cout << "ra growing!" << endl;
    }
    //
    // Calculate pdot, edot:
    //
    pdot = 0.5*(radot*ome*ome + rpdot*ope*ope);
    edot = 0.5*(radot*ome*ome2 - rpdot*ope*ome2)/p;
  }
}

Real GKR::p_LSO(Real csi) 
{
  const Real r_ipro = Kerr::isco_pro(kerr_a);
  const Real r_iret = Kerr::isco_ret(kerr_a);
  const Real r_horz = Kerr::rplus(kerr_a);
  //
  Real p_l, p_h, p_g,cosiota_g;
  int i;
  IEKG *iekg;
  LB *lb;
  Real *x_l; x_l = Realvector(1, 3);
  Real *x_h; x_h = Realvector(1, 3);
  Real *x_g; x_g = Realvector(1, 3);
  const Real tolx = 1.e-14, tolf = 1.e-14;
  //
  if (e == 0.) {
    //
    // If circular and equatorial, it's pretty trivial.
    //
    /* NB This can cause problems if the offset in the ORBIT_I direction is small and we are centred on an equatorial orbit.
       The computer then takes a positive and negative offset to be the same as the zero offset, so we do not get 
       off-diagonal elements that are 2*gxx, but actually just gxx. Must be wary of this in the future.... */
    
    if (fabs(csi - 1.0) < 1.e-38) {
      p_l = p_h = r_ipro;
      csi = 1.0;
    } else if (fabs(csi + 1.0) < 1.e-38) {
      p_l = p_h = r_iret;
      csi = -1.0;
    } else {
      //
      // If circular and inclined, use the circular LB code.
      //
      // Inner bound: prograde orbit
      //
      p_l = r_ipro;
      x_l[1] = Kerr::Eeqpro(p_l, kerr_a);
      x_l[2] = Kerr::Lzeqpro(p_l, kerr_a);
      x_l[3] = 0.;
      //
      // Outer bound: retrograde orbit
      //
      p_h = r_iret;
      x_h[1] = Kerr::Eeqret(p_h, kerr_a);
      x_h[2] = Kerr::Lzeqret(p_h, kerr_a);
      x_h[3] = 0.;
      //
      while (p_h - p_l > 1.e-12) {
	p_g = 0.5*(p_h + p_l);
	x_g[1] = 0.5*(x_l[1] + x_h[1]);
	x_g[2] = 0.5*(x_l[2] + x_h[2]);
	x_g[3] = 0.5*(x_l[3] + x_h[3]);
	lb=new LB(p_g, kerr_a);
	lb->mnewt(100, x_g, 3, tolx, tolf);
	delete lb;
	cosiota_g = x_g[2]/sqrt(x_g[2]*x_g[2] + x_g[3]);
	if (csi < cosiota_g) { // Move farther out (larger iota)
	  p_l = p_g;
	  x_l[1] = x_g[1];
	  x_l[2] = x_g[2];
	  x_l[3] = x_g[3];
	} else {
	  p_h = p_g;
	  x_h[1] = x_g[1];
	  x_h[2] = x_g[2];
	  x_h[3] = x_g[3];
	}
      }
      p_g = 0.5*(p_h + p_l);
      p_l = p_h = p_g;
    }
  } else {
    // Set p_l to slightly below previous answer, p_h to above.
    p_l = r_ipro;
    p_h = 5.*r_iret;
  }
  // Now, bisect
  while (p_h - p_l > 1.e-14) {
    p_g = 0.5*(p_h + p_l);
    iekg=new IEKG(p_g, e, csi, kerr_a);
    if (iekg->Stable == 1) {
      // Guess is too far out: move p_h in
      p_h = p_g;
    } else if (iekg->Stable == -1) {
      // Guess is too far in: move p_l out
      p_l = p_g;
    } else if (iekg->Stable == 0) {
      // Land right on it ... it could happen!
      p_h = p_l = p_g;
    }
    delete iekg;
  }
  p_g = 0.5*(p_h + p_l);
  return(p_g);
}

Real GKR::Qdot_GG()
{
  Real k,dQdt,circpiece,dLzdEovrcosi,c1a,c1b,c1c,c2a,c2b,c2c,c3a,c3b,c3c,c4a,c4b,c4c,c5a,c5b,c5c,c6a,c6b,c6c,
    c7a,c7b,c7c,c8a,c8b,c8c,c9a,c9b,c9c,c10a,c10b,c10c,c11a,c11b,c11c,higherorderfit;

  Real ic3a,ic3b,ic3c,ic4a,ic4b,ic4c;

  Real csqi=cosiota*cosiota;

  /* Based on a fit to Lzdot and idot from Scott's circular Teukolsky results. */

  ic3a=-0.03093408;
  ic3b=-22.24163077;
  ic3c=7.55265467;
  ic4a=-3.33475762;
  ic4b=22.70130573;
  ic4c=-12.470005617;

  c3a=-28.15174147;
  c3b=60.9607071973;
  c3c=40.99984205;
  c4a=-0.348161211;
  c4b=2.37258476;
  c4c=-66.65840948;
  c5a=-0.715392387;
  c5b=3.21592568;
  c5c=5.28887649;
  c6a=-7.6103411;
  c6b=128.87778309;
  c6c=-475.4650442;
  c7a=12.290783385;
  c7b=-113.1250548;
  c7c=306.11883292;
  c8a=40.9258725;
  c8b=-347.2713496;
  c8c=886.50332051;
  c9a=-25.48313727;
  c9b=224.22721861;
  c9c=-490.98212316;
  c10a=-9.006337706;
  c10b=91.17666278;
  c10c=-297.001939215;
  c11a=-0.64500047;
  c11b=-5.13591989;
  c11c=47.19818628;  

  Real correction,idtcorrection;
  
  correction=kerr_a*((736.2086781-283.9553066*sqrt(p)+kerr_a*(-1325.1852209+483.266206498*sqrt(p)+kerr_a*(634.49936445-219.223848944*sqrt(p))))
    +(82.07804475-25.82025864*sqrt(p)+kerr_a*(-904.16109275+301.477789146*sqrt(p)+kerr_a*(827.31891826-271.9659423*sqrt(p))))*
	      csqi)/(p*sqrt(p));

  idtcorrection=kerr_a*kerr_a*(247.1682656/sqrt(p)-162.2684644+csqi*(-267.5529723/sqrt(p)+184.4645976)
			       +kerr_a*(-182.165263315/sqrt(p)+152.125216225+csqi*(254.0668915/sqrt(p)-188.131613584)))/(p*p);

  higherorderfit=(c3a+(c3b+c3c/sqrt(p))/sqrt(p)+kerr_a*kerr_a*((c4a+(c4b+c4c/sqrt(p))/sqrt(p))+kerr_a*kerr_a*(c5a+(c5b+c5c/sqrt(p))/sqrt(p))))+cosiota*kerr_a*(c6a+(c6b+c6c/sqrt(p))/sqrt(p)+kerr_a*kerr_a*(c7a+(c7b+c7c/sqrt(p))/sqrt(p)))+csqi*kerr_a*kerr_a*((c8a+(c8b+c8c/sqrt(p))/sqrt(p))+kerr_a*kerr_a*(c9a+(c9b+c9c/sqrt(p))/sqrt(p)))+csqi*cosiota*kerr_a*kerr_a*kerr_a*((c10a+(c10b+c10c/sqrt(p))/sqrt(p))+cosiota*kerr_a*(c11a+(c11b+c11c/sqrt(p))/sqrt(p)))+correction;

  circpiece=1. - kerr_a/(pow(p, 1.5))*( cosiota*61./8.) -  
    1247./(336.*p) + 4.*M_PI/pow(p,1.5)-44711./(9072.*p*p)+
    kerr_a*kerr_a*(-57./16.+45.*csqi/8.)/(p*p)+higherorderfit/(p*p*sqrt(p));

  circpiece-=kerr_a*kerr_a*(ic3a+ic3b/p+ic3c/(p*sqrt(p))+kerr_a*cosiota*(ic4a/(sqrt(p))+ic4b/p+ic4c/(p*sqrt(p)))
			    +idtcorrection)/(p*p);

  Qdotcirc=-12.8*Lovercosicirc*siniotasqr*circpiece/pow(p,3.5);

  Real Lzesqcorr=0.;
  dQdt=-12.8*Lovercosi*siniotasqr*(1.-e*e)*sqrt(1.-e*e)*
    (circpiece+7.*e*e/8.-(e*e*425./336.)/p+M_PI*(97.*e*e/8.)/(p*sqrt(p))
     -(302893.*e*e/6048.)/(p*p)-kerr_a*cosiota*e*e*(91./4.+461.*e*e/64.)/(p*sqrt(p))
     +kerr_a*kerr_a*(95.*e*e/16.)/(p*p)+Lzesqcorr)/pow(p,3.5);

  return (dQdt);
}

Real GKR::Edot_GG()
{
  const Real pref= -6.4*pow(ome2, 1.5)/pow(p, 5.);
  const Real term1=e*e*(73./24. + e*e*37./96.);
  const Real term2=e*e*(823./24.+ e*e*(949./32. + e*e*491./192.));
  const Real term3=9181.*e*e/672.;
  const Real term4=1375.*e*e/48.;
  const Real term5=172157*e*e/2592.;
  const Real term6=359.*e*e/32.;
  const Real csqi=cosiota*cosiota;

  Real esqcorr=0.;

  Real fact,dEdLz,dEdQ;

  fact=pow(ome2,1.5);

  dEdLz=(2.*((p*p-2.*p)*Lzcirc+2.*kerr_a*p*Ecirc))/(4.*kerr_a*p*(Lzcirc-kerr_a*Ecirc)-2.*p*p*Ecirc*(p*p+kerr_a*kerr_a));
  dEdQ=(-2.*p+p*p+kerr_a*kerr_a)/(4.*kerr_a*p*(Lzcirc-kerr_a*Ecirc)-2.*p*p*Ecirc*(p*p+kerr_a*kerr_a));

  const Real circbit=-(dEdLz*Ldotcirc+dEdQ*Qdotcirc);

  return(pref*(term1 - kerr_a*cosiota*term2/(pow(p, 1.5))-term3/p+M_PI*term4/pow(p,1.5)
	       -term5/(p*p)+kerr_a*kerr_a*term6/(p*p)+esqcorr)+fact*circbit);
}

Real GKR::Lzdot_GG()
{
  const Real pref = -6.4*pow(ome2, 1.5)/pow(p, 3.5);
  const Real term1 = 0.875*e*e;
  const Real term2 = e*e*(63./8. + e*e*95./64);
  const Real term3 = e*e*(91./4. + e*e*461./64.);
  const Real term4 = 425.*e*e/336.;
  const Real term5 = 97.*e*e/8.;
  const Real term6 = 302893.*e*e/6048.;
  const Real term7 = 95.*e*e/16.;
  const Real term8 = 109.*e*e/16.;
  const Real csqi = cosiota * cosiota;
  
  Real esqcorr=0.;
  Real c1a,c1b,c1c,c2a,c2b,c2c,c3a,c3b,c3c,c4a,c4b,c4c,c5a,c5b,c5c,c6a,c6b,c6c,c7a,c7b,c7c,c8a,c8b,c8c,c9a,c9b,c9c,c10a,c10b,c10c,
    c11a,c11b,c11c,higherorderfit,c1,c2;

  c1a=-10.741956;
  c1b=28.5942157;
  c1c=-9.077378144;
  c1=c1a+(c1b+c1c/sqrt(p))/sqrt(p);

  c2a=-1.428362761;
  c2b=10.70029768;
  c2c=-33.70903016;
  c2=(c2a+(c2b+c2c/sqrt(p))/sqrt(p));

  c3a=-28.15174147;
  c3b=60.9607071973;
  c3c=40.99984205;
  c4a=-0.348161211;
  c4b=2.37258476;
  c4c=-66.65840948;
  c5a=-0.715392387;
  c5b=3.21592568;
  c5c=5.28887649;
  c6a=-7.6103411;
  c6b=128.87778309;
  c6c=-475.4650442;
  c7a=12.290783385;
  c7b=-113.1250548;
  c7c=306.11883292;
  c8a=40.9258725;
  c8b=-347.2713496;
  c8c=886.50332051;
  c9a=-25.48313727;
  c9b=224.22721861;
  c9c=-490.98212316;
  c10a=-9.006337706;
  c10b=91.17666278;
  c10c=-297.001939215;
  c11a=-0.64500047;
  c11b=-5.13591989;
  c11c=47.19818628;  

  Real correction;
  
  correction=kerr_a*((736.2086781-283.9553066*sqrt(p)+kerr_a*(-1325.1852209+483.266206498*sqrt(p)+kerr_a*(634.49936445-219.223848944*sqrt(p))))*cosiota
	      +(82.07804475-25.82025864*sqrt(p)+kerr_a*(-904.16109275+301.477789146*sqrt(p)+kerr_a*(827.31891826-271.9659423*sqrt(p))))*
	      cosiota*csqi)/(p*sqrt(p));

  higherorderfit=kerr_a*(c1+kerr_a*kerr_a*c2)+cosiota*(c3a+(c3b+c3c/sqrt(p))/sqrt(p)+kerr_a*kerr_a*((c4a+(c4b+c4c/sqrt(p))/sqrt(p))+kerr_a*kerr_a*(c5a+(c5b+c5c/sqrt(p))/sqrt(p))))+csqi*kerr_a*(c6a+(c6b+c6c/sqrt(p))/sqrt(p)+kerr_a*kerr_a*(c7a+(c7b+c7c/sqrt(p))/sqrt(p)))+cosiota*csqi*kerr_a*kerr_a*((c8a+(c8b+c8c/sqrt(p))/sqrt(p))+kerr_a*kerr_a*(c9a+(c9b+c9c/sqrt(p))/sqrt(p)))+csqi*csqi*kerr_a*kerr_a*kerr_a*((c10a+(c10b+c10c/sqrt(p))/sqrt(p))+cosiota*kerr_a*(c11a+(c11b+c11c/sqrt(p))/sqrt(p)))+correction;

  const Real circbit=cosiota + kerr_a/(pow(p, 1.5))*(61./24. - csqi*61./8.) -  
    cosiota*1247./(336.*p) + 4.*M_PI*cosiota/pow(p,1.5)-44711.*cosiota/(9072.*p*p)+
    kerr_a*kerr_a*cosiota*(-57./16.+45.*csqi/8.)/(p*p)+higherorderfit/(p*p*sqrt(p));

  Ldotcirc=-6.4*circbit/pow(p,3.5);

  return(pref*(cosiota*term1 + kerr_a/(pow(p, 1.5))*(term2 - csqi*term3) -  
	       cosiota*term4/p +M_PI*term5*cosiota/pow(p,1.5)-term6*cosiota/(p*p)+
	       kerr_a*kerr_a*term7*cosiota/(p*p) +esqcorr + circbit ));
}

//
// The following 4 functions are used for the matrix which
// goes from energy, angmom and Carter constant fluxes to d(ra)/dt, d(rp)/dt.
//

//
// Denominator function.
//
Real GKR::Denom(const Real r)
{
  const Real lzmaE = Lz - kerr_a*E;
  const Real omE2 = (1. - E)*(1. + E);

  const Real tmp1 = Q + lzmaE*lzmaE;
  const Real tmp2 = kerr_a*kerr_a*omE2 + Lz*Lz + Q;

  return(tmp1 + r*(-tmp2 + r*(3. + r*(-2.*omE2))));
}

//
// Numerator 1 function.
//
Real GKR::Numer1(const Real r)
{
  const Real lzmaE = Lz - kerr_a*E;

  return(r*(-2.*kerr_a*lzmaE + r*(kerr_a*kerr_a*E + r*r*E)));
}

//
// Numerator 2 function.
//
Real GKR::Numer2(const Real r)
{
  const Real lzmaE = Lz - kerr_a*E;
  const Real tmp = Lz*Lz + Q;

  return(-kerr_a*kerr_a*Q + r*(2.*(tmp - kerr_a*E*Lz) +
			       r*(-tmp)));
}

//
// Numerator 3 function.
//
Real GKR::Numer3(const Real r)
{
  return(Lz*Lz*sini*(r*(2.-r)-kerr_a*kerr_a)/(cosiota*cosiota*cosiota));
}
