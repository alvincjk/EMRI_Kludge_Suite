#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Globals.h"
#include "NRUtil.h"
#include "WDBinary.h"

#define SOLARMASSINSEC (G*Msun/(C*C*C)) /* 1 solar mass in units of seconds. */
#define SOLARMASSINM (G*Msun/(C*C)) /* 1 solar mass in units of metres. */

WDBinary::WDBinary(const Real freq, const Real totmass, const Real massrat, const Real csthetaL, const Real phiL, const Real csthetaS, const Real phiS, 
		   const Real dist, const Real phase):f(freq),M(totmass),eta(massrat),csthl(csthetaL),pL(phiL),csthsky(csthetaS),pS(phiS),d(dist),phi0(phase)
{
  /* Compute various useful quantities. */

  /* We assume supplied frequency is the GW frequency. */
  forb=f/2.;
  ainM=pow(M_PI*f*M*SOLARMASSINSEC,-2./3.);
  a=ainM*M*SOLARMASSINM;

  /* Assume d was supplied in kpc and convert to physical units, i.e., m */
  dphys=d*1000.*pc;

  cphis=cos(pS);
  sphis=sin(pS);
  sinthsky=sqrt(1.-csthsky*csthsky);
  cphil=cos(pL);
  sphil=sin(pL);
  sinthl=sqrt(1.-csthl*csthl);

  cosi=csthl*csthsky+sinthl*sinthsky*(cphil*cphis+sphil*sphis);
  amp=2.*eta*M*SOLARMASSINM/((1.+eta)*(1.+eta)*ainM*dphys);
}

WDBinary::~WDBinary()
{

}

/* Currently returns plus and cross polarisations with respect to the principal axes of the binary. */

void WDBinary::hpluscross(Real *hp, Real *hx, Real dt, int npts)
{
  Real Ap,Ax,t,pls,crs;
  int i;

  Ap=amp*(1.+cosi*cosi);
  Ax=-2.*amp*cosi;

  t=0.;
  for (i=0;i<npts;i++) {
    pls=Ap*cos(2.*M_PI*f*t+phi0);
    crs=Ax*sin(2.*M_PI*f*t+phi0);
    hp[i]=pls;
    hx[i]=crs;
    //hp[i]=c2psi*pls-s2psi*crs;
    //hx[i]=s2psi*pls+c2psi*crs;
    t+=dt;
  }

  return;
}

void WDBinary::hLISA(Real *hI, Real *hII, Real dt, int npts, Real t0)
{
  Real Ap,Ax,t;
  Real orbitalt,cost,sint,csth,tdphi,tphi,c2phi,s2phi,tpsi,c2psi,s2psi,fonex,fonep,ftwox,ftwop,phiD,doppf,rtthr;
  int i;

  t=t0;
  doppf=2.*M_PI*f*AU*sinthsky/C;
  rtthr=sqrt(3.);

  Ap=0.5*rtthr*amp*(1.+cosi*cosi);
  Ax=-rtthr*amp*cosi;

  for (i=0;i<npts;i++) {
    t+=dt;
    orbitalt=t/YEAR;
    cost=cos(2.*M_PI*orbitalt);
    sint=sin(2.*M_PI*orbitalt);
    csth=0.5*csthsky-0.5*rtthr*sinthsky*(cost*cphis+sint*sphis);
    tdphi=-0.5*(cost*cphis+sint*sphis+rtthr*csthsky/sinthsky)/(sint*cphis-cost*sphis);
    tphi=(sint-cost*tdphi)/(cost+sint*tdphi);
    c2phi=2./(1.+tphi*tphi)-1.;
    s2phi=2.*tphi/(1.+tphi*tphi);
    tpsi= ( 0.5*csthl-0.5*rtthr*sinthl*(cost*cphil+sint*sphil)
	      -csth*(csthl*csthsky+sinthl*sinthsky*(cphil+cphis+sphil*sphis) ) )/
		(0.5*sinthl*sinthsky*(sphil*cphis-sphis*cphil)-0.5*rtthr*
		 cost*(csthl*sinthsky*sphis-csthsky*sinthl*sphil)
		 -0.5*rtthr*sint*(csthsky*sinthl*cphil-csthl*sinthsky*cphis));
    c2psi=2./(1.+tpsi*tpsi)-1.;
    s2psi=2.*tpsi/(1.+tpsi*tpsi);
    phiD=doppf*(cost*cphis+sint*sphis);
    if (hI) {
      fonep=0.5*(1.+csth*csth)*c2phi*c2psi-csth*s2phi*s2psi;
      fonex=0.5*(1.+csth*csth)*c2phi*s2psi +csth*s2phi*c2psi;
      hI[i]=Ap*fonep*cos(2.*M_PI*f*t+phi0+phiD)+Ax*fonex*sin(2.*M_PI*f*t+phi0+phiD);
    }
    if (hII) {
      ftwop=0.5*(1.+csth*csth)*s2phi*c2psi+csth*c2phi*s2psi;
      ftwox=0.5*(1.+csth*csth)*s2phi*s2psi-csth*c2phi*c2psi;
      hII[i]=Ap*ftwop*cos(2.*M_PI*f*t+phi0+phiD)+Ax*ftwox*sin(2.*M_PI*f*t+phi0+phiD);
    }
  }

  return;
}
