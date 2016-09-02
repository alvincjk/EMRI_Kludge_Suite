#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "Inspiral.h"
#include "Waveform.h"

#define degrees         1.74532925199e-2           /* radians per degree */

Waveform::Waveform(const Real *params, Inspiral *ins) : pE(params),insp(ins)
{
  int i,j;
  Real numfact;
  UseRedShiftedMass=false;

  rtthr=sqrt(3.);
  // No factor of 0.5, since I define my Ax and A+ in a slightly different way to C+L
  numfact=rtthr*pow(G*Msun/(C*C),3)/(1.e6*C*C*pc);

  norm=numfact*exp((*insp).p0[(*insp).LN_ORBIT_ETA])*pow(exp((*insp).p0[(*insp).LN_SYSTEM_M]),3);

  sinthsky=sqrt(1.-pE[COS_THETA_SKY]*pE[COS_THETA_SKY]);
  csthsky=pE[COS_THETA_SKY];
  phik=pE[PHI_SPIN];
  phis=pE[PHI_SKY];
  cphik=cos(phik);
  cphis=cos(phis);
  sphik=sin(phik);
  sphis=sin(phis);
  csthk=pE[COS_THETA_SPIN];
  sinthk=sqrt(1.-csthk*csthk);

  /* Compute inclination of black-hole spin to line of sight as dot product of
     source position and spin direction in ecliptic coordinate system. */

  csinc=sinthsky*sinthk*(cphis*cphik+sphis*sphik)+csthsky*csthk;

  /* Generate default cosmology. */
  H0=70.;
  OmegaM=0.3;
  OmegaL=0.7;

  if (pE[REDSHIFT])
    LCDM=new Cosmology(H0,OmegaM,OmegaL);
}

Waveform::~Waveform()
{

}


void Waveform::ProjectQuadMomTensor(int Nmin, int Nmax, Real **QuadTensor)
{
  int i,j,k;

  Real x0,x1,x2,y0,y1,y2;
  Real cphi0,sphi0,sini;
  Real coeffs[2][5];

  cphi0=cos(pE[ORBIT_PHI0]);
  sphi0=sin(pE[ORBIT_PHI0]);

  if (fabs(csinc) < 1.) 
    sini = sqrt(1.-csinc*csinc); /* precomputation */
  else
    sini=0.;

  /* Compute projections into transverse plane and onto actual x, y axes. */
  /* The plus component is given by 0.5 * d^2 (X^2 - Y^2 )/dt^2 and the cross component by d^2(XY)/dt^2, 
     where X and Y are components relative to an X, Y basis on the sky. This basis is taken to be right handed with 
     respect to the direction of propagation of the gravitational wave (i.e. it is a left handed set with respect 
     to the direction to the source, n ).
     
     We take the X axis to be perpendicular to the spin of the black hole and to our line of sight, with positive
     x in the direction of n^S (this is also how we define phi0). The x[i] coordinates of the inspiral are computed 
     taking phi(0) = 0, so we can construct X[i] = x[i]*cos(phi0)-y[i]*sin(phi0) and 
     Y[i] = -(y[i]*cos(phi0)+x[i]*sin(phi0))*cos(i) + z[i]*sin(i), where i is the inclination angle between the line 
     of sight and the spin axis of the source. 
     
     We may use these definitions and the reduced quadrupole moment tensor in the black hole frame to generate the
     plus and cross components on the sky. We can precompute the necessary coefficients for the rotation. */

  //std::cout << sini << " " << csinc << std::endl;
  //exit(0);
  coeffs[0][0]=0.5*(1.+cphi0*cphi0-csinc*csinc*(1.+sphi0*sphi0));
  coeffs[0][1]=0.5*(1.+sphi0*sphi0-csinc*csinc*(1.+cphi0*cphi0));
  coeffs[0][2]=-cphi0*sphi0*(1.+csinc*csinc);
  coeffs[0][3]=sini*csinc*sphi0;
  coeffs[0][4]=sini*csinc*cphi0;

  coeffs[1][0]=-cphi0*sphi0*csinc;
  coeffs[1][1]=cphi0*sphi0*csinc;
  coeffs[1][2]=(sphi0*sphi0-cphi0*cphi0)*csinc;
  coeffs[1][3]=cphi0*sini;
  coeffs[1][4]=-sphi0*sini;

  for ( i = Nmin; i < Nmax; i++ ){
    for (j=0;j<2;j++) {
      QuadTensor[i-Nmin][j]=0.;
      for (k=0;k<5;k++)
	QuadTensor[i-Nmin][j]+=coeffs[j][k]*((*insp).Quad[k][i]);
    }
  }
  
  return;
}

/* Nmin, Nmax are the start and end points of the trajectory, while t0 is the time at the 0th element of the array, 
   needed only for the modulation by the detector motion. */

int Waveform::hpluscross(int Nmin, int Nmax, Real *hp, Real *hx)
{
  int i,j,k,l,ct,num;
  Real **q;
  q=(Real **)malloc(2*sizeof(Real *));
  for (i=0;i<2;i++)
    q[i]=(Real *)malloc(2*sizeof(Real));

  Real dist,amp;

  /* See comment in hLISA for a note about the factors here. */

  if (pE[REDSHIFT]) {
    if (UseRedShiftedMass)
      dist=LCDM->LumD(pE[REDSHIFT]);
    else
      dist=LCDM->ProperD(pE[REDSHIFT]);
  } 
  else 
    dist=pE[SYSTEM_DP];
  
  /* NB We included a factor of sqrt(3)/2 in the definition of norm which really forms part of the detector
        response, so need to get rid of this now if we only want the h+ and hx waveforms. */

  amp=2.*norm/(dist*sqrt(3.));

  num=Nmax-Nmin;
  /*  cout << Nmin << endl;
  for (ct=0;ct<num;ct++)
  cout << hp[ct+Nmin] << " " << hx[ct+Nmin] << endl;*/

  Real t,csth,phi,psi,fonep,ftwop,fonex,ftwox,orbitalt;
  Real **QuadTensor;
  QuadTensor=(Real **)malloc((num)*sizeof(Real *));
  for (i=0;i<(num);i++)
    QuadTensor[i]=(Real *)malloc(2*sizeof(Real));
  //cout << "Projecting" << endl;
  ProjectQuadMomTensor(Nmin,Nmax,QuadTensor);

  //cout << "Filling waveforms" << endl;
  for (ct=0;ct<num;ct++) {
    hp[ct+Nmin]=amp*QuadTensor[ct][0];
    hx[ct+Nmin]=amp*QuadTensor[ct][1];
  }

  Real xst,yst,zst,xfin,yfin,zfin,p,e,inc,r,deltaphi;
  xst=((*insp).x[Nmin]);
  yst=((*insp).y[Nmin]);
  zst=((*insp).z[Nmin]);
  //cout << (*insp).x[0] << " " << (*insp).y[0] << " " << (*insp).z[0] << endl;
  //cout << xst << " " << yst << " " << zst << endl;

  p=(*insp).p_st;
  e=(*insp).e_st;
  inc=(*insp).i_st;
  
  r=sqrt(xst*xst+yst*yst+zst*zst);
  deltaphi=acos(xst/sqrt(xst*xst+yst*yst));
  if (yst < 0.)
    deltaphi*=-1.;
  Dphi0st=-deltaphi;
  
  if (e==0.)
    chi0_st=0.;
  else {
    chi0_st=acos(zst/(sqrt((*insp).zmin_st)*r));
    if ((*insp).z[Nmin+1]>(*insp).z[Nmin]){
      chi0_st*=-1.;
      //cout << "Changing sign of chi0" << endl;
    }
  }
  if (e==0.)
    psi0_st=0.;
  else {
    psi0_st=acos((p/r -1.)/e);
    if (sqrt((*insp).x[Nmin+1]*(*insp).x[Nmin+1]+(*insp).y[Nmin+1]*(*insp).y[Nmin+1]+(*insp).z[Nmin+1]*(*insp).z[Nmin+1]) <
	sqrt((*insp).x[Nmin]*(*insp).x[Nmin]+(*insp).y[Nmin]*(*insp).y[Nmin]+(*insp).z[Nmin]*(*insp).z[Nmin])) {
      psi0_st*=-1.;
      //cout << "Changing sign of psi0" << endl;
    }
  }

  p_st=p;
  e_st=e;
  i_st=inc;
  
  p=(*insp).p_end;
  e=(*insp).e_end;
  inc=(*insp).i_end;
  
  r=sqrt(xfin*xfin+yfin*yfin+zfin*zfin);
  deltaphi=acos(xfin/sqrt(xfin*xfin+yfin*yfin));
  if (yfin < 0.)
    deltaphi*=-1.;
  Dphi0=deltaphi;
  if (e==0.)
    chi0_end=0.;
  else {
    chi0_end=acos(zfin/(sqrt((*insp).zmin_end)*r));
    if ((*insp).z[num-1]>(*insp).z[num-2]){
      chi0_end*=-1.;
      //cout << "Changing sign of chi0" << endl;
    }
  }
  if (e==0.)
    psi0_end=0.;
  else {
    psi0_end=acos((p/r -1.)/e);
    if (sqrt((*insp).x[num-1]*(*insp).x[num-1]+(*insp).y[num-1]*(*insp).y[num-1]+(*insp).z[num-1]*(*insp).z[num-1]) <
	sqrt((*insp).x[num-2]*(*insp).x[num-2]+(*insp).y[num-2]*(*insp).y[num-2]+(*insp).z[num-2]*(*insp).z[num-2])) {
      psi0_end*=-1.;
      //cout << "Changing sign of psi0" << endl;
    }
  }

  p_end=p;
  e_end=e;
  i_end=inc;

  for (i=0;i<num;i++)
    free(QuadTensor[i]);
  free(QuadTensor);

  return num;
}

void Waveform::ChangeCos(Real Hubb0, Real OmM, Real OmL)
{
  if (pE[REDSHIFT]) {
    delete LCDM;
    LCDM=new Cosmology(Hubb0, OmM, OmL);
    H0=Hubb0;
    OmegaM=OmM;
    OmegaL=OmL;
  }
}
