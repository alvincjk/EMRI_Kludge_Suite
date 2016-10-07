#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IEKG.h"
#include "KSParMap.h"
#include "AAK.h"

using namespace std;

// code is modified from AK.cc

/* ************************************************************
   FUNCTION J0
   ************************************************************ */

double J0(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
					    +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
					  +y*(59272.64853+y*(267.8532712+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
				    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			       +y*(-0.6911147651e-5+y*(0.7621095161e-6
						       -y*0.934935152e-7)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}

/* ************************************************************
   FUNCTION J1
   ************************************************************ */

double J1(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
					      +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
					   +y*(99447.43394+y*(376.9991397+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-2.356194491;
    ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			       +y*(0.2457520174e-5+y*(-0.240337019e-6))));
    ans2=0.04687499995+y*(-0.2002690873e-3
			  +y*(0.8449199096e-5+y*(-0.88228987e-6
						 +y*0.105787412e-6)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    if (x < 0.0) ans = -ans;
  }
  return ans;
}

/* ************************************************************
   FUNCTION Jn
   ************************************************************ */

double Jn(int n, double x) {
  const int IACC=160;
  const double BIGNO=1.e10;
  const double BIGNI=1.e-10;
  int j,jsum,m;
  double ax,bj,bjm,bjp,sum,tox,Jnval;
  if(n == 0)
    return J0(x);
  if(n == 1)
    return J1(x);
  ax=fabs(x);
  if(ax == 0.)
    Jnval=0.;
  else if(ax > ((double)n)) {
    tox=2./ax;
    bjm=J0(ax);
    bj=J1(ax);
    for(j=1;j<n;j++) {
      bjp=j*tox*bj-bjm;
      bjm=bj;
      bj=bjp;
    }
    Jnval=bj;
  } else {
    tox=2./ax;
    m=2*((n+((int)(sqrt((double)(IACC*n)))))/2);
    Jnval=0.;
    jsum=0;
    sum=0.;
    bjp=0.;
    bj=1.;
    for (j=m;j>0;j--) {
      bjm=j*tox*bj-bjp;
      bjp=bj;
      bj=bjm;
      if(fabs(bj) > BIGNO) {
	bj=bj*BIGNI;
	bjp=bjp*BIGNI;
	Jnval*=BIGNI;
	sum*=BIGNI;
      }
      if(jsum != 0)
	sum+=bj;
      jsum=1-jsum;
      if(j == n)
	Jnval = bjp;
    }
    sum=2.*sum-bj;
    Jnval/=sum;
  }
  if((x < 0.) && (n & 1)) 
    Jnval*=-1.;
  if(n == -1)
    Jnval*=-1.;
  return Jnval;
}


/* ************************************************************
   SUBROUTINE PNEVOLUTION
   ************************************************************ */

void PNevolution(int vlength, double timestep, double *par, double v_map[], double *gimdotvec, double *e, double *nu, double *Phi,
		 double *gim, double *alp, double *v, double *M, double *S, double e_traj[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, int *i_max) {
  int i,i0;
  double tend;
  double gimdot;
  double edot,vdot,Phidot,alpdot,Mdot,Sdot,edot0,vdot0;
  double edotp,vdotp,Phidotp,gimdotp,alpdotp,Mdotp,Sdotp;
  double edotm,vdotm,Phidotm,gimdotm,alpdotm,Mdotm,Sdotm;
  double t0,mu,M0,e0,Phi0,qS,phiS,gim0;
  double lam,alp0,S0,qK,phiK;
  double Z,Y,hour;
  double cosqS,sinqS,coslam,sinlam,cosqK,sinqK;
  double cosalp0,sinalp0;
  double SdotN,L0dotN,NcrossSdotL0,kappa0;
  double e1,e2,v0,v1,v2,M1,M2,S1,S2,timestep2;

  hour=3600.;
  
  t0=par[1];
  mu=par[2];
  M0=par[3];
  e0=par[4];
  gim0=par[5];
  Phi0=par[6];
  qS=par[7];
  phiS=par[8];
  lam=par[9];
  alp0=par[10];
  S0=par[11];
  qK=par[12];
  phiK=par[13];
  v0=v_map[0];

  cosqS=cos(qS);
  sinqS=sin(qS);
  coslam=cos(lam);
  sinlam=sin(lam);
  cosqK=cos(qK);
  sinqK=sin(qK);
  cosalp0=cos(alp0);
  sinalp0=sin(alp0);

  i0=(int)(t0/timestep);
  timestep2=timestep*timestep;

  // ----- evolve (e_AK,v_AK) -----
  e[i0]=e0;
  v[i0]=v0;
  for(i=i0;i<vlength;i++){
    edotm=edot;
    vdotm=vdot;
    edot=dedt(v[i],e[i],coslam,mu,M0,S0);
    vdot=dvdt(v[i],e[i],coslam,mu,M0,S0);
    if(i==i0){
      edotm=edot;
      vdotm=vdot;
    }
    e[i+1]=e[i]+(1.5*edot-.5*edotm)*timestep;
    v[i+1]=v[i]+(1.5*vdot-.5*vdotm)*timestep;
    if(e[i+1]<0. || isnan(e[i+1]) || isnan(v[i+1])){
      if(*i_max==vlength) *i_max=i;
      e[i+1]=e[*i_max];
      v[i+1]=v[*i_max];
    }
  }
  for(i=i0;i>0;i--){
    edotp=edot;
    vdotp=vdot;
    edot=dedt(v[i],e[i],coslam,mu,M0,S0);
    vdot=dvdt(v[i],e[i],coslam,mu,M0,S0);
    if(i==i0){
      edotp=edot;
      vdotp=vdot;
    }
    e[i-1]=e[i]-(1.5*edot-.5*edotp)*timestep;
    v[i-1]=v[i]-(1.5*vdot-.5*vdotp)*timestep;
  }
  // ----------

  // ----- fit to (e_NK,v_map) -----
  int i_neg=*i_max;
  int points=steps;
  while(timestep*i_neg<dt_map[points-1]) points--;
  int di=(int)(dt_map[points-1]/timestep/(points-1))+1;
  if(points<5){fprintf(stderr,"Fitting error: Reduce duration of local fit\n"); exit(EXIT_FAILURE);}

  double *t_in,*e_in,*v_in,*t_out,*e_out,*v_out,*e_coeff,*v_coeff,*M_coeff,*S_coeff;
  t_in=(double*)malloc(points*sizeof(double));
  e_in=(double*)malloc(points*sizeof(double));
  v_in=(double*)malloc(points*sizeof(double));
  t_out=(double*)malloc(points*sizeof(double));
  e_out=(double*)malloc(points*sizeof(double));
  v_out=(double*)malloc(points*sizeof(double));
  e_coeff=(double*)malloc(4*sizeof(double));
  v_coeff=(double*)malloc(4*sizeof(double));
  M_coeff=(double*)malloc(4*sizeof(double));
  S_coeff=(double*)malloc(4*sizeof(double));

  for(i=0;i<points;i++){
    t_in[i]=timestep*i*di;
    e_in[i]=e[i*di];
    v_in[i]=v[i*di];
    t_out[i]=dt_map[i];
  }
  Interp(t_in,e_in,t_out,e_out,points);
  Interp(t_in,v_in,t_out,v_out,points);

  for(i=0;i<points;i++){
    e_out[i]=e_traj[i]-e_out[i];
    v_out[i]=v_map[i]-v_out[i];
  }
  PolyFit(e_coeff,t_out,e_out,points);
  PolyFit(v_coeff,t_out,v_out,points);

  for(i=0;i<=vlength;i++){
    e[i]=e[i]+e_coeff[0]*timestep*(i-i0)+e_coeff[1]*timestep2*(i-i0)*(i-i0)+e_coeff[2]*timestep*timestep2*(i-i0)*(i-i0)*(i-i0)+e_coeff[3]*timestep2*timestep2*(i-i0)*(i-i0)*(i-i0)*(i-i0);
    v[i]=v[i]+v_coeff[0]*timestep*(i-i0)+v_coeff[1]*timestep2*(i-i0)*(i-i0)+v_coeff[2]*timestep*timestep2*(i-i0)*(i-i0)*(i-i0)+v_coeff[3]*timestep2*timestep2*(i-i0)*(i-i0)*(i-i0)*(i-i0);
    if(i>*i_max){
      e[i]=e[*i_max];
      v[i]=v[*i_max];
    }
    else if(e[i]<0.){
      *i_max=i-1;
      e[i]=e[*i_max];
      v[i]=v[*i_max];
    }
  }
  i_neg=*i_max;
  // ----------

  // ----- evolve (M_map,S_map) -----
  PolyFit(M_coeff,dt_map,M_map,points);
  PolyFit(S_coeff,dt_map,S_map,points);

  for(i=0;i<=vlength;i++){
    M[i]=M0+(M_coeff[0]*timestep*(i-i0)+M_coeff[1]*timestep2*(i-i0)*(i-i0)+M_coeff[2]*timestep*timestep2*(i-i0)*(i-i0)*(i-i0)+M_coeff[3]*timestep2*timestep2*(i-i0)*(i-i0)*(i-i0)*(i-i0))*SOLARMASSINSEC;
    S[i]=S0+S_coeff[0]*timestep*(i-i0)+S_coeff[1]*timestep2*(i-i0)*(i-i0)+S_coeff[2]*timestep*timestep2*(i-i0)*(i-i0)*(i-i0)+S_coeff[3]*timestep2*timestep2*(i-i0)*(i-i0)*(i-i0)*(i-i0);
    if(i>i_neg){
      M[i]=M[i_neg];
      S[i]=S[i_neg];
    }
  }
  // ----------

  // ----- check for plunge -----
  int i_min=0;
  int i_RR=(int)((SOLARMASSINSEC*M_phys*SOLARMASSINSEC*M_phys/mu)/timestep)+1;
  double z1=1.+pow(1.-S_phys*S_phys,1./3.)*(pow(1.+S_phys,1./3.)+pow(1.-S_phys,1./3.));
  double z2=sqrt(3.*S_phys*S_phys+z1*z1);
  double LSO_min=3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2));
  double LSO_max=3.+z2+sqrt((3.-z1)*(3.+z1+2.*z2));

  for(i=i0;i<=vlength;i++){
    if(i>*i_max) break;
    if((1.-e[i]*e[i])*pow(dphidm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i]*M_phys*SOLARMASSINSEC,-2./3.)<LSO_min && *i_max==i_neg){
      i_min=max(i-i_RR,0);
      int i_buffer=(int)(10./(dphidm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/2./M_PI/M[i])/timestep)+1; // 10 extra orbits
      *i_max=min(i+i_buffer,i_neg);
    }
    if(i>i0 && (i-i0)%i_RR==0 && (1.-e[i]*e[i])*pow(dphidm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i]*M_phys*SOLARMASSINSEC,-2./3.)<LSO_max && *i_max==i_neg){
      i_min=max(i-i_RR,0);
      IEKG iekg((1.-e[i]*e[i])*pow(dphidm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i]*M_phys*SOLARMASSINSEC,-2./3.),e[i],coslam,S_phys);
      if(iekg.Stable==-1){
        int i_buffer=(int)(10./(dphidm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/2./M_PI/M[i])/timestep)+1; // 10 extra orbits
        *i_max=min(i+i_buffer,i_neg);
      }
    }
    if(e[i]>e[max(i0,i-1)] || v[i]<v[max(i0,i-1)]){fprintf(stderr,"Fitting error: Adjust duration of local fit\n"); exit(EXIT_FAILURE);}
  }

  while(*i_max-i_min>1){
    int i_mid=(*i_max+i_min)/2;
    IEKG iekg((1.-e[i_mid]*e[i_mid])*pow(dphidm(v[i_mid],e[i_mid],coslam,S[i_mid])/dtdm(v[i_mid],e[i_mid],coslam,S[i_mid])/M[i_mid]*M_phys*SOLARMASSINSEC,-2./3.),e[i_mid],coslam,S_phys);
    if(iekg.Stable==-1) *i_max=i_mid;
    else i_min=i_mid;
  }
  // ----------

  free(t_in);
  free(e_in);
  free(v_in);
  free(t_out);
  free(e_out);
  free(v_out);
  free(e_coeff);
  free(v_coeff);
  free(M_coeff);
  free(S_coeff);

/*
  // OLD QUADRATIC FIT USING DERIVATIVES (NOW DEFUNCT)

  // ----- fit to (e_NK,v_map) -----
  edot0=dedt(v0,e0,coslam,mu,M0,S0);
  vdot0=dvdt(v0,e0,coslam,mu,M0,S0);
  e1=(e_traj[1]-e_traj[0])/dt_map-edot0;
  e2=(e_traj[2]-2.*e_traj[1]+e_traj[0])/2./dt_map/dt_map-dedt2(v0,e0,coslam,mu,M0,S0,vdot0,edot0);
  v1=(v_map[1]-v_map[0])/dt_map-vdot0;
  v2=(v_map[2]-2.*v_map[1]+v_map[0])/2./dt_map/dt_map-dvdt2(v0,e0,coslam,mu,M0,S0,vdot0,edot0);
  for(i=0;i<=vlength;i++){
    e[i]=e[i]+e1*timestep*(i-i0)+e2*timestep2*(i-i0)*(i-i0);
    v[i]=v[i]+v1*timestep*(i-i0)+v2*timestep2*(i-i0)*(i-i0);
  }
  // ----------

  // ----- evolve (M_map,S_map) -----
  M1=(M_map[1]-M_map[0])/dt_map*SOLARMASSINSEC;
  M2=(M_map[2]-2.*M_map[1]+M_map[0])/2./dt_map/dt_map*SOLARMASSINSEC;
  S1=(S_map[1]-S_map[0])/dt_map;
  S2=(S_map[2]-2.*S_map[1]+S_map[0])/2./dt_map/dt_map;
  for(i=0;i<=vlength;i++){
    M[i]=M0+M1*timestep*(i-i0)+M2*timestep2*(i-i0)*(i-i0);
    S[i]=S0+S1*timestep*(i-i0)+S2*timestep2*(i-i0)*(i-i0);
  }
  // ----------
*/

  Phi[i0]=Phi0;
  gim[i0]=gim0;
  alp[i0]=alp0;

  //     <<< EVOLVE FORWARD FROM t0 to tend >>>

  for (i=i0;i<vlength;i++) {
    Phidotm=Phidot;
    gimdotm=gimdot;
    alpdotm=alpdot;

    Phidot=drdm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i];
    alpdot=(dphidm(v[i],e[i],coslam,S[i])-dthetadm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    gimdot=(dthetadm(v[i],e[i],coslam,S[i])-drdm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];

    if (i == i0) {
      Phidotm=Phidot;
      gimdotm=gimdot;
      alpdotm=alpdot;
    }
    nu[i]=Phidot/2./M_PI;
    gimdotvec[i]=gimdot;
    Phi[i+1]=Phi[i]+(1.5*Phidot-.5*Phidotm)*timestep;
    gim[i+1]=gim[i]+(1.5*gimdot-.5*gimdotm)*timestep;
    alp[i+1]=alp[i]+(1.5*alpdot-.5*alpdotm)*timestep;
    if(i>=i_neg){
      Phi[i+1]=Phi[i_neg];
      gim[i+1]=gim[i_neg];
      alp[i+1]=alp[i_neg];
    }
  }
  nu[vlength]=drdm(v[vlength],e[vlength],coslam,S[vlength])/dtdm(v[vlength],e[vlength],coslam,S[vlength])/(2.*M_PI*M[vlength]);
  gimdotvec[vlength]=(dthetadm(v[vlength],e[vlength],coslam,S[vlength])-drdm(v[vlength],e[vlength],coslam,S[vlength]))/dtdm(v[vlength],e[vlength],coslam,S[vlength])/M[vlength];

  //     <<< EVOLVE BACKWARD FROM t0 to t=0 >>>

  for (i=i0;i>0;i--) {
    Phidotp=Phidot;
    gimdotp=gimdot;
    alpdotp=alpdot;

    Phidot=drdm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i];
    alpdot=(dphidm(v[i],e[i],coslam,S[i])-dthetadm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    gimdot=(dthetadm(v[i],e[i],coslam,S[i])-drdm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];

    if (i == i0) {
      Phidotp=Phidot;
      gimdotp=gimdot;
      alpdotp=alpdot;
    }
    nu[i]=Phidot/2./M_PI;
    gimdotvec[i]=gimdot;
    Phi[i-1]=Phi[i]-(1.5*Phidot-.5*Phidotp)*timestep;
    gim[i-1]=gim[i]-(1.5*gimdot-.5*gimdotp)*timestep;
    alp[i-1]=alp[i]-(1.5*alpdot-.5*alpdotp)*timestep;
  }
  nu[0]=drdm(v[0],e[0],coslam,S[0])/dtdm(v[0],e[0],coslam,S[0])/(2.*M_PI*M[0]);
  gimdotvec[0]=(dthetadm(v[0],e[0],coslam,S[0])-drdm(v[0],e[0],coslam,S[0]))/dtdm(v[0],e[0],coslam,S[0])/M[0];

  return;
}


/* ************************************************************
   SUBROUTINE WAVEFORM
   ************************************************************ */

void waveform(double tend,double *par, double v_map[], int vlength, double timestep, double *hI, double *hII, int nmodes, double zeta, double e_traj[],
	      double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, bool mich, bool traj) {
  int i,i0,n;
  double t0,mu,M,qS,phiS,lam,S,qK,phiK,Sn;
  double year,fn,invsqrtS[modes+1];
  double e,v,nu,Phi,gim,alp,t,ne,nPhi;
  double a,b,c,An,Amp,Jm2,Jm1,Jm0,Jp1,Jp2,it1,tint;
  double cosqL,sinqL,Ldotn,phiL,BB,CC,Ldotn2;
  double gam2,Aplus,Acros,Sdotn,PhiT,PhiTdot,alpdot;
  double halfsqrt3,prefact,cosq,cosqS,sinqS,coslam,sinlam;
  double cosqK,sinqK,cosalp,sinalp,cosphiK,sinphiK;
  double orbphs,cosorbphs,sinorbphs,phiw,psi,psidown,psiup;
  double cos2phi,sin2phi,cos2psi,sin2psi,cosq1;
  double FplusI,FcrosI,FplusII,FcrosII,AUsec,Doppler;
  double hnI,hnII,cos2gam,sin2gam,h2I,h3I;
  double betaup,betadown,beta;

  double *evec,*nuvec,*alpvec,*Phivec,*gimvec,*gimdotvec,*vvec,*Mvec,*Svec;
  evec=(double *)malloc((vlength+1)*sizeof(double));
  nuvec=(double *)malloc((vlength+1)*sizeof(double));
  alpvec=(double *)malloc((vlength+1)*sizeof(double));
  Phivec=(double *)malloc((vlength+1)*sizeof(double));
  gimvec=(double *)malloc((vlength+1)*sizeof(double));
  gimdotvec=(double *)malloc((vlength+1)*sizeof(double));
  vvec=(double *)malloc((vlength+1)*sizeof(double));
  Mvec=(double *)malloc((vlength+1)*sizeof(double));
  Svec=(double *)malloc((vlength+1)*sizeof(double));

  halfsqrt3=0.5*sqrt(3.);

  if (mich)
    prefact=halfsqrt3;
  else
    prefact=1.;

  year=31536000.;
  AUsec=499.004783702731;
  
  t0=par[1];
  mu=par[2];
  qS=par[7];
  phiS=par[8];
  lam=par[9];
  qK=par[12];
  phiK=par[13];

  coslam=cos(lam);
  sinlam=sin(lam);
  cosqS=cos(qS);
  sinqS=sin(qS);
  cosqK=cos(qK);
  sinqK=sin(qK);
  cosphiK=cos(phiK);
  sinphiK=sin(phiK);

  int *i_max=(int*)malloc(sizeof(int));
  *i_max=vlength;
  PNevolution(vlength,timestep,par,v_map,gimdotvec,evec,nuvec,Phivec,gimvec,alpvec,vvec,Mvec,Svec,e_traj,M_phys,M_map,S_phys,S_map,dt_map,steps,i_max);
  int i_buffer=(int)(10./(dphidm(vvec[*i_max],evec[*i_max],coslam,Svec[*i_max])/dtdm(vvec[*i_max],evec[*i_max],coslam,Svec[*i_max])/2./M_PI/Mvec[*i_max])/timestep)+1; // 10 orbits after LSO

  // ----- output trajectory -----
  if(traj==true){
    for(i=0;i<*i_max;i++){
      e=evec[i];
      v=vvec[i];
      M=Mvec[i];
      S=Svec[i];
      hI[i]=(1.-e*e)*pow(dphidm(v,e,coslam,S)/dtdm(v,e,coslam,S)/M*M_phys*SOLARMASSINSEC,-2./3.);
      hII[i]=e;
/*
      // ACCURATE RECOVERY OF P FROM UNPHYSICAL TRAJECTORY (SLOW AND MAY BREAK NEAR PLUNGE, NOT USED TO TEST FOR STABILITY)

      double Omega_map[3],inv_map[3];
      Omega_map[0]=drdm(v,e,coslam,S)/dtdm(v,e,coslam,S)/2./M_PI;
      Omega_map[1]=dthetadm(v,e,coslam,S)/dtdm(v,e,coslam,S)/2./M_PI;
      Omega_map[2]=dphidm(v,e,coslam,S)/dtdm(v,e,coslam,S)/2./M_PI;
      ParInvMap(inv_map,Omega_map,1./v/v,M/SOLARMASSINSEC,S,e,lam);
      hI[i]=1./inv_map[0]/inv_map[0];
      hII[i]=e;
*/
    }
    for(i=*i_max;i<=vlength;i++){
      hI[i]=-1.;
      hII[i]=-1.;
    }
    free(evec);
    free(nuvec);
    free(alpvec);
    free(Phivec);
    free(gimvec);
    free(gimdotvec);
    free(vvec);
    free(Mvec);
    free(Svec);
    free(i_max);
    return;
  }
  // ----------

  i0=(int)(t0/timestep);
  h2I=h3I=0.;

  for(i=i0;i<=vlength;i++) {

    hI[i]=0.;
    hII[i]=0.;

    t=timestep*(double)i;

    e=evec[i];
    v=vvec[i];
    nu=nuvec[i];
    Phi=Phivec[i];
    gim=gimvec[i];
    alp=alpvec[i];
    M=Mvec[i];
    S=Svec[i];
    cosalp=cos(alp);
    sinalp=sin(alp);

    cosqL=cosqK*coslam+sinqK*sinlam*cosalp;
    sinqL=sqrt(1-cosqL*cosqL);
    BB=sinqK*cosphiK*coslam+sinphiK*sinlam*sinalp-cosqK*cosphiK*sinlam*cosalp;
    CC=sinqK*sinphiK*coslam-cosphiK*sinlam*sinalp-cosqK*sinphiK*sinlam*cosalp;
    phiL=atan2(CC,BB);
    Ldotn=cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS);
    Ldotn2=Ldotn*Ldotn;

    if (mich) {
      orbphs=2.*M_PI*t/year;
      cosorbphs=cos(orbphs-phiS);
      sinorbphs=sin(orbphs-phiS);
      cosq=.5*cosqS-halfsqrt3*sinqS*cosorbphs;
      phiw=orbphs+atan2(halfsqrt3*cosqS+.5*sinqS*cosorbphs,sinqS*sinorbphs);
      psiup=.5*cosqL-halfsqrt3*sinqL*cos(orbphs-phiL)-cosq*(cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS));
      psidown=.5*sinqL*sinqS*sin(phiL-phiS)-halfsqrt3*cos(orbphs)*(cosqL*sinqS*sin(phiS)-cosqS*sinqL*sin(phiL))-halfsqrt3*sin(orbphs)*(cosqS*sinqL*cos(phiL)-cosqL*sinqS*cos(phiS));
      psi=atan2(psiup,psidown);
      cosq1=.5*(1+cosq*cosq);
      cos2phi=cos(2.*phiw);
      sin2phi=sin(2.*phiw);
      cos2psi=cos(2.*psi);
      sin2psi=sin(2.*psi);

      FplusI=cosq1*cos2phi*cos2psi-cosq*sin2phi*sin2psi;
      FcrosI=cosq1*cos2phi*sin2psi+cosq*sin2phi*cos2psi;
      FplusII=cosq1*sin2phi*cos2psi+cosq*cos2phi*sin2psi;
      FcrosII=cosq1*sin2phi*sin2psi-cosq*cos2phi*cos2psi;
    } else {
      FplusI=1.;
      FcrosI=0.;
      FplusII=0.;
      FcrosII=1.;
    }

    Sdotn=cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS);
    betaup=-Sdotn+coslam*Ldotn;
    betadown=sinqS*sin(phiK-phiS)*sinlam*cosalp+(cosqK*Sdotn-cosqS)/sinqK*sinlam*sinalp;
    beta=atan2(betaup,betadown);
    gam2=2.*(gim+beta);
    cos2gam=cos(gam2);
    sin2gam=sin(gam2);
    
    Amp=pow(dphidm(v,e,coslam,S)/dtdm(v,e,coslam,S)/M*M_phys*SOLARMASSINSEC,2./3.)*zeta; // amplitude uses azimuthal frequency and physical mass rather than radial frequency and mapped mass

    for(n=1;n<nmodes+1;n++) {
      fn=n*nu+gimdotvec[i]/M_PI;
      Doppler=2.*M_PI*fn*AUsec*sinqS*cosorbphs;
      if (mich)
	      nPhi=n*Phi+Doppler;
      else
	      nPhi=n*Phi;
      ne=n*e;
      Jm2=Jn(n-2,ne);
      Jm1=Jn(n-1,ne);
      Jm0=Jn(n,ne);
      Jp1=Jn(n+1,ne);
      Jp2=Jn(n+2,ne);   
      a=-n*Amp*(Jm2-2.*e*Jm1+2./n*Jm0+2.*e*Jp1-Jp2)*cos(nPhi);
      b=-n*Amp*sqrt(1-e*e)*(Jm2-2.*Jm0+Jp2)*sin(nPhi);
      c=2.*Amp*Jm0*cos(nPhi);

      Aplus=-(1.+Ldotn2)*(a*cos2gam-b*sin2gam)+c*(1-Ldotn2);
      Acros=2.*Ldotn*(b*cos2gam+a*sin2gam);

      hnI=prefact*(FplusI *Aplus+FcrosI *Acros);
      hnII=prefact*(FplusII*Aplus+FcrosII*Acros);
      if (n==2)
	h2I+=(hnI*hnI);
      if (n==3)
	h3I+=(hnI*hnI);
      hI[i]=hI[i]+hnI;
      hII[i]=hII[i]+hnII;
    }

    if(i==*i_max+i_buffer) break;

  }

  // ----- Planck-taper window -----
  double t_LSO=timestep*(*i_max);
  double t_end=timestep*(*i_max+i_buffer);
  for(i=*i_max+1;i<=vlength;i++){
    if(i<*i_max+i_buffer){
      hI[i]=hI[i]/(exp((t_LSO-t_end)/(timestep*i-t_LSO)+(t_LSO-t_end)/(timestep*i-t_end))+1.);
      hII[i]=hII[i]/(exp((t_LSO-t_end)/(timestep*i-t_LSO)+(t_LSO-t_end)/(timestep*i-t_end))+1.);
    }
    else{
      hI[i]=0.;
      hII[i]=0.;
    }
  }
  // ----------

  for(i=i0-1;i>=0;i--) {

    hI[i]=0.;
    hII[i]=0.;

    t=timestep*(double)i;

    e=evec[i];
    v=vvec[i];
    nu=nuvec[i];
    Phi=Phivec[i];
    gim=gimvec[i];
    alp=alpvec[i];
    M=Mvec[i];
    S=Svec[i];
    cosalp=cos(alp);
    sinalp=sin(alp);

    cosqL=cosqK*coslam+sinqK*sinlam*cosalp;
    sinqL=sqrt(1-cosqL*cosqL);
    BB=sinqK*cosphiK*coslam+sinphiK*sinlam*sinalp-cosqK*cosphiK*sinlam*cosalp;
    CC=sinqK*sinphiK*coslam-cosphiK*sinlam*sinalp-cosqK*sinphiK*sinlam*cosalp;
    phiL=atan2(CC,BB);
    Ldotn=cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS);
    Ldotn2=Ldotn*Ldotn;

    if (mich) {
      orbphs=2.*M_PI*t/year;
      cosorbphs=cos(orbphs-phiS);
      sinorbphs=sin(orbphs-phiS);
      cosq=.5*cosqS-halfsqrt3*sinqS*cosorbphs;
      phiw=orbphs+atan2(halfsqrt3*cosqS+.5*sinqS*cosorbphs,sinqS*sinorbphs);
      psiup=.5*cosqL-halfsqrt3*sinqL*cos(orbphs-phiL)-cosq*(cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS));
      psidown=.5*sinqL*sinqS*sin(phiL-phiS)-halfsqrt3*cos(orbphs)*(cosqL*sinqS*sin(phiS)-cosqS*sinqL*sin(phiL))
        -halfsqrt3*sin(orbphs)*(cosqS*sinqL*cos(phiL)-cosqL*sinqS*cos(phiS));
      psi=atan2(psiup,psidown);
      cosq1=.5*(1+cosq*cosq);
      cos2phi=cos(2.*phiw);
      sin2phi=sin(2.*phiw);
      cos2psi=cos(2.*psi);
      sin2psi=sin(2.*psi);

      FplusI=cosq1*cos2phi*cos2psi-cosq*sin2phi*sin2psi;
      FcrosI=cosq1*cos2phi*sin2psi+cosq*sin2phi*cos2psi;
      FplusII=cosq1*sin2phi*cos2psi+cosq*cos2phi*sin2psi;
      FcrosII=cosq1*sin2phi*sin2psi-cosq*cos2phi*cos2psi;
    } else {
      FplusI=1.;
      FcrosI=0.;
      FplusII=0.;
      FcrosII=1.;
    }
      
    Sdotn=cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS);
    betaup=-Sdotn+coslam*Ldotn;
    betadown=sinqS*sin(phiK-phiS)*sinlam*cosalp+(cosqK*Sdotn-cosqS)/sinqK*sinlam*sinalp;
    beta=atan2(betaup,betadown);
    gam2=2.*(gim+beta);
    cos2gam=cos(gam2);
    sin2gam=sin(gam2);

    Amp=pow(dphidm(v,e,coslam,S)/dtdm(v,e,coslam,S)/M*M_phys*SOLARMASSINSEC,2./3.)*zeta; // amplitude uses azimuthal frequency and physical mass rather than radial frequency and mapped mass

    for(n=1;n<nmodes+1;n++) {
      fn=n*nu+gimdotvec[i]/M_PI;
      Doppler=2.*M_PI*fn*AUsec*sinqS*cosorbphs;
      if (mich)
	      nPhi=n*Phi+Doppler;
      else
	      nPhi=n*Phi;
      ne=n*e;
      Jm2=Jn(n-2,ne);
      Jm1=Jn(n-1,ne);
      Jm0=Jn(n,ne);
      Jp1=Jn(n+1,ne);
      Jp2=Jn(n+2,ne);
      a=-n*Amp*(Jm2-2.*e*Jm1+2./n*Jm0+2.*e*Jp1-Jp2)*cos(nPhi);
      b=-n*Amp*sqrt(1-e*e)*(Jm2-2.*Jm0+Jp2)*sin(nPhi);
      c=2.*Amp*Jm0*cos(nPhi);

      Aplus=-(1.+Ldotn2)*(a*cos2gam-b*sin2gam)+c*(1-Ldotn2);
      Acros=2.*Ldotn*(b*cos2gam+a*sin2gam);

      hnI=prefact*(FplusI *Aplus+FcrosI *Acros);
      hnII=prefact*(FplusII*Aplus+FcrosII*Acros);

      hI[i]=hI[i]+hnI;
      hII[i]=hII[i]+hnII;
    }
  }

  free(evec);
  free(nuvec);
  free(alpvec);
  free(Phivec);
  free(gimvec);
  free(gimdotvec);
  free(vvec);
  free(Mvec);
  free(Svec);
  free(i_max);
  return;

}


void GenBCWave(double *hI, double *hII, double deltat, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_phys, double S_map[], double dist, double inc, double gam0, double Phi0,
	       double qS, double phiS, double alp0, double qK, double phiK, double dt_map[], int steps, bool mich, bool traj) {
  int i,n,nmodes,p,pp,m,k,iv;
  int dim,direction,comp,pointM,pointe0;
  int timestep1,nsteps1;
  double par[14],dt,t,spinres,dirres;

  double bigest,count,sqrtdet;
  double integral,StoN,eisco,nuisco,gimdot0;
  double delta,deltahalf,tend,q[14],element;
  double parplus[14],parminus[14],det;
  double normI,normII,norm,dtdays;
  double t0res,mures,Mres,e0res,gim0res,Phi0res;
  double lamres,alp0res,OmegaSres,Dres,zeta,invD,Gpsc;

  /* -------------------------------------------------------------
     SETUP OF PHYSICAL AND NUMERICAL PARAMETERS
     -------------------------------------------------------------*/

  dt=deltat*((double)vlength);
  tend=dt;

  //     t0
  par[1]=250.; // constant time offset to ensure far-field agreement with NK at t=0 - need to track this down
  //     mu
  par[2]=mu*SOLARMASSINSEC;
  //     M
  par[3]=M_map[0]*SOLARMASSINSEC;
  //     e0
  par[4]=e_traj[0];
  //     S/M^2
  par[11]=S_map[0];
  //     inc
  par[9]=inc;

  invD=1./(dist*Gpc);
  zeta=par[2]*invD;

  //     gam0
  par[5]=gam0;
  //     Phi0
  par[6]=Phi0;
  //     qS
  par[7]=qS;
  //     phiS
  par[8]=phiS;
  //     alp0
  par[10]=alp0;
  //     qK
  par[12]=qK;
  //     phiK
  par[13]=phiK;

  /* -------------------------------------------------------------
     NUMBER OF MODES TO BE CONSIDERED (as a function of e0)
     ------------------------------------------------------------- */

  nmodes=(int)(30*par[4]);
  if (par[4] < 0.135) nmodes=4;

  /* --------------------------------------------------------------
     CALCULATING WAVEFORM & DERIVATIVES
     -------------------------------------------------------------- */

  waveform(tend,par,v_map,vlength,deltat,hI,hII,nmodes,zeta,e_traj,M_phys,M_map,S_phys,S_map,dt_map,steps,mich,traj);
  return;
}

