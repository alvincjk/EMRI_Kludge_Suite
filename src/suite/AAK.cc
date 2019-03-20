#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#include "IEKG.h"
#include "KSParMap.h"
#include "AAK.h"

using namespace std;


// ----- magnitude of azimuthal angular frequency for prograde/retrograde orbits -----
double OmegaPhi(double v, double e, double cosiota, double s, double M){

  double omegaphi;
  if(cosiota>0) omegaphi=dphidm(v,e,cosiota,s)/dtdm(v,e,cosiota,s)/M;
  else omegaphi=dphidm(v,e,-cosiota,-s)/dtdm(v,e,-cosiota,-s)/M;

  return omegaphi;

}
// ----------


void PNevolution(double *t, double *e, double *v, double *M, double *S, double *gim, double *Phi, double *alp, double *nu, double *gimdotvec, double timestep, int vlength, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, int *i_plunge, int *i_buffer, bool backint){

  double mu=par[0];
  double M0=par[1];
  double S0=par[2];
  double e0=par[3];
  double v0=v_map[0];
  double coslam=cos(par[4]);
  double gim0=par[5];
  double Phi0=par[6];
  double alp0=par[11];

  double edot,vdot,Phidot,gimdot,alpdot;
  double edotp,vdotp,Phidotp,gimdotp,alpdotp;
  double edotm,vdotm,Phidotm,gimdotm,alpdotm;

  double dt_large=dt_map[1];
  double dt_small=dt_map[steps-1]-dt_map[steps-2];
  int j0=(int)((vlength*timestep-dt_map[steps-1])/dt_large+steps+1);
  int j_start=j0;
  int j_min=j_start;
  int j_end=2*j0-1;
  int j_max=j_end;
  int j_max_temp=j_max;
  double t_end=vlength*timestep;

  double *e_AK,*v_AK,*t_fit,*e_fit,*v_fit,*M_fit,*S_fit;
  e_AK=(double*)malloc(2*j0*sizeof(double));
  v_AK=(double*)malloc(2*j0*sizeof(double));
  t_fit=(double*)malloc(2*j0*sizeof(double));
  e_fit=(double*)malloc(2*j0*sizeof(double));
  v_fit=(double*)malloc(2*j0*sizeof(double));
  M_fit=(double*)malloc(2*j0*sizeof(double));
  S_fit=(double*)malloc(2*j0*sizeof(double));

  // ----- evolve AK from t_0 to T_fit -----
  e_AK[j0]=e0;
  v_AK[j0]=v0;
  for(int j=j0;j<j0+steps-1;j++){
    edotm=edot;
    vdotm=vdot;
    edot=dedt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    vdot=dvdt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    if(j==j0){
      edotm=edot;
      vdotm=vdot;
    }
    e_AK[j+1]=e_AK[j]+(1.5*edot-.5*edotm)*(dt_map[j-j0+1]-dt_map[j-j0]);
    v_AK[j+1]=v_AK[j]+(1.5*vdot-.5*vdotm)*(dt_map[j-j0+1]-dt_map[j-j0]);
    if(e_AK[j+1]<0. || v_AK[j+1]<v_AK[j] || isnan(e_AK[j+1]) || isnan(v_AK[j+1])){
      if(j_max==j_max_temp) j_max=j;
      e_AK[j+1]=e_AK[j_max];
      v_AK[j+1]=v_AK[j_max];
    }
  }
  j_max_temp=j_max;
  // ----------

  // ----- compute fit coefficients -----
  int points=min(steps,j_max-j0+1);
  while(e_traj[points-1]>e_traj[points-2]) points--;
  if(points<3){fprintf(stderr,"Fitting error: t_0 too close to plunge\n"); exit(EXIT_FAILURE);}

  double *e_diff,*v_diff,*e_coeff,*v_coeff,*M_coeff,*S_coeff;
  e_diff=(double*)malloc(points*sizeof(double));
  v_diff=(double*)malloc(points*sizeof(double));
  e_coeff=(double*)malloc(2*sizeof(double));
  v_coeff=(double*)malloc(2*sizeof(double));
  M_coeff=(double*)malloc(2*sizeof(double));
  S_coeff=(double*)malloc(2*sizeof(double));

  for(int i=0;i<points;i++){
    e_diff[i]=e_traj[i]-e_AK[i+j0];
    v_diff[i]=v_map[i]-v_AK[i+j0];
  }

  PolyFit(e_coeff,dt_map,e_diff,points);
  PolyFit(v_coeff,dt_map,v_diff,points);
  PolyFit(M_coeff,dt_map,M_map,points);
  PolyFit(S_coeff,dt_map,S_map,points);
  // ----------

  // ----- evolve AK from T_fit to t_end -----
  for(int j=j0+steps-1;j<j_end;j++){
    edotm=edot;
    vdotm=vdot;
    edot=dedt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    vdot=dvdt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    e_AK[j+1]=e_AK[j]+(1.5*edot-.5*edotm)*dt_large;
    v_AK[j+1]=v_AK[j]+(1.5*vdot-.5*vdotm)*dt_large;
    if(e_AK[j+1]<0. || v_AK[j+1]<v_AK[j] || isnan(e_AK[j+1]) || isnan(v_AK[j+1])){
      if(j_max==j_max_temp) j_max=j;
      e_AK[j+1]=e_AK[j_max];
      v_AK[j+1]=v_AK[j_max];
    }
  }
  j_max_temp=j_max;
  // ----------

  // ----- fit AK from t_0 to t_end -----
  for(int j=j0;j<=j_end;j++){
  	double dt;
  	if(j<j0+steps) dt=dt_map[j-j0];
    else dt=dt_map[steps-1]+dt_large*(j-j0-steps+1);
  	double dt2=dt*dt;
    t_fit[j]=dt;
    e_fit[j]=e_AK[j]+e_coeff[0]*dt+e_coeff[1]*dt2;
    v_fit[j]=v_AK[j]+v_coeff[0]*dt+v_coeff[1]*dt2;
    M_fit[j]=M0+(M_coeff[0]*dt+M_coeff[1]*dt2)*SOLARMASSINSEC;
    S_fit[j]=S0+S_coeff[0]*dt+S_coeff[1]*dt2;
    if(e_fit[j]<0. || v_fit[j]<v_fit[max(j0,j-1)]){
      if(j_max==j_max_temp) j_max=j-1;
      e_fit[j]=e_fit[j_max];
      v_fit[j]=v_fit[j_max];
      M_fit[j]=M_fit[j_max];
      S_fit[j]=S_fit[j_max];
    }
  }
  j_max_temp=j_max;
  // ----------

  // ----- check for plunge -----
  int j_RR=(int)((SOLARMASSINSEC*M_phys*SOLARMASSINSEC*M_phys/mu)/dt_large)+1;
  double z1=1.+pow(1.-S_phys*S_phys,1./3.)*(pow(1.+S_phys,1./3.)+pow(1.-S_phys,1./3.));
  double z2=sqrt(3.*S_phys*S_phys+z1*z1);
  double LSO_min=3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2));
  double LSO_max=3.+z2+sqrt((3.-z1)*(3.+z1+2.*z2));

  for(int j=j0;j<=j_max_temp;j++){
    if(j>j_max) break;
    if((1.-e_fit[j]*e_fit[j])*pow(OmegaPhi(v_fit[j],e_fit[j],coslam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.)<LSO_min && j_max==j_max_temp){
      j_min=max(j-j_RR,j0);
      j_max=min(j+1,j_max_temp);
    }
    if(j>j0 && (j-j0)%j_RR==0 && (1.-e_fit[j]*e_fit[j])*pow(OmegaPhi(v_fit[j],e_fit[j],coslam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.)<LSO_max && j_max==j_max_temp){
      j_min=max(j-j_RR,j0);
      IEKG iekg((1.-e_fit[j]*e_fit[j])*pow(OmegaPhi(v_fit[j],e_fit[j],coslam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.),e_fit[j],coslam,S_phys);
      if(iekg.Stable==-1 || iekg.E>1) j_max=min(j+1,j_max_temp);
    }
  }

  while(j_max-j_min>1){
    int j_mid=(j_max+j_min)/2;
    IEKG iekg((1.-e_fit[j_mid]*e_fit[j_mid])*pow(OmegaPhi(v_fit[j_mid],e_fit[j_mid],coslam,S_fit[j_mid],M_fit[j_mid])*M_phys*SOLARMASSINSEC,-2./3.),e_fit[j_mid],coslam,S_phys);
    if(iekg.Stable==-1 || iekg.E>1) j_max=j_mid;
    else j_min=j_mid;
  }
  // ----------

  // ----- evolve and fit AK from t_0 to t_start -----
  if(backint){

  	j_end=j_max;
    t_end=min(t_fit[j_end],vlength*timestep);
    j_start=j0+floor((t_end-vlength*timestep)/dt_large);

    for(int j=j0;j>j_start;j--){
      edotp=edot;
      vdotp=vdot;
      edot=dedt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
      vdot=dvdt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
      if(j==j0){
        edotp=edot;
        vdotp=vdot;
      }
      e_AK[j-1]=e_AK[j]-(1.5*edot-.5*edotp)*dt_large;
      v_AK[j-1]=v_AK[j]-(1.5*vdot-.5*vdotp)*dt_large;
    }

    for(int j=j_start;j<j0;j++){
      double dt_undecayed=dt_large*(j-j0);
      double dt;
      if(j>j0-points) dt=(dt_undecayed+dt_undecayed*(j-j0+points)/points)/2.; // decaying C1 fit
      else dt=-dt_large*points/2.;
      double dt2=dt*dt;
      t_fit[j]=dt_undecayed;
      e_fit[j]=e_AK[j]+e_coeff[0]*dt+e_coeff[1]*dt2;
      v_fit[j]=v_AK[j]+v_coeff[0]*dt+v_coeff[1]*dt2;
      M_fit[j]=M0+(M_coeff[0]*dt+M_coeff[1]*dt2)*SOLARMASSINSEC;
      S_fit[j]=S0+S_coeff[0]*dt+S_coeff[1]*dt2;
    }

  }
  // ----------

  // ----- interpolate AK from t_start to t_end -----
  double *t_in,*e_in,*v_in,*M_in,*S_in;
  t_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  e_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  v_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  M_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  S_in=(double*)malloc((j_end-j_start+1)*sizeof(double));

  for(int j=j_start;j<=j_end;j++){
    t_in[j-j_start]=t_fit[j];
    e_in[j-j_start]=e_fit[j];
    v_in[j-j_start]=v_fit[j];
    M_in[j-j_start]=M_fit[j];
    S_in[j-j_start]=S_fit[j];    
  }

  int i0=vlength-(int)(t_end/timestep);
  t[0]=t_in[0];
  for(int i=1;i<vlength;i++) t[i]=(i-i0)*timestep;
  t[vlength]=t_in[j_end-j_start];

  Interp(t_in,e_in,j_end-j_start+1,t,e,vlength+1);
  Interp(t_in,v_in,j_end-j_start+1,t,v,vlength+1);
  Interp(t_in,M_in,j_end-j_start+1,t,M,vlength+1);
  Interp(t_in,S_in,j_end-j_start+1,t,S,vlength+1);
  // ----------

  // ----- evolve phases from t_0 to t_end -----
  if(j_max==j_max_temp) *i_plunge=vlength-1;
  else *i_plunge=i0+(int)(t_fit[j_min]/timestep);
  *i_buffer=(int)(10./(OmegaPhi(v[*i_plunge],e[*i_plunge],coslam,S[*i_plunge],M[*i_plunge])/2./M_PI)/timestep)+1; // 10 orbits after plunge

  gim[i0]=gim0;
  Phi[i0]=Phi0;
  alp[i0]=alp0;
  for(int i=i0;i<vlength;i++){
    gimdotm=gimdot;
    Phidotm=Phidot;
    alpdotm=alpdot;
    gimdot=(dthetadm(v[i],e[i],coslam,S[i])-drdm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    Phidot=drdm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i];
    alpdot=(dphidm(v[i],e[i],coslam,S[i])-dthetadm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    if(i==i0){
      gimdotm=gimdot;
      Phidotm=Phidot;
      alpdotm=alpdot;
    }
    nu[i]=Phidot/2./M_PI;
    gimdotvec[i]=gimdot;
    gim[i+1]=gim[i]+(1.5*gimdot-.5*gimdotm)*timestep;
    Phi[i+1]=Phi[i]+(1.5*Phidot-.5*Phidotm)*timestep;
    alp[i+1]=alp[i]+(1.5*alpdot-.5*alpdotm)*timestep;
    if(i>*i_plunge+*i_buffer){
      gim[i+1]=gim[*i_plunge];
      Phi[i+1]=Phi[*i_plunge];
      alp[i+1]=alp[*i_plunge];
    }
  }
  nu[vlength]=drdm(v[vlength],e[vlength],coslam,S[vlength])/dtdm(v[vlength],e[vlength],coslam,S[vlength])/(2.*M_PI*M[vlength]);
  gimdotvec[vlength]=(dthetadm(v[vlength],e[vlength],coslam,S[vlength])-drdm(v[vlength],e[vlength],coslam,S[vlength]))/dtdm(v[vlength],e[vlength],coslam,S[vlength])/M[vlength];
  // ----------

  // ----- evolve phases from t_0 to t_start -----
  for(int i=i0;i>0;i--){
    gimdotp=gimdot;
    Phidotp=Phidot;
    alpdotp=alpdot;
    gimdot=(dthetadm(v[i],e[i],coslam,S[i])-drdm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    Phidot=drdm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i];
    alpdot=(dphidm(v[i],e[i],coslam,S[i])-dthetadm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    if(i==i0){
      gimdotp=gimdot;
      Phidotp=Phidot;
      alpdotp=alpdot;
    }
    nu[i]=Phidot/2./M_PI;
    gimdotvec[i]=gimdot;
    gim[i-1]=gim[i]-(1.5*gimdot-.5*gimdotp)*timestep;
    Phi[i-1]=Phi[i]-(1.5*Phidot-.5*Phidotp)*timestep;
    alp[i-1]=alp[i]-(1.5*alpdot-.5*alpdotp)*timestep;
  }
  nu[0]=drdm(v[0],e[0],coslam,S[0])/dtdm(v[0],e[0],coslam,S[0])/(2.*M_PI*M[0]);
  gimdotvec[0]=(dthetadm(v[0],e[0],coslam,S[0])-drdm(v[0],e[0],coslam,S[0]))/dtdm(v[0],e[0],coslam,S[0])/M[0];
  // ----------

  free(e_AK);
  free(v_AK);
  free(t_fit);
  free(e_fit);
  free(v_fit);
  free(M_fit);
  free(S_fit);
  free(e_diff);
  free(v_diff);
  free(e_coeff);
  free(v_coeff);
  free(M_coeff);
  free(S_coeff);
  free(t_in);
  free(e_in);
  free(v_in);
  free(M_in);
  free(S_in);
  return;

}


void waveform(double *t, double *hI, double *hII, double timestep, int vlength, int nmodes, double zeta, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, bool backint, bool mich, bool traj){

  double mu=par[0];
  double lam=par[4];
  double qS=par[7];
  double phiS=par[8];
  double qK=par[9];
  double phiK=par[10];

  double coslam=cos(lam);
  double sinlam=sin(lam);
  double cosqS=cos(qS);
  double sinqS=sin(qS);
  double cosqK=cos(qK);
  double sinqK=sin(qK);
  double cosphiK=cos(phiK);
  double sinphiK=sin(phiK);
  double halfsqrt3=sqrt(3.)/2.;

  double *tvec,*evec,*vvec,*Mvec,*Svec,*gimvec,*Phivec,*alpvec,*nuvec,*gimdotvec;
  tvec=(double *)malloc((vlength+1)*sizeof(double));
  evec=(double *)malloc((vlength+1)*sizeof(double));
  vvec=(double *)malloc((vlength+1)*sizeof(double));
  Mvec=(double *)malloc((vlength+1)*sizeof(double));
  Svec=(double *)malloc((vlength+1)*sizeof(double));
  gimvec=(double *)malloc((vlength+1)*sizeof(double));
  Phivec=(double *)malloc((vlength+1)*sizeof(double));
  alpvec=(double *)malloc((vlength+1)*sizeof(double));
  nuvec=(double *)malloc((vlength+1)*sizeof(double));
  gimdotvec=(double *)malloc((vlength+1)*sizeof(double));

  int i_plunge,i_buffer;
  PNevolution(tvec,evec,vvec,Mvec,Svec,gimvec,Phivec,alpvec,nuvec,gimdotvec,timestep,vlength,par,e_traj,v_map,M_phys,M_map,S_phys,S_map,dt_map,steps,&i_plunge,&i_buffer,backint);

  // ----- output trajectory -----
  if(traj){

    for(int i=0;i<=i_plunge;i++){
      double e=evec[i];
      double v=vvec[i];
      double M=Mvec[i];
      double S=Svec[i];
      t[i]=tvec[i];
      hI[i]=(1.-e*e)*pow(OmegaPhi(v,e,coslam,S,M)*M_phys*SOLARMASSINSEC,-2./3.);
      hII[i]=e;

/*
      // ACCURATE RECOVERY OF P FROM UNPHYSICAL TRAJECTORY (SLOW AND MAY BREAK NEAR PLUNGE, NOT USED TO TEST FOR STABILITY)
      double Omega_map[3],inv_map[3];
      if(coslam>0){
        Omega_map[0]=drdm(v,e,coslam,S)/dtdm(v,e,coslam,S)/2./M_PI;
        Omega_map[1]=dthetadm(v,e,coslam,S)/dtdm(v,e,coslam,S)/2./M_PI;
        Omega_map[2]=dphidm(v,e,coslam,S)/dtdm(v,e,coslam,S)/2./M_PI;
      }
      else{
        Omega_map[0]=drdm(v,e,-coslam,-S)/dtdm(v,e,-coslam,-S)/2./M_PI;
        Omega_map[1]=dthetadm(v,e,-coslam,-S)/dtdm(v,e,-coslam,-S)/2./M_PI;
        Omega_map[2]=-dphidm(v,e,-coslam,-S)/dtdm(v,e,-coslam,-S)/2./M_PI;
      }
      ParInvMap(inv_map,Omega_map,1./v/v,M/SOLARMASSINSEC,S,e,lam);
      hI[i]=1./inv_map[0]/inv_map[0];
      hII[i]=e;
*/

    }

    for(int i=i_plunge+1;i<vlength;i++){
      t[i]=-1.;
      hI[i]=-1.;
      hII[i]=-1.;
    }

    free(tvec);
    free(evec);
    free(vvec);
    free(Mvec);
    free(Svec);
    free(gimvec);
    free(Phivec);
    free(alpvec);
    free(nuvec);
    free(gimdotvec);
    return;

  }
  // ----------

  // ----- compute waveform from t_start to t_end -----
  for(int i=0;i<vlength;i++){

    t[i]=tvec[i];
    hI[i]=0.;
    hII[i]=0.;

    double e=evec[i];
    double v=vvec[i];
    double M=Mvec[i];
    double S=Svec[i];
    double gim=gimvec[i];
    double Phi=Phivec[i];
    double alp=alpvec[i];
    double nu=nuvec[i];
    double gimdot=gimdotvec[i];

    double cosalp=cos(alp);
    double sinalp=sin(alp);
    double cosqL=cosqK*coslam+sinqK*sinlam*cosalp;
    double sinqL=sqrt(1.-cosqL*cosqL);
    double phiLup=sinqK*sinphiK*coslam-cosphiK*sinlam*sinalp-cosqK*sinphiK*sinlam*cosalp;
    double phiLdown=sinqK*cosphiK*coslam+sinphiK*sinlam*sinalp-cosqK*cosphiK*sinlam*cosalp;
    double phiL=atan2(phiLup,phiLdown);
    double Ldotn=cosqL*cosqS+sinqL*sinqS*cos(phiL-phiS);
    double Ldotn2=Ldotn*Ldotn;
    double Sdotn=cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS);
    double betaup=-Sdotn+coslam*Ldotn;
    double betadown=sinqS*sin(phiK-phiS)*sinlam*cosalp+(cosqK*Sdotn-cosqS)/sinqK*sinlam*sinalp;
    double beta=atan2(betaup,betadown);
    double gam=2.*(gim+beta);
    double cos2gam=cos(gam);
    double sin2gam=sin(gam);

    double orbphs,cosorbphs,sinorbphs,FplusI,FcrosI,FplusII,FcrosII;
    if(mich){
      orbphs=2.*M_PI*t[i]/year;
      cosorbphs=cos(orbphs-phiS);
      sinorbphs=sin(orbphs-phiS);
      double cosq=.5*cosqS-halfsqrt3*sinqS*cosorbphs;
      double phiw=orbphs+atan2(halfsqrt3*cosqS+.5*sinqS*cosorbphs,sinqS*sinorbphs);
      double psiup=.5*cosqK-halfsqrt3*sinqK*cos(orbphs-phiK)-cosq*(cosqK*cosqS+sinqK*sinqS*cos(phiK-phiS));
      double psidown=.5*sinqK*sinqS*sin(phiK-phiS)-halfsqrt3*cos(orbphs)*(cosqK*sinqS*sin(phiS)-cosqS*sinqK*sin(phiK))-halfsqrt3*sin(orbphs)*(cosqS*sinqK*cos(phiK)-cosqK*sinqS*cos(phiS));
      double psi=atan2(psiup,psidown);
      double cosq1=.5*(1.+cosq*cosq);
      double cos2phi=cos(2.*phiw);
      double sin2phi=sin(2.*phiw);
      double cos2psi=cos(2.*psi);
      double sin2psi=sin(2.*psi);
      FplusI=cosq1*cos2phi*cos2psi-cosq*sin2phi*sin2psi;
      FcrosI=cosq1*cos2phi*sin2psi+cosq*sin2phi*cos2psi;
      FplusII=cosq1*sin2phi*cos2psi+cosq*cos2phi*sin2psi;
      FcrosII=cosq1*sin2phi*sin2psi-cosq*cos2phi*cos2psi;
    }
    else{
      FplusI=1.;
      FcrosI=0.;
      FplusII=0.;
      FcrosII=1.;
    }
    
    double Amp=pow(OmegaPhi(v,e,coslam,S,M)*M_phys*SOLARMASSINSEC,2./3.)*zeta;

    for(int n=1;n<=nmodes;n++){

      double fn,Doppler,nPhi;
      if(mich){
        fn=n*nu+gimdot/M_PI;
        Doppler=2.*M_PI*fn*AUsec*sinqS*cosorbphs;
        nPhi=n*Phi+Doppler;
      }
      else nPhi=n*Phi;

      double ne=n*e,J[5];
      if(n==1){
        J[0]=-gsl_sf_bessel_J1(ne);
        J[1]=gsl_sf_bessel_J0(ne);
        J[2]=gsl_sf_bessel_J1(ne);
        J[3]=gsl_sf_bessel_Jn(2,ne);
        J[4]=gsl_sf_bessel_Jn(3,ne);
      }
      else gsl_sf_bessel_Jn_array(n-2,n+2,ne,J);
      double a=-n*Amp*(J[0]-2.*e*J[1]+2./n*J[2]+2.*e*J[3]-J[4])*cos(nPhi);
      double b=-n*Amp*sqrt(1-e*e)*(J[0]-2.*J[2]+J[4])*sin(nPhi);
      double c=2.*Amp*J[2]*cos(nPhi);
      double Aplus=-(1.+Ldotn2)*(a*cos2gam-b*sin2gam)+c*(1-Ldotn2);
      double Acros=2.*Ldotn*(b*cos2gam+a*sin2gam);

      // ----- rotate to NK wave frame -----
      double rot[4],Aplusold=Aplus,Acrosold=Acros;
      RotCoeff(rot,lam,qS,phiS,qK,phiK,alp);
      Aplus=Aplusold*rot[0]+Acrosold*rot[1];
      Acros=Aplusold*rot[2]+Acrosold*rot[3];
      // ----------

      double hnI,hnII;
      if(mich){
      	hnI=halfsqrt3*(FplusI*Aplus+FcrosI*Acros);
        hnII=halfsqrt3*(FplusII*Aplus+FcrosII*Acros);
      }
      else{
      	hnI=FplusI*Aplus+FcrosI*Acros;
        hnII=FplusII*Aplus+FcrosII*Acros;
      }

      hI[i]+=hnI;
      hII[i]+=hnII;

    }

    if(i==i_plunge+i_buffer) break;

  }
  // ----------

  // ----- Planck-taper window -----
  double t_plunge=t[i_plunge],t_zero=t_plunge+timestep*i_buffer;
  for(int i=i_plunge+1;i<vlength;i++){
  	t[i]=tvec[i];
    if(i<i_plunge+i_buffer){
      hI[i]=hI[i]/(exp((t_plunge-t_zero)/(t[i]-t_plunge)+(t_plunge-t_zero)/(t[i]-t_zero))+1.);
      hII[i]=hII[i]/(exp((t_plunge-t_zero)/(t[i]-t_plunge)+(t_plunge-t_zero)/(t[i]-t_zero))+1.);
    }
    else{
      hI[i]=0.;
      hII[i]=0.;
    }
  }
  // ----------

  free(tvec);
  free(evec);
  free(vvec);
  free(Mvec);
  free(Svec);
  free(gimvec);
  free(Phivec);
  free(alpvec);
  free(nuvec);
  free(gimdotvec);
  return;

}


void GenWave(double *t, double *hI, double *hII, double timestep, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_phys, double S_map[], double dist, double inc, double gim0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, double dt_map[], int steps, bool backint, bool mich, bool traj){

  double par[12];
  par[0]=mu*SOLARMASSINSEC;
  par[1]=M_map[0]*SOLARMASSINSEC;
  par[2]=S_map[0];
  par[3]=e_traj[0];
  par[4]=inc;
  par[5]=gim0;
  par[6]=Phi0;
  par[7]=qS;
  par[8]=phiS;
  par[9]=qK;
  par[10]=phiK;
  par[11]=alp0;

  // ----- number of modes summed -----
  int nmodes=(int)(30*par[3]);
  if (par[3]<0.135) nmodes=4;
  // ----------

  double zeta=par[0]/dist/Gpc; // M/D

  waveform(t,hI,hII,timestep,vlength,nmodes,zeta,par,e_traj,v_map,M_phys,M_map,S_phys,S_map,dt_map,steps,backint,mich,traj);
  return;

}

