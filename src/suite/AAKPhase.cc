#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IEKG.h"
#include "KSParMap.h"
#include "AAK.h"
#include "AAKPhase.h"

using namespace std;

void PNevolutionPhase(double *t, double *e, double *v, double *M, double *S, double *gim, double *Phi, double *alp, double *nu, double *gimdotvec, double *Phidotvec, double *alpdotvec, double timestep, int vlength, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, int *i_plunge, int *i_buffer, bool backint){

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
    if(Phidot<=0.) Phidot=Phidotm;
    if(i==i0){
      gimdotm=gimdot;
      Phidotm=Phidot;
      alpdotm=alpdot;
    }
    nu[i]=Phidot/2./M_PI;
    gimdotvec[i]=gimdot;
    Phidotvec[i]=Phidot;
    alpdotvec[i]=alpdot;
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
  Phidotvec[vlength]=drdm(v[vlength],e[vlength],coslam,S[vlength])/dtdm(v[vlength],e[vlength],coslam,S[vlength])/M[vlength];
  alpdotvec[vlength]=(dphidm(v[vlength],e[vlength],coslam,S[vlength])-dthetadm(v[vlength],e[vlength],coslam,S[vlength]))/dtdm(v[vlength],e[vlength],coslam,S[vlength])/M[vlength];
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
    Phidotvec[i]=Phidot;
    alpdotvec[i]=alpdot;
    gim[i-1]=gim[i]-(1.5*gimdot-.5*gimdotp)*timestep;
    Phi[i-1]=Phi[i]-(1.5*Phidot-.5*Phidotp)*timestep;
    alp[i-1]=alp[i]-(1.5*alpdot-.5*alpdotp)*timestep;
  }
  nu[0]=drdm(v[0],e[0],coslam,S[0])/dtdm(v[0],e[0],coslam,S[0])/(2.*M_PI*M[0]);
  gimdotvec[0]=(dthetadm(v[0],e[0],coslam,S[0])-drdm(v[0],e[0],coslam,S[0]))/dtdm(v[0],e[0],coslam,S[0])/M[0];
  Phidotvec[0]=drdm(v[0],e[0],coslam,S[0])/dtdm(v[0],e[0],coslam,S[0])/M[0];
  alpdotvec[0]=(dphidm(v[0],e[0],coslam,S[0])-dthetadm(v[0],e[0],coslam,S[0]))/dtdm(v[0],e[0],coslam,S[0])/M[0];
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


void GenPhase(double *t, double *phase_r, double *phase_theta, double *phase_phi, double *omega_r, double *omega_theta, double *omega_phi, double *eccentricity, double timestep, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_phys, double S_map[], double dist, double inc, double gim0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, double dt_map[], int steps, bool backint, bool mich, bool traj){

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

  double *tvec,*evec,*vvec,*Mvec,*Svec,*gimvec,*Phivec,*alpvec,*nuvec,*gimdotvec,*Phidotvec,*alpdotvec;
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
  Phidotvec=(double *)malloc((vlength+1)*sizeof(double));
  alpdotvec=(double *)malloc((vlength+1)*sizeof(double));

  int i_plunge,i_buffer;
  PNevolutionPhase(tvec,evec,vvec,Mvec,Svec,gimvec,Phivec,alpvec,nuvec,gimdotvec,Phidotvec,alpdotvec,timestep,vlength,par,e_traj,v_map,M_phys,M_map,S_phys,S_map,dt_map,steps,&i_plunge,&i_buffer,backint);

  for(int i=0;i<vlength;i++){
  	t[i]=tvec[i];
  	phase_r[i]=Phivec[i];
  	phase_theta[i]=Phivec[i]+gimvec[i];
  	phase_phi[i]=Phivec[i]+gimvec[i]+alpvec[i];
  	omega_r[i]=Phidotvec[i];
  	omega_theta[i]=Phidotvec[i]+gimdotvec[i];
  	omega_phi[i]=Phidotvec[i]+gimdotvec[i]+alpdotvec[i];
  	eccentricity[i]=evec[i];
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
  free(Phidotvec);
  free(alpdotvec);
  return;

}

