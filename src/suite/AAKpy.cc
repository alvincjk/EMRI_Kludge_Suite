#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAK.h"
#include "AAKPhase.h"
#include "AAKTDI.h"
#include "AKTDI.h"

double AAKwave(SetPar &AAK, double *t, double *hI, double *hII){

  clock_t ticks=clock();

  GKTrajFast gktraj3(cos(AAK.iota),AAK.s);
  gktraj3.p=AAK.p;
  gktraj3.ecc=AAK.e;
  int maxsteps=100;
  int steps=0;
  double dt_fit=min(AAK.T_fit,AAK.length*AAK.dt/SOLARMASSINSEC/AAK.M/AAK.M*AAK.mu)/(maxsteps-1);
  TrajData *traj3;
  traj3=(TrajData*)malloc((size_t)((maxsteps+1)*sizeof(TrajData)));
  gktraj3.Eccentric(dt_fit,traj3,maxsteps,steps);
  double Omega_t[3],ang[3],map_t[3],e_traj[steps],v_map[steps],M_map[steps],s_map[steps],dt_map[steps];
  double Phi;
  for(int i=1;i<=steps;i++){
    IEKG geodesic_t(traj3[i].p,traj3[i].ecc,traj3[i].cosiota,AAK.s);
    geodesic_t.Frequencies(Omega_t);
    if(i==1){
      ParAng(ang,AAK.e,AAK.iota,AAK.gamma,AAK.psi,AAK.theta_S,AAK.phi_S,AAK.theta_K,AAK.phi_K,AAK.alpha,geodesic_t.zedminus);
      Phi=ang[0]; // initial mean anomaly
    }
    ParMap(map_t,Omega_t,traj3[i].p,AAK.M,AAK.s,traj3[i].ecc,AAK.iota);
    e_traj[i-1]=traj3[i].ecc;
    v_map[i-1]=map_t[0]; // mapped initial velocity in c
    M_map[i-1]=map_t[1]; // mapped BH mass in solar masses
    s_map[i-1]=map_t[2]; // mapped spin parameter a/M = S/M^2
    dt_map[i-1]=traj3[i].t*SOLARMASSINSEC*AAK.M*AAK.M/AAK.mu;
  }

  GenWave(t,hI,hII,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.backint,AAK.LISA,false);

  ticks=clock()-ticks;
  return ((double)ticks)/CLOCKS_PER_SEC;

}

double AAKphase(SetPar &AAK, double *t, double *phase_r, double *phase_theta, double *phase_phi, double *omega_r, double *omega_theta, double *omega_phi, double *eccentricity){

  clock_t ticks=clock();

  GKTrajFast gktraj3(cos(AAK.iota),AAK.s);
  gktraj3.p=AAK.p;
  gktraj3.ecc=AAK.e;
  int maxsteps=100;
  int steps=0;
  double dt_fit=min(AAK.T_fit,AAK.length*AAK.dt/SOLARMASSINSEC/AAK.M/AAK.M*AAK.mu)/(maxsteps-1);
  TrajData *traj3;
  traj3=(TrajData*)malloc((size_t)((maxsteps+1)*sizeof(TrajData)));
  gktraj3.Eccentric(dt_fit,traj3,maxsteps,steps);
  double Omega_t[3],ang[3],map_t[3],e_traj[steps],v_map[steps],M_map[steps],s_map[steps],dt_map[steps];
  double Phi;
  for(int i=1;i<=steps;i++){
    IEKG geodesic_t(traj3[i].p,traj3[i].ecc,traj3[i].cosiota,AAK.s);
    geodesic_t.Frequencies(Omega_t);
    if(i==1){
      ParAng(ang,AAK.e,AAK.iota,AAK.gamma,AAK.psi,AAK.theta_S,AAK.phi_S,AAK.theta_K,AAK.phi_K,AAK.alpha,geodesic_t.zedminus);
      Phi=ang[0]; // initial mean anomaly
    }
    ParMap(map_t,Omega_t,traj3[i].p,AAK.M,AAK.s,traj3[i].ecc,AAK.iota);
    e_traj[i-1]=traj3[i].ecc;
    v_map[i-1]=map_t[0]; // mapped initial velocity in c
    M_map[i-1]=map_t[1]; // mapped BH mass in solar masses
    s_map[i-1]=map_t[2]; // mapped spin parameter a/M = S/M^2
    dt_map[i-1]=traj3[i].t*SOLARMASSINSEC*AAK.M*AAK.M/AAK.mu;
  }

  GenPhase(t,phase_r,phase_theta,phase_phi,omega_r,omega_theta,omega_phi,eccentricity,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.backint,false,false);

  ticks=clock()-ticks;
  return ((double)ticks)/CLOCKS_PER_SEC;

}

double AAKTDI(SetPar &AAK, double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im){

  clock_t ticks=clock();

  GKTrajFast gktraj3(cos(AAK.iota),AAK.s);
  gktraj3.p=AAK.p;
  gktraj3.ecc=AAK.e;
  int maxsteps=100;
  int steps=0;
  double dt_fit=min(AAK.T_fit,AAK.length*AAK.dt/SOLARMASSINSEC/AAK.M/AAK.M*AAK.mu)/(maxsteps-1);
  TrajData *traj3;
  traj3=(TrajData*)malloc((size_t)((maxsteps+1)*sizeof(TrajData)));
  gktraj3.Eccentric(dt_fit,traj3,maxsteps,steps);
  double Omega_t[3],ang[3],map_t[3],e_traj[steps],v_map[steps],M_map[steps],s_map[steps],dt_map[steps];
  double Phi;
  for(int i=1;i<=steps;i++){
    IEKG geodesic_t(traj3[i].p,traj3[i].ecc,traj3[i].cosiota,AAK.s);
    geodesic_t.Frequencies(Omega_t);
    if(i==1){
      ParAng(ang,AAK.e,AAK.iota,AAK.gamma,AAK.psi,AAK.theta_S,AAK.phi_S,AAK.theta_K,AAK.phi_K,AAK.alpha,geodesic_t.zedminus);
      Phi=ang[0]; // initial mean anomaly
    }
    ParMap(map_t,Omega_t,traj3[i].p,AAK.M,AAK.s,traj3[i].ecc,AAK.iota);
    e_traj[i-1]=traj3[i].ecc;
    v_map[i-1]=map_t[0]; // mapped initial velocity in c
    M_map[i-1]=map_t[1]; // mapped BH mass in solar masses
    s_map[i-1]=map_t[2]; // mapped spin parameter a/M = S/M^2
    dt_map[i-1]=traj3[i].t*SOLARMASSINSEC*AAK.M*AAK.M/AAK.mu;
  }

  GenTDI(f,Xf_r,Xf_im,Yf_r,Yf_im,Zf_r,Zf_im,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.backint,false,false);

  ticks=clock()-ticks;
  return ((double)ticks)/CLOCKS_PER_SEC;

}

double AKTDI(SetPar &AK, double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im){

  clock_t ticks=clock();

  double ang[3];
  IEKG geodesic(AK.p,AK.e,cos(AK.iota),AK.s);
  ParAng(ang,AK.e,AK.iota,AK.gamma,AK.psi,AK.theta_S,AK.phi_S,AK.theta_K,AK.phi_K,AK.alpha,geodesic.zedminus);
  double Phi=ang[0]; // initial mean anomaly
  double nu=pow((1.-AK.e*AK.e)/AK.p,3./2.)/(2.*M_PI*AK.M*SOLARMASSINSEC); // initial Keplerian frequency in Hz

  GenAKTDI(f,Xf_r,Xf_im,Yf_r,Yf_im,Zf_r,Zf_im,AK.dt,AK.length,AK.e,nu,AK.M,AK.mu,AK.s,AK.D,AK.iota,AK.gamma,Phi,AK.theta_S,AK.phi_S,AK.alpha,AK.theta_K,AK.phi_K,false,false);

  ticks=clock()-ticks;
  return ((double)ticks)/CLOCKS_PER_SEC;

}
