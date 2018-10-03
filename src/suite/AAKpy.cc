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

double AAKphase(SetPar &AAK, double *t, double *phase_r, double *phase_theta, double *phase_phi, double *omega_r, double *omega_theta, double *omega_phi){

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

  GenPhase(t,phase_r,phase_theta,phase_phi,omega_r,omega_theta,omega_phi,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.backint,AAK.LISA,false);

  ticks=clock()-ticks;
  return ((double)ticks)/CLOCKS_PER_SEC;

}