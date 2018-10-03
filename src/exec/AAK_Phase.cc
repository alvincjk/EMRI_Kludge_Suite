#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAKPhase.h"

int main(int argc, char *argv[]){

  // ----- interface -----
  if(argc<2){
    fprintf(stderr,"Input error: Specify path to settings/parameters file\n");
    return 0;
  }
  char *AAK_par=argv[1];
  if(CheckFile(AAK_par)==0){
    fprintf(stderr,"Input error: Cannot read settings/parameters file\n");
    return 0;
  }
  // ----------

  // ----- load settings and parameters -----
  SetPar AAK;
  if(LoadSetPar(&AAK,AAK_par)==0) return 0;
  // ----------

  // ----- generate AAK phases and frequencies -----
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

  double *t,*phase_r,*phase_theta,*phase_phi,*omega_r,*omega_theta,*omega_phi;
  t=(double*)malloc(AAK.length*sizeof(double));
  phase_r=(double*)malloc(AAK.length*sizeof(double));
  phase_theta=(double*)malloc(AAK.length*sizeof(double));
  phase_phi=(double*)malloc(AAK.length*sizeof(double));
  omega_r=(double*)malloc(AAK.length*sizeof(double));
  omega_theta=(double*)malloc(AAK.length*sizeof(double));
  omega_phi=(double*)malloc(AAK.length*sizeof(double));
  GenPhase(t,phase_r,phase_theta,phase_phi,omega_r,omega_theta,omega_phi,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.backint,AAK.LISA,false);

  ticks=clock()-ticks;
  double secs=((double)ticks)/CLOCKS_PER_SEC;
  // ----------

  // ----- output to file -----
  FILE *file;
  char filename[100];
  strcpy(filename,AAK.path);
  strcat(filename,"_wave.dat");
  if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
  file=fopen(filename,"w");
  for(int i=0;i<AAK.length;i++) fprintf(file,"%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n",t[i],phase_r[i],phase_theta[i],phase_phi[i],omega_r[i],omega_theta[i],omega_phi[i]);
  fclose(file);

  if(AAK.timing==true){
    strcpy(filename,AAK.path);
    strcat(filename,"_info.txt");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    fprintf(file,"Clock ticks: %.2e\nSeconds: %.6fs\n\n",(double)ticks,secs);
    fclose(file);
  }
  // ----------

  free(traj3);
  free(phase_r);
  free(phase_theta);
  free(phase_phi);
  free(omega_r);
  free(omega_theta);
  free(omega_phi);

}

