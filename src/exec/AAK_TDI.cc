#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAKTDI.h"

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

  // ----- generate AAK waveform -----
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

  double *f, *Xf_r,*Xf_im,*Yf_r,*Yf_im,*Zf_r,*Zf_im;
  f=(double*)calloc((AAK.length+1)/2,sizeof(double));
  Xf_r=(double*)calloc((AAK.length+1)/2,sizeof(double));
  Xf_im=(double*)calloc((AAK.length+1)/2,sizeof(double));
  Yf_r=(double*)calloc((AAK.length+1)/2,sizeof(double));
  Yf_im=(double*)calloc((AAK.length+1)/2,sizeof(double));
  Zf_r=(double*)calloc((AAK.length+1)/2,sizeof(double));
  Zf_im=(double*)calloc((AAK.length+1)/2,sizeof(double));
  GenTDI(f,Xf_r,Xf_im,Yf_r,Yf_im,Zf_r,Zf_im,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.backint,false,false);

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
  for(int i=0;i<(AAK.length+1)/2;i++) fprintf(file,"%14.12e %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e\n",f[i],Xf_r[i],Xf_im[i],Yf_r[i],Yf_im[i],Zf_r[i],Zf_im[i]);
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
  free(f);
  free(Xf_r);
  free(Xf_im);
  free(Yf_r);
  free(Yf_im);
  free(Zf_r);
  free(Zf_im);

}
