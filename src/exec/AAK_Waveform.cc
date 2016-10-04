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
#include "AAK.h"

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
  int maxsteps=10;
  int steps=0;
  double dt_fit=AAK.T_fit/(SOLARMASSINSEC*AAK.M*AAK.M/AAK.mu)/(maxsteps-1);
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

  double *hI,*hII;
  hI=(double*)fftw_malloc(AAK.length*sizeof(double));
  hII=(double*)fftw_malloc(AAK.length*sizeof(double));
  GenBCWave(hI,hII,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.LISA,false);

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
  double t=0.;
  for(int i=0;i<AAK.length;i++){
    fprintf(file,"%8.6e %14.12e %14.12e\n",t,hI[i],hII[i]);
    t+=AAK.dt;
  }
  fclose(file);

  if(AAK.traj==true){
    double *pvec,*evec;
    pvec=(double*)malloc(AAK.length*sizeof(double));
    evec=(double*)malloc(AAK.length*sizeof(double));
    GenBCWave(pvec,evec,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.LISA,true);
    strcpy(filename,AAK.path);
    strcat(filename,"_traj.dat");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    t=0.;
    double dt_RR=0.001; // radiation-reaction timestep for downsampling
    int i_RR=(int)(dt_RR*(SOLARMASSINSEC*AAK.M*AAK.M/AAK.mu)/AAK.dt);
    int i_max=0;
    while(pvec[i_max]>0.) i_max++;
    for(int i=0;i<i_max;i++){
      if(i%i_RR==0 || i+i_RR>=i_max){
        IEKG geodesic_t(pvec[i],evec[i],cos(AAK.iota),AAK.s);
        fprintf(file,"%8.6e %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e\n",t,pvec[i],evec[i],AAK.iota,geodesic_t.E,geodesic_t.Lz,geodesic_t.Q);
      }
      t+=AAK.dt;
    }
    fclose(file);
    free(pvec);
    free(evec);
  }

  if(AAK.SNR==true||AAK.timing==true){
    strcpy(filename,AAK.path);
    strcat(filename,"_info.txt");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    if(AAK.SNR==true){
      fftw_plan planhI,planhII;
      planhI=fftw_plan_r2r_1d(AAK.length,hI,hI,FFTW_R2HC,FFTW_ESTIMATE);
      planhII=fftw_plan_r2r_1d(AAK.length,hII,hII,FFTW_R2HC,FFTW_ESTIMATE);
      fftw_execute(planhI);
      fftw_execute(planhII);
      double hI2=InnProd(hI,hI,AAK.dt,AAK.length);
      double hII2=InnProd(hII,hII,AAK.dt,AAK.length);
      fprintf(file,"SNR_I: %.1f\nSNR_II: %.1f\nSNR: %.1f\n\n",sqrt(hI2),sqrt(hII2),sqrt(hI2+hII2));
      fftw_destroy_plan(planhI);
      fftw_destroy_plan(planhII);
    }
    if(AAK.timing==true){
      fprintf(file,"Clock ticks: %.2e\nSeconds: %.2fs\n\n",(double)ticks,secs);
    }
    fclose(file);
  }
  // ----------

  free(traj3);
  fftw_free(hI);
  fftw_free(hII);

}

