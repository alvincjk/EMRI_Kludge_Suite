#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>

#include "Globals.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "DopplerShiftedWaveform.h"

int main(int argc, char *argv[]){

  // ----- interface -----
  if(argc<2){
    fprintf(stderr,"Input error: Specify path to settings/parameters file\n");
    return 0;
  }
  char *NK_par=argv[1];
  if(CheckFile(NK_par)==0){
    fprintf(stderr,"Input error: Cannot read settings/parameters file\n");
    return 0;
  }
  // ----------

  // ----- load settings and parameters -----
  SetPar NK;
  if(LoadSetPar(&NK,NK_par)==0) return 0;
  // ----------

  // ----- generate NK waveform -----
  clock_t ticks=clock();

  double ang[3];
  IEKG geodesic(NK.p,NK.e,cos(NK.iota),NK.s);
  ParAng(ang,NK.e,NK.iota,NK.gamma,NK.psi,NK.theta_S,NK.phi_S,NK.theta_K,NK.phi_K,NK.alpha,geodesic.zedminus);
  double chi=ang[1]; // initial polar angle phase in BL coordinates
  double phi=ang[2]; // initial azimuthal angle from nxS in BL coordinates

  Real *hI,*hII;
  hI=(Real*)fftw_malloc(NK.length*sizeof(Real));
  hII=(Real*)fftw_malloc(NK.length*sizeof(Real));
  Real Intrinsparams[]={NK.p,log(NK.mu/NK.M),log(NK.M),NK.e,NK.psi,NK.iota,chi,NK.s};
  Real Extrinsparams[]={cos(NK.theta_S),NK.phi_S,cos(NK.theta_K),NK.phi_K,phi,1000.*NK.D,0.};
  DopplerShiftedWaveform wave(Intrinsparams,Extrinsparams,NK.dt,0,NK.length,0.,70.,0.3,0.7);
  if(NK.LISA==true) wave.hLISA(0,NK.length,hI,hII,0.);
  else wave.hpluscross(0,NK.length,hI,hII);
  ticks=clock()-ticks;
  double secs=((double)ticks)/CLOCKS_PER_SEC;
  // ----------

  // ----- output to file -----
  FILE *file;
  char filename[100];
  strcpy(filename,NK.path);
  strcat(filename,"_wave.dat");
  if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
  file=fopen(filename,"w");
  double t=0.;
  int nancount=0;
  for(int i=0;i<NK.length;i++){
    // ----- check for and fix isolated NaNs -----
    if(isnan(hI[i])){
      nancount++;
      fprintf(stderr,"Waveform warning: NaN found, fixing\n");
      hI[i]=(hI[max(i-1,0)]+hI[min(i+1,NK.length-1)])/2.;
      if(isnan(hI[i])) hI[i]=0.;
    }
    if(isnan(hII[i])){
      nancount++;
      fprintf(stderr,"Waveform warning: NaN found, fixing\n");
      hII[i]=(hII[max(i-1,0)]+hII[min(i+1,NK.length-1)])/2.;
      if(isnan(hII[i])) hII[i]=0.;
    }
    if(nancount>10){
      fprintf(stderr,"Waveform error: Too many NaNs found, stopping\n");
      return 0;
    }
    // ----------
    fprintf(file,"%8.6e %14.12e %14.12e\n",t,hI[i],hII[i]);
    t+=NK.dt;
  }
  fclose(file);

  if(NK.traj==true){
    GKTraj gktraj(cos(NK.iota),NK.s);
    gktraj.p=NK.p;
    gktraj.ecc=NK.e;
    int maxsteps=(int)(NK.length*NK.dt*NK.mu/SOLARMASSINSEC/NK.M/NK.M/0.001)+50000;
    int steps=0;
    TrajData *trajdata;
    trajdata=(TrajData*)malloc((size_t)((maxsteps+1)*sizeof(TrajData)));
    gktraj.Eccentric(0.001,trajdata,maxsteps,steps);
    strcpy(filename,NK.path);
    strcat(filename,"_traj.dat");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    for(int i=1;i<=steps;i++){
      fprintf(file,"%8.6e %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e\n",trajdata[i].t*SOLARMASSINSEC*NK.M*NK.M/NK.mu,trajdata[i].p,trajdata[i].ecc,acos(trajdata[i].cosiota),trajdata[i].E,trajdata[i].Lz,trajdata[i].Q);
    }
    fclose(file);
    free(trajdata);
  }

  if(NK.SNR==true||NK.timing==true){
    strcpy(filename,NK.path);
    strcat(filename,"_info.txt");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    if(NK.SNR==true){
      fftw_plan planhI,planhII;
      planhI=fftw_plan_r2r_1d(NK.length,hI,hI,FFTW_R2HC,FFTW_ESTIMATE);
      planhII=fftw_plan_r2r_1d(NK.length,hII,hII,FFTW_R2HC,FFTW_ESTIMATE);
      fftw_execute(planhI);
      fftw_execute(planhII);
      double hI2=InnProd(hI,hI,NK.dt,NK.length);
      double hII2=InnProd(hII,hII,NK.dt,NK.length);
      fprintf(file,"SNR_I: %.1f\nSNR_II: %.1f\nSNR: %.1f\n\n",sqrt(hI2),sqrt(hII2),sqrt(hI2+hII2));
      fftw_destroy_plan(planhI);
      fftw_destroy_plan(planhII);
    }
    if(NK.timing==true){
      fprintf(file,"Clock ticks: %.2e\nSeconds: %.2fs\n\n",(double)ticks,secs);
    }
    fclose(file);
  }
  // ----------

  free(hI);
  free(hII);

}

