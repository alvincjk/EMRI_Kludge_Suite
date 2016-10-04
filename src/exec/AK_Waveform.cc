#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>

#include "Globals.h"
#include "IEKG.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AK.h"

int main(int argc, char *argv[]){

  // ----- interface -----
  if(argc<2){
    fprintf(stderr,"Input error: Specify path to settings/parameters file\n");
    return 0;
  }
  char *AK_par=argv[1];
  if(CheckFile(AK_par)==0){
    fprintf(stderr,"Input error: Cannot read settings/parameters file\n");
    return 0;
  }
  // ----------

  // ----- load settings and parameters -----
  SetPar AK;
  if(LoadSetPar(&AK,AK_par)==0) return 0;
  // ----------

  // ----- generate AK waveform -----
  clock_t ticks=clock();

  double ang[3];
  IEKG geodesic(AK.p,AK.e,cos(AK.iota),AK.s);
  ParAng(ang,AK.e,AK.iota,AK.gamma,AK.psi,AK.theta_S,AK.phi_S,AK.theta_K,AK.phi_K,AK.alpha,geodesic.zedminus);
  double Phi=ang[0]; // initial mean anomaly
  double nu=pow((1.-AK.e*AK.e)/AK.p,3./2.)/(2.*M_PI*AK.M*SOLARMASSINSEC); // initial Keplerian frequency in Hz

  double *hI,*hII;
  hI=(double*)fftw_malloc(AK.length*sizeof(double));
  hII=(double*)fftw_malloc(AK.length*sizeof(double));
  GenBCWave(hI,hII,AK.dt,AK.length,AK.e,nu,AK.M,AK.mu,AK.s,AK.D,AK.iota,AK.gamma,Phi,AK.theta_S,AK.phi_S,AK.alpha,AK.theta_K,AK.phi_K,AK.LISA,false);

  ticks=clock()-ticks;
  double secs=((double)ticks)/CLOCKS_PER_SEC;
  // ----------

  // ----- output to file -----
  FILE *file;
  char filename[100];
  strcpy(filename,AK.path);
  strcat(filename,"_wave.dat");
  if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
  file=fopen(filename,"w");
  double t=0.;
  for(int i=0;i<AK.length;i++){
    fprintf(file,"%8.6e %14.12e %14.12e\n",t,hI[i],hII[i]);
    t+=AK.dt;
  }
  fclose(file);

  if(AK.traj==true){
    double *nuvec,*evec;
    nuvec=(double*)malloc(AK.length*sizeof(double));
    evec=(double*)malloc(AK.length*sizeof(double));
    GenBCWave(nuvec,evec,AK.dt,AK.length,AK.e,nu,AK.M,AK.mu,AK.s,AK.D,AK.iota,AK.gamma,Phi,AK.theta_S,AK.phi_S,AK.alpha,AK.theta_K,AK.phi_K,AK.LISA,true);
    strcpy(filename,AK.path);
    strcat(filename,"_traj.dat");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    t=0.;
    double dt_RR=0.001; // radiation-reaction timestep for downsampling
    int i_RR=(int)(dt_RR*(SOLARMASSINSEC*AK.M*AK.M/AK.mu)/AK.dt);
    for(int i=0;i<AK.length;i++){
      if(i%i_RR==0 || i+i_RR>=AK.length){
        double p_t=(1.-evec[i]*evec[i])/pow(2.*M_PI*AK.M*SOLARMASSINSEC*nuvec[i],2./3.);
        IEKG geodesic_t(p_t,evec[i],cos(AK.iota),AK.s);
        fprintf(file,"%8.6e %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e\n",t,p_t,evec[i],AK.iota,geodesic_t.E,geodesic_t.Lz,geodesic_t.Q);
      }
      t+=AK.dt;
    }
    fclose(file);
    free(nuvec);
    free(evec);
  }

  if(AK.SNR==true||AK.timing==true){
    strcpy(filename,AK.path);
    strcat(filename,"_info.txt");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    if(AK.SNR==true){
      fftw_plan planhI,planhII;
      planhI=fftw_plan_r2r_1d(AK.length,hI,hI,FFTW_R2HC,FFTW_ESTIMATE);
      planhII=fftw_plan_r2r_1d(AK.length,hII,hII,FFTW_R2HC,FFTW_ESTIMATE);
      fftw_execute(planhI);
      fftw_execute(planhII);
      double hI2=InnProd(hI,hI,AK.dt,AK.length);
      double hII2=InnProd(hII,hII,AK.dt,AK.length);
      fprintf(file,"SNR_I: %.1f\nSNR_II: %.1f\nSNR: %.1f\n\n",sqrt(hI2),sqrt(hII2),sqrt(hI2+hII2));
      fftw_destroy_plan(planhI);
      fftw_destroy_plan(planhII);
    }
    if(AK.timing==true){
      fprintf(file,"Clock ticks: %.2e\nSeconds: %.2fs\n\n",(double)ticks,secs);
    }
    fclose(file);
  }
  // ----------

  fftw_free(hI);
  fftw_free(hII);

}

