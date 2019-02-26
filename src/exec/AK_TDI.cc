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
#include "AKTDI.h"

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

  double *f, *Xf_r,*Xf_im,*Yf_r,*Yf_im,*Zf_r,*Zf_im;
  f=(double*)calloc((AK.length+1)/2,sizeof(double));
  Xf_r=(double*)calloc((AK.length+1)/2,sizeof(double));
  Xf_im=(double*)calloc((AK.length+1)/2,sizeof(double));
  Yf_r=(double*)calloc((AK.length+1)/2,sizeof(double));
  Yf_im=(double*)calloc((AK.length+1)/2,sizeof(double));
  Zf_r=(double*)calloc((AK.length+1)/2,sizeof(double));
  Zf_im=(double*)calloc((AK.length+1)/2,sizeof(double));
  GenAKTDI(f,Xf_r,Xf_im,Yf_r,Yf_im,Zf_r,Zf_im,AK.dt,AK.length,AK.e,nu,AK.M,AK.mu,AK.s,AK.D,AK.iota,AK.gamma,Phi,AK.theta_S,AK.phi_S,AK.alpha,AK.theta_K,AK.phi_K,false,false);

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
  for(int i=0;i<(AK.length+1)/2;i++) fprintf(file,"%14.12e %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e\n",f[i],Xf_r[i],Xf_im[i],Yf_r[i],Yf_im[i],Zf_r[i],Zf_im[i]);
  fclose(file);

  if(AK.timing==true){
    strcpy(filename,AK.path);
    strcat(filename,"_info.txt");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    fprintf(file,"Clock ticks: %.2e\nSeconds: %.6fs\n\n",(double)ticks,secs);
    fclose(file);
  }
  // ----------

  free(f);
  free(Xf_r);
  free(Xf_im);
  free(Yf_r);
  free(Yf_im);
  free(Zf_r);
  free(Zf_im);

}

