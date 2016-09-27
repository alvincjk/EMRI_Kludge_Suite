#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <sstream>
#include <string.h>

#include "KSTools.h"

using namespace std;

int CheckFile(const char *filename){

  FILE *file=fopen(filename,"r");
  if(file!=NULL){
    fclose(file);
    return 1;
  }
  return 0;

}

int LoadSetPar(SetPar *par,const char *filename){

  char str_par[100];
  bool bool_par[4];
  double num_par[19];

  ifstream file(filename);
  string line;
  int str_count=0,bool_count=0,num_count=0;

  while(getline(file,line)){
    if(line.find("#")==string::npos&&!line.empty()){
      istringstream line_stream(line);
      if(str_count<1){
        line_stream>>str_par;
//        fprintf(stderr,"%s\n",str_par);
        if(line_stream.fail()){
          fprintf(stderr,"Parsing error: Check entry %d in settings/parameters file\n",str_count+1);
          return 0;
        }
        str_count++;
      }
      else if(bool_count<4){
        line_stream>>boolalpha>>bool_par[bool_count];
//        fprintf(stderr,"%s\n",bool_par[bool_count]?"true":"false");
        if(line_stream.fail()){
          fprintf(stderr,"Parsing error: Check entry %d in settings/parameters file\n",bool_count+2);
          return 0;
        }
        bool_count++;
      }
      else if(num_count<19){
        line_stream>>num_par[num_count];
//        fprintf(stderr,"%f\n",num_par[num_count]);
        if(line_stream.fail()){
          fprintf(stderr,"Parsing error: Check entry %d in settings/parameters file\n",num_count+6);
          return 0;
        }
        num_count++;
      }
      else{
        fprintf(stderr,"Parsing error: Too many entries in settings/parameters file\n");
        return 0;
      }
    }
  }
  if(str_count+bool_count+num_count<24){
    fprintf(stderr,"Parsing error: Too few entries in settings/parameters file\n");
    return 0;
  }

  strcpy(par->path,str_par);

  par->LISA=bool_par[0];
  par->traj=bool_par[1];
  par->SNR=bool_par[2];
  par->timing=bool_par[3];

  par->length=(int)num_par[0];
  par->dt=num_par[1];
  par->p=num_par[2];
  par->T=num_par[3];
  par->f=num_par[4];
  par->T_fit=num_par[5];

  par->mu=num_par[6];
  par->M=num_par[7];
  par->s=num_par[8];
  par->e=num_par[9];
  par->iota=num_par[10];
  par->gamma=num_par[11];
  par->psi=num_par[12];

  par->theta_S=num_par[13];
  par->phi_S=num_par[14];
  par->theta_K=num_par[15];
  par->phi_K=num_par[16];
  par->alpha=num_par[17];
  par->D=num_par[18];

  if(par->dt<0) par->dt=par->T*31536000./par->length;
  if(par->p<0) par->p=(1.-par->e*par->e)/pow(M_PI*par->M*SOLARMASSINSEC*par->f,2./3.);

  return 1;

}

double eLISA_Noise(double f){

  double L=1.e9;
  double S_acc=1.37e-32*(1.+1.e-4/f)/pow(f,4.);
  double S_sn=5.25e-23;
  double S_omn=6.28e-23;
  double S_n=20./3.*(4.*S_acc+S_sn+S_omn)/L/L*(1.+pow(2.*L*f/0.41/C,2.));
  return S_n;

}

double InnProd(double *a,double *b,double dt,int length){

  double f;
  double innprod=0.;
  for(int i=1;i<(int)((length+1)/2);i++){ 
    f=i/dt/length;
    if(f>F_MIN&&f<F_MAX){
      innprod=innprod+4.*(a[i]*b[i]+a[length-i]*b[length-i])/eLISA_Noise(f);
    }
  }
  innprod=innprod*dt/length;
  return innprod;

}

