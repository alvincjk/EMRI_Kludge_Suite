// EMRI Kludge Suite: Tools

#ifndef _KSTOOLS_H
#define _KSTOOLS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <fstream>
#include <sstream>
#include <string.h>

#define G (6.6726e-11)
#define C (299792458.)
#define Msun (1.9889e30)
#define SOLARMASSINSEC (G*Msun/(C*C*C))
#define F_MIN (1.e-5)
#define F_MAX (3.)

struct SetPar {

  char path[100];

  bool backint;
  bool LISA;
  bool traj;
  bool SNR;
  bool timing;

  int length; // waveform points
  double dt; // time step in seconds
  double p; // initial semi-latus rectum in M
  double T; // waveform duration in years
  double f; // initial GW reference frequency in Hz
  double T_fit; // duration of local fit in radiation-reaction time steps M^2/mu

  double mu; // CO mass in solar masses
  double M; // BH mass in solar masses
  double s; // spin parameter (=a/M=S/M^2)
  double e; // initial eccentricity
  double iota; // inclination angle of L from S
  double gamma; // initial angle of periapsis from LxS
  double psi; // initial true anomaly

  double theta_S; // initial source polar angle in ecliptic coordinates
  double phi_S; // initial source azimuthal angle in ecliptic coordinates
  double theta_K; // initial BH spin polar angle in ecliptic coordinates
  double phi_K; // initial BH spin azimuthal angle in ecliptic coordinates
  double alpha; // initial azimuthal orientation
  double D; // source distance in Gpc

};

int CheckFile(const char *filename);

int LoadSetPar(SetPar *par,const char *filename);

double eLISA_Noise(double f);

double InnProd(double *a,double *b,double dt,int length);

#endif

