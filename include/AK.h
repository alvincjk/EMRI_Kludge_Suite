// AK: Waveform

#ifndef _AK_H
#define _AK_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Constants.h"

const int nt = 620000;
const int ndim = 12;
const int modes = 40;

double ArcT(double down, double up);

double J0(double x);

double J1(double x);

double Jn(int n, double x);

void PNevolution(int vlength, double timestep, double *par, double nu0, double *gimdotvec, double *e, double *nu, double *Phi, double *gim, double *alp);

void waveform(double tend,double *par, double nu0, int vlength, double timestep, double *hI, double *hII, int nmodes, double zeta, int WCN, bool mich, bool traj);

void GenBCWave(double *hI, double *hII, double deltat, int vlength, double e0, double nu0, double M, double mu, double S, double dist, double inc, double gam0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, bool mich, bool traj);

#endif
