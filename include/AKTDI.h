// AK: TDIs

#ifndef _AKTDI_H
#define _AKTDI_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>

#include "Constants.h"

using namespace std;

void PNevolutionAKTDI(complex<double> **Xp, complex<double> **Xc, complex<double> **Yp, complex<double> **Yc, complex<double> **Zp, complex<double> **Zc,\
 double *Ampvec, double *tvec, double *evec, double *gimvec, double *Phivec, double *alpvec, double *nuvec, double *gimdotvec, double *alpdotvec,\
  double *gimddotvec, double *Phiddotvec, double *alpddotvec, double dt_ph, int Nph, double *par, double nu0, double zeta, int nmodes, double arm, int *i_plunge);

void waveformAKTDI(double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, int nmodes, double zeta, double *par, double nu0);

void GenAKTDI(double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, double e0, double nu0, double M, double mu, double S, double dist, double inc, double gam0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, bool mich, bool traj);

#endif
