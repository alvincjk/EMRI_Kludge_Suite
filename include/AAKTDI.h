// AAK: TDIs

#ifndef _AAKTDI_H
#define _AAKTDI_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_trig.h>

#include "Constants.h"

using namespace std;

double OmegaPhiTDI(double v, double e, double cosiota, double s, double M);

void PNevolutionTDI(complex<double> **Xp, complex<double> **Xc, complex<double> **Yp, complex<double> **Yc, complex<double> **Zp, complex<double> **Zc,\
 double *Ampvec, double *tvec, double *evec, double *vvec, double *Mvec, double *Svec, double *gimvec, double *Phivec, double *alpvec,\
  double *nuvec, double *gimdotvec, double *alpdotvec, double *gimddotvec, double *Phiddotvec, double *alpddotvec, double timestep, int vlength,\
   double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps,\
    double zeta, int nmodes, double arm, int *i_plunge, bool backint);

void waveformTDI(double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, int nmodes, double zeta, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, bool backint);

void GenTDI(double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im, double timestep, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_phys, double S_map[], double dist, double inc, double gam0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, double dt_map[], int steps, bool backint, bool mich, bool traj);

#endif
