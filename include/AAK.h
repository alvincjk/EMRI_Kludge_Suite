// AAK: Waveform

#ifndef _AAK_H
#define _AAK_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Constants.h"

double OmegaPhi(double v, double e, double cosiota, double s, double M);

void PNevolution(double *t, double *e, double *v, double *M, double *S, double *gim, double *Phi, double *alp, double *nu, double *gimdotvec, double timestep, int vlength, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, int *i_plunge, int *i_buffer, bool backint);

void waveform(double *t, double *hI, double *hII, double timestep, int vlength, int nmodes, double zeta, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, bool backint, bool mich, bool traj);

void GenWave(double *t, double *hI, double *hII, double timestep, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_phys, double S_map[], double dist, double inc, double gim0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, double dt_map[], int steps, bool backint, bool mich, bool traj);

#endif
