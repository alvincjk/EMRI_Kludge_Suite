// AAK: Phases and frequencies

#ifndef _AAKPHASE_H
#define _AAKPHASE_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Constants.h"
#include "AAK.h"

void PNevolutionPhase(double *t, double *e, double *v, double *M, double *S, double *gim, double *Phi, double *alp, double *nu, double *gimdotvec, double *Phidotvec, double *alpdotvec, double timestep, int vlength, double *par, double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[], int steps, int *i_plunge, int *i_buffer, bool backint);

void GenPhase(double *t, double *phase_r, double *phase_theta, double *phase_phi, double *omega_r, double *omega_theta, double *omega_phi, double timestep, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_phys, double S_map[], double dist, double inc, double gim0, double Phi0, double qS, double phiS, double alp0, double qK, double phiK, double dt_map[], int steps, bool backint, bool mich, bool traj);

#endif
