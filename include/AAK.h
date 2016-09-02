// AAK: Waveform

#ifndef _AAK_H
#define _AAK_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G (6.6726e-11)
#define C (299792458.)
#define Msun (1.9889e30)
#define SOLARMASSINSEC (G*Msun/(C*C*C))
#define Gpc (1.02938e17)
#define GPCINSEC (3.08567818585e25/C)
#define nt (620000)
#define ndim (12)
#define modes (40)

double J0(double x);

double J1(double x);

double Jn(int n, double x);

void PNevolution(int vlength, double timestep, double *par, double v_map[], double *gimdotvec, double *e, double *nu, double *Phi,
		 double *gim, double *alp, double *v, double *M, double *S, double e_traj[], double M_map[], double S_map[], double dt_map);

void waveform(double tend,double *par, double v_map[], int vlength, double timestep, double *hI, double *hII, int nmodes, double zeta, double e_traj[],
	      double M_phys, double M_map[], double S_map[], double dt_map, bool mich, bool traj);

void GenBCWave(double *hI, double *hII, double deltat, int vlength, double e_traj[], double v_map[], double M_phys, double M_map[], double mu, double S_map[], double dist, double inc, double gam0, double Phi0,
	       double qS, double phiS, double alp0, double qK, double phiK, double dt_map, bool mich, bool traj);

#endif
