#ifndef __KERNEL_HH__
#define __KERNEL_HH__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cuComplex.h>
#include "interpolate.hh"
#include "Globals.h"

void cumsum(double *data, double phase0, int n);

__global__
void produce_phasing(double *e_out, double *v_out, double *M_out, double *S_out, double *gimdot_out, double *nu_out, double *alpdot_out,
                    double *gim_out, double *Phi_out, double *alp_out,
                     double *tvec, InterpArrayContainer evec, InterpArrayContainer vvec, InterpArrayContainer Mvec, InterpArrayContainer Svec,
                            double lam,
                            int init_length,
                             double init_dt, double timestep, double t_clip, int run_length);

__global__
void kernel_create_waveform(double *t, double *hI, double *hII,
                         double *tvec, double *evec, double *vvec, double *Mvec, double *Svec,
                         double *gimvec, double *Phivec, double *alpvec,
                         double *nuvec, double *gimdotvec, double lam,
                         double qS, double phiS, double qK, double phiK,
                         bool mich, int init_length, int vlength,int nmodes,
                         int i_plunge, int i_buffer, double zeta, double M_phys,
                         double init_dt, double timestep, int run_length);


__global__
void likelihood_prep(cuDoubleComplex *template_channel1, cuDoubleComplex *template_channel2, double *noise_channel1_inv, double *noise_channel2_inv, int length);

#endif //__KERNEL_HH__
