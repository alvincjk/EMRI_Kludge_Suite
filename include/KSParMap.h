// EMRI Kludge Suite: Parameter maps

#ifndef _KSPARMAP_H
#define _KSPARMAP_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#undef C
#include <gsl/gsl_blas.h>
#define C (299792458.)
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_spline.h>
#undef G
#include <gsl/gsl_multifit.h>
#define G (6.6726e-11)

#define EulerGamma (0.5772156649015329)
#define PI2 (9.86960440109)

struct sol_par{

  double Omega_r;
  double Omega_theta;
  double Omega_phi;
  double M;
  double e;
  double iota;

};

void cross(const gsl_vector *u,const gsl_vector *v,gsl_vector *w);

void ParAng(double ang[],double e,double iota,double gamma,double psi,double theta_S,double phi_S,double theta_K,double phi_K,double alpha,double zm);

// ----- 3PN O(e^6) equations (Sago & Fujita, 2015) -----
double dvdt(double v,double e,double Y,double m,double M,double q);
double dedt(double v,double e,double Y,double m,double M,double q);
double dtdm(double v,double e,double Y,double q);
double drdm(double v,double e,double Y,double q);
double dthetadm(double v,double e,double Y,double q);
double dphidm(double v,double e,double Y,double q);
double dvdt2(double v,double e,double Y,double m,double M,double q,double dv,double de);
double dedt2(double v,double e,double Y,double m,double M,double q,double dv,double de);
// ----------

int sol_fun(const gsl_vector *x,void *p,gsl_vector *f);

int sol_inv(const gsl_vector *x,void *p,gsl_vector *f);

void print_state(size_t i, gsl_multiroot_fsolver *sol);

void ParMap(double map[],double Omega[],double p,double M,double s,double e,double iota);

void ParInvMap(double map[],double Omega[],double p,double M,double s,double e,double iota);

void Interp(double *x_in,double *y_in,int n_in,double *x_out,double *y_out,int n_out);

void PolyFit(double *coeff,double *x,double *y,int n);

void RotCoeff(double rot[],double iota,double theta_S,double phi_S,double theta_K,double phi_K,double alpha);

#endif

