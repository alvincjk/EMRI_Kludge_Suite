#ifndef _AAKPY_H
#define _AAKPY_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAK.h"
#include "AAKPhase.h"
#include "AAKTDI.h"
#include "AKTDI.h"

double AAKwave(SetPar &AAK, double *t, double *hI, double *hII);

double AAKphase(SetPar &AAK, double *t, double *phase_r, double *phase_theta, double *phase_phi, double *omega_r, double *omega_theta, double *omega_phi, double *eccentricity);

double AAKTDI(SetPar &AAK, double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im);

double AKTDI(SetPar &AK, double *f, double *Xf_r, double *Xf_im, double *Yf_r, double *Yf_im, double *Zf_r, double *Zf_im);

#endif
