// AAK: Waveform

#ifndef _AAK_H
#define _AAK_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Constants.h"

double OmegaPhi(double v, double e, double cosiota, double s, double M);

void PNevolution(double *t_in, double *e_in, double *v_in, double *M_in, double *S_in,
                double timestep, int vlength, double *par,
                double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[],
                int steps, int *i_plunge, int *i_buffer, bool backint, double *t_clip, double *interp_timestep, int *init_length);

#endif
