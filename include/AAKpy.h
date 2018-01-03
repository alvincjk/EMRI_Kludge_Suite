#ifndef _AAKPY_H
#define _AAKPY_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAK.h"

double AAKwave(SetPar &AAK, double *t, double *hI, double *hII);

#endif