%module AAKwrapper

// these are the C includes needed in the interface C++ file
%{
#include "AAKpy.h"
%}

// these are the numpy includes that will be wrapper
%include "include/KSTools.h"

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{ import_array(); %}

// see https://docs.scipy.org/doc/numpy/reference/swig.interface-file.html, 
// and especially discussion "Beyond the provided typemaps"

%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* vec1),(int len2, double* vec2),(int len3, double* vec3)}
%exception swig_AAKwave {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
double swig_AAKwave(struct SetPar* setpar, int len1, double* vec1, int len2, double* vec2, int len3, double* vec3) {
    if (len1 != setpar->length || len2 != setpar->length || len3 != setpar->length) {
        PyErr_Format(PyExc_ValueError, "Arrays of lengths (%d,%d,%d) given, need %d", len1, len2, len3, setpar->length);
        return 0.0;
    }
    return AAKwave(*setpar, vec1, vec2, vec3);
}
%}

%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* vec1),(int len2, double* vec2),(int len3, double* vec3),(int len4, double* vec4),(int len5, double* vec5),(int len6, double* vec6),(int len7, double* vec7)}
%exception swig_AAKphase {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
double swig_AAKphase(struct SetPar* setpar, int len1, double* vec1, int len2, double* vec2, int len3, double* vec3, int len4, double* vec4, int len5, double* vec5, int len6, double* vec6, int len7, double* vec7) {
    if (len1 != setpar->length || len2 != setpar->length || len3 != setpar->length || len4 != setpar->length || len5 != setpar->length || len6 != setpar->length || len7 != setpar->length) {
        PyErr_Format(PyExc_ValueError, "Arrays of lengths (%d,%d,%d,%d,%d,%d,%d) given, need %d", len1, len2, len3, len4, len5, len6, len7, setpar->length);
        return 0.0;
    }
    return AAKphase(*setpar, vec1, vec2, vec3, vec4, vec5, vec6, vec7);
}
%}

%pythoncode %{
import numpy

def wave(pars = {}):
    setpar = SetPar()

    setpar.backint = pars['backint']
    setpar.LISA    = pars['LISA']
    setpar.traj    = False
    setpar.SNR     = False
    setpar.timing  = False

    setpar.length = pars['length']

    setpar.dt = pars['dt']
    setpar.p = pars['p']
    setpar.T = pars['T']
    setpar.f = pars['f']
    setpar.T_fit = pars['T_fit']

    setpar.mu = pars['mu']
    setpar.M = pars['M']
    setpar.s = pars['s']
    setpar.e = pars['e']
    setpar.iota = pars['iota']
    setpar.gamma = pars['gamma']
    setpar.psi = pars['psi']

    setpar.theta_S = pars['theta_S']
    setpar.phi_S = pars['phi_S']
    setpar.theta_K = pars['theta_K']
    setpar.phi_K = pars['phi_K']
    setpar.alpha = pars['alpha']
    setpar.D = pars['D']

    t   = numpy.zeros(setpar.length, 'd')
    hI  = numpy.zeros(setpar.length, 'd')
    hII = numpy.zeros(setpar.length, 'd')

    timing = swig_AAKwave(setpar, t, hI, hII)

    return t, hI, hII, timing
%}

%pythoncode %{
import numpy

def phase(pars = {}):
    setpar = SetPar()

    setpar.backint = pars['backint']
    setpar.LISA    = pars['LISA']
    setpar.traj    = False
    setpar.SNR     = False
    setpar.timing  = False

    setpar.length = pars['length']

    setpar.dt = pars['dt']
    setpar.p = pars['p']
    setpar.T = pars['T']
    setpar.f = pars['f']
    setpar.T_fit = pars['T_fit']

    setpar.mu = pars['mu']
    setpar.M = pars['M']
    setpar.s = pars['s']
    setpar.e = pars['e']
    setpar.iota = pars['iota']
    setpar.gamma = pars['gamma']
    setpar.psi = pars['psi']

    setpar.theta_S = pars['theta_S']
    setpar.phi_S = pars['phi_S']
    setpar.theta_K = pars['theta_K']
    setpar.phi_K = pars['phi_K']
    setpar.alpha = pars['alpha']
    setpar.D = pars['D']

    t           = numpy.zeros(setpar.length, 'd')
    phase_r     = numpy.zeros(setpar.length, 'd')
    phase_theta = numpy.zeros(setpar.length, 'd')
    phase_phi   = numpy.zeros(setpar.length, 'd')
    omega_r     = numpy.zeros(setpar.length, 'd')
    omega_theta = numpy.zeros(setpar.length, 'd')
    omega_phi   = numpy.zeros(setpar.length, 'd')

    timing = swig_AAKphase(setpar, t, phase_r, phase_theta, phase_phi, omega_r, omega_theta, omega_phi)

    return t, phase_r, phase_theta, phase_phi, omega_r, omega_theta, omega_phi, timing
%}
