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

%pythoncode %{
import numpy

def AAK(pars = {}):
    setpar = SetPar()

    setpar.backint = True
    setpar.LISA    = False
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

    swig_AAKwave(setpar, t, hI, hII)

    return t, hI, hII
%}
