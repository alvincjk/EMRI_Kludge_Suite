from libcpp cimport bool

import numpy
cimport numpy

cdef extern from "KSTools.h":
    ctypedef struct SetPar:
        char path[100]

        bool backint
        bool LISA
        bool traj
        bool SNR
        bool timing

        int length      # waveform points
        double dt       # time step in seconds
        double p        # initial semi-latus rectum in M
        double T        # waveform duration in years
        double f        # initial GW reference frequency in Hz
        double T_fit    # duration of local fit in radiation-reaction time steps M^2/mu

        double mu       # CO mass in solar masses
        double M        # BH mass in solar masses
        double s        # spin parameter (=a/M=S/M^2)
        double e        # initial eccentricity
        double iota     # inclination angle of L from S
        double gamma    # initial angle of periapsis from LxS
        double psi      # initial true anomaly

        double theta_S  # initial source polar angle in ecliptic coordinates
        double phi_S    # initial source azimuthal angle in ecliptic coordinates
        double theta_K  # initial BH spin polar angle in ecliptic coordinates
        double phi_K    # initial BH spin azimuthal angle in ecliptic coordinates
        double alpha    # initial azimuthal orientation
        double D        # source distance in Gpc

cdef extern from "AAKwave.h":
    double AAKwave(SetPar& AAK, double *t, double *hI, double *hII)

def AAKdemo():
    pars = {'length': 1000000, 'dt': 5., 'p': 5., 'T': 1., 'f': 2.e-3, 'T_fit': 1., 'mu': 1.e1, 'M': 1.e6, 's': 0.5, 'e': 0.1, 'iota': 0.524, 'gamma': 0., 'psi': 0., 'theta_S': 0.785, 'phi_S': 0.785, 'theta_K': 1.05, 'phi_K': 1.05, 'alpha': 0., 'D': 1.}

    return AAK(pars)

def AAK(pars = {}):
    cdef SetPar setpar

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

    cdef numpy.ndarray[double,ndim=1] t   = numpy.zeros(setpar.length, 'd')
    cdef numpy.ndarray[double,ndim=1] hI  = numpy.zeros(setpar.length, 'd')
    cdef numpy.ndarray[double,ndim=1] hII = numpy.zeros(setpar.length, 'd')

    AAKwave(setpar,&t[0], &hI[0], &hII[0])

    return t, hI, hII
