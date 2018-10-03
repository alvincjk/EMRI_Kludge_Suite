import AAKwrapper

print("Loading AAK from {}.".format(AAKwrapper.__path__))

pars = {'backint': True,
        'LISA': False,
        'length': 1000000,
        'dt': 5.184,
        'p': 6.,
        'T': 1.,
        'f': 2.e-3,
        'T_fit': 1.,
        'mu': 1.e1,
        'M': 1.e6,
        's': 0.5,
        'e': 0.1,
        'iota': 0.524,
        'gamma': 0.,
        'psi': 0.,
        'theta_S': 0.785,
        'phi_S': 0.785,
        'theta_K': 1.05,
        'phi_K': 1.05,
        'alpha': 0.,
        'D': 1.}

print("Computing waveform...")

t, hI, hII, timing = AAKwrapper.wave(pars)

print("Time taken: {}".format(timing))

pars = {'backint': True,
        'LISA': False,
        'length': 1000,
        'dt': 63072.,
        'p': 6.,
        'T': 1.,
        'f': 2.e-3,
        'T_fit': 1.,
        'mu': 1.e1,
        'M': 1.e6,
        's': 0.5,
        'e': 0.1,
        'iota': 0.524,
        'gamma': 0.,
        'psi': 0.,
        'theta_S': 0.785,
        'phi_S': 0.785,
        'theta_K': 1.05,
        'phi_K': 1.05,
        'alpha': 0.,
        'D': 1.}

print("Computing phases...")

t, phase_r, phase_theta, phase_phi, omega_r, omega_theta, omega_phi, timing = AAKwrapper.phase(pars)

print("Time taken: {}".format(timing))
