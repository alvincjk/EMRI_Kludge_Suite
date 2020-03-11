import AAKwrapper


##### GENERATE AAK WAVEFORM #####

pars = {
    "backint": True,
    "LISA": False,
    "length": 1000000,
    "dt": 10.0,
    "p": 6.5,
    "T": 1.0,
    "f": 2.0e-3,
    "T_fit": 1.0,
    "mu": 1.0e1,
    "M": 1.0e6,
    "s": 0.5,
    "e": 0.2,
    "iota": 0.5,
    "gamma": 0.0,
    "psi": 0.0,
    "theta_S": 1.05,
    "phi_S": 0.0,
    "theta_K":0.8,
    "phi_K": 0.8,
    "alpha": 0.0,
    "D": 1.8,
}

print("Computing waveform...")

t, hI, hII, timing = AAKwrapper.wave(pars)
import pdb

pdb.set_trace()
print("Time taken: {}".format(timing))


##### GENERATE AAK TDIS #####

pars = {
    "backint": False,
    "LISA": False,  # NOT USED
    "length": 1000000,
    "dt": 10.0,
    "p": 8.0,
    "T": 1.0,
    "f": 2.0e-3,
    "T_fit": 1.0,
    "mu": 1.0e1,
    "M": 1.0e6,
    "s": 0.5,
    "e": 0.1,
    "iota": 0.524,
    "gamma": 0.0,
    "psi": 0.0,
    "theta_S": 0.785,
    "phi_S": 0.785,
    "theta_K": 1.05,
    "phi_K": 1.05,
    "alpha": 0.0,
    "D": 1.0,
}

print("Computing TDIs (AAK)...")

f, Xf_r, Xf_im, Yf_r, Yf_im, Zf_r, Zf_im, timing = AAKwrapper.tdi(pars)

print("Time taken: {}".format(timing))


##### GENERATE AAK PHASES #####

pars = {
    "backint": True,
    "LISA": False,  # NOT USED
    "length": 1000,
    "dt": 63072.0,
    "p": 6.0,
    "T": 1.0,
    "f": 2.0e-3,
    "T_fit": 1.0,
    "mu": 1.0e1,
    "M": 1.0e6,
    "s": 0.5,
    "e": 0.1,
    "iota": 0.524,
    "gamma": 0.0,
    "psi": 0.0,
    "theta_S": 0.785,
    "phi_S": 0.785,
    "theta_K": 1.05,
    "phi_K": 1.05,
    "alpha": 0.0,
    "D": 1.0,
}

print("Computing phases...")

t, phase_r, phase_theta, phase_phi, omega_r, omega_theta, omega_phi, eccentricity, timing = AAKwrapper.phase(
    pars
)

print("Time taken: {}".format(timing))


##### GENERATE AK TDIS #####

pars = {
    "backint": False,  # NOT USED
    "LISA": False,  # NOT USED
    "length": 1000000,
    "dt": 5.184,
    "p": 8.0,
    "T": 1.0,
    "f": 2.0e-3,
    "T_fit": 1.0,
    "mu": 1.0e1,
    "M": 1.0e6,
    "s": 0.5,
    "e": 0.1,
    "iota": 0.524,
    "gamma": 0.0,
    "psi": 0.0,
    "theta_S": 0.785,
    "phi_S": 0.785,
    "theta_K": 1.05,
    "phi_K": 1.05,
    "alpha": 0.0,
    "D": 1.0,
}

print("Computing TDIs (AK)...")

f, Xf_r, Xf_im, Yf_r, Yf_im, Zf_r, Zf_im, timing = AAKwrapper.aktdi(pars)

print("Time taken: {}".format(timing))
