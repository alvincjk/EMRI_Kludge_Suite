import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
from pygpuAAK import pygpuAAK
import tdi

import time


def test():
    iota = 0.2
    s = 0.8
    p = 10.0
    e = 0.7
    T_fit = 1.0
    init_length = 14000
    length = 1400000
    init_dt = 1000.0
    dt = 10.0
    M = 1e6
    mu = 1e1
    gamma = 0.0
    psi = 0.0
    alph = 0.0
    theta_S = 0.785
    phi_S = 0.785
    theta_K = 1.05
    phi_K = 1.05
    D = 1.0
    LISA = True
    backint = True

    pars = {'backint': True,
        'LISA': True,
        'length': length,
        'dt': dt,
        'p': p,
        'T': 1.,
        'f': 2.e-3,
        'T_fit': T_fit,
        'mu': mu,
        'M': M,
        's': s,
        'e': e,
        'iota': iota,
        'gamma': gamma,
        'psi': psi,
        'theta_S': theta_S,
        'phi_S': phi_S,
        'theta_K': theta_K,
        'phi_K': phi_K,
        'alpha': alph,
        'D': D}

    print('Running cpu waveform.')
    st = time.perf_counter()
    tc, hIc, hIIc, timing = cpu_AAK.wave(pars)
    et = time.perf_counter()
    print('CPU Waveform complete: {} seconds'.format(et - st))

    hI_f, hII_f = np.fft.rfft(hIc), np.fft.rfft(hIIc)
    data_stream = {'channel1': hI_f, 'channel2': hII_f}

    freqs = np.fft.rfftfreq(len(hIc), d=dt)
    freqs[0] = 1e-8
    deltaF = 1/dt
    AE_ASD = tdi.noisepsd_AE(freqs, model='SciRDv1', includewd=3)
    noise_channels = {'channel1': AE_ASD, 'channel2': AE_ASD}

    kwargs = {
        'T_fit': 1.0,
        'LISA': True,
        'backint': True
    }

    like_class =  pygpuAAK.pyGPUAAK(data_stream, noise_channels, length, dt, init_dt, **kwargs)
    likelihood = like_class.NLL(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)
    num = 1000
    st = time.perf_counter()
    for i in range(num):
        likelihood = like_class.NLL(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                                phi_S, theta_K, phi_K, D)
    et = time.perf_counter()
    print('{} likelihood calculations with {} time points.'.format(num, length))
    print('Time per likelihood calculation:', (et - st)/num)

    snr = like_class.NLL(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D, return_snr=True)

    print('SNR: ', snr, 'Likelihood:', likelihood)

if __name__ == "__main__":
    test()
