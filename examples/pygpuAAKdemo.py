import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
from pygpuAAK import pygpuAAK

import time

fit_pars = {
    0.5:{
        'alpha': 0.133,
        'beta': 243,
        'kappa': 482,
        'gamma': 917,
        'f_k': 2.58e-3
    },
    1:{
        'alpha': 0.171,
        'beta': 292,
        'kappa': 1020,
        'gamma': 1680,
        'f_k': 2.15e-3
    },
    2:{
        'alpha': 0.165,
        'beta': 299,
        'kappa': 611,
        'gamma': 1340,
        'f_k': 1.73e-3
    },
    4:{
        'alpha': 0.138,
        'beta': -221,
        'kappa': 521,
        'gamma': 1680,
        'f_k': 1.13e-3
    }
}

def P_OMS(f):
    return (1.5e-11)**2*(1 + (2e-3/f)**4)

def P_acc(f):
    return (3e-15)**2*(1 + (0.4e-3/f)**2)*(1+(f/8e-3)**4)

def S_c(f, dur=4):
    if dur not in [0.5, 1, 2, 4]:
        raise ValueError("dur needs to be 0.5, 1, 2, or 4 years.")

    alpha = fit_pars[dur]['alpha']
    beta = fit_pars[dur]['beta']
    kappa = fit_pars[dur]['kappa']
    gamma = fit_pars[dur]['gamma']
    f_k = fit_pars[dur]['f_k']
    A = 9e-45
    return A* f**(-7/3)*np.exp(-f**alpha + beta*f*np.sin(kappa*f))*(1+np.tanh(gamma*(f_k - f)))


def LISA_Noise(f, L=2.5e9, f_star=19.09e-3, dur=4):
    S_n=20./(3.*L**2)*(P_OMS(f) + 4*P_acc(f)/(2*np.pi*f)**4)*(1+6/10*(f/f_star)**2)+S_c(f, dur=4)
    return S_n;


def test():
    iota = 0.2
    s = 0.8
    p = 10.0
    e = 0.7
    T_fit = 1.0
    init_length = 100000
    length = 1000000
    init_dt = 100.0
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
    D = 5.0
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
    ASD = np.sqrt(LISA_Noise(freqs, dur=4))
    noise_channels = {'channel1': ASD**2, 'channel2': ASD**2}

    kwargs = {
        'T_fit': 1.0,
        'LISA': True,
        'backint': True
    }

    like_class =  pygpuAAK.pyGPUAAK(data_stream, noise_channels, length, dt, init_dt, **kwargs)
    likelihood = like_class.NLL(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)
    num = 10
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
    import pdb; pdb.set_trace()

if __name__ == "__main__":
    test()
