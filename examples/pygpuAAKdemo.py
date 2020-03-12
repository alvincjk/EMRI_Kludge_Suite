import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
from pygpuAAK import pygpuAAK

import time

fit_pars = {
    0.5: {"alpha": 0.133, "beta": 243, "kappa": 482, "gamma": 917, "f_k": 2.58e-3},
    1: {"alpha": 0.171, "beta": 292, "kappa": 1020, "gamma": 1680, "f_k": 2.15e-3},
    2: {"alpha": 0.165, "beta": 299, "kappa": 611, "gamma": 1340, "f_k": 1.73e-3},
    4: {"alpha": 0.138, "beta": -221, "kappa": 521, "gamma": 1680, "f_k": 1.13e-3},
}


def P_OMS(f):
    return (1.5e-11) ** 2 * (1 + (2e-3 / f) ** 4)


def P_acc(f):
    return (3e-15) ** 2 * (1 + (0.4e-3 / f) ** 2) * (1 + (f / 8e-3) ** 4)


def S_c(f, dur=4):
    if dur not in [0.5, 1, 2, 4]:
        raise ValueError("dur needs to be 0.5, 1, 2, or 4 years.")

    alpha = fit_pars[dur]["alpha"]
    beta = fit_pars[dur]["beta"]
    kappa = fit_pars[dur]["kappa"]
    gamma = fit_pars[dur]["gamma"]
    f_k = fit_pars[dur]["f_k"]
    A = 9e-45
    return (
        A
        * f ** (-7 / 3)
        * np.exp(-(f ** alpha) + beta * f * np.sin(kappa * f))
        * (1 + np.tanh(gamma * (f_k - f)))
    )


def LISA_Noise(f, L=2.5e9, f_star=19.09e-3, dur=4):
    S_n = 20.0 / (3.0 * L ** 2) * (P_OMS(f) + 4 * P_acc(f) / (2 * np.pi * f) ** 4) * (
        1 + 6 / 10 * (f / f_star) ** 2
    ) + S_c(f, dur=4)
    return S_n


def test():
    iota = 0.5
    s = 0.9697
    p = 7.0
    e = 0.22865665220266215
    T_fit = 1.0
    init_length = 100000
    length = 4194304
    init_dt = 100.0
    dt = 15.0
    M = 1134944.869275098
    mu = 29.489999547765798
    gamma = 0.0
    psi = 0.0
    alph = 1.1
    theta_S = np.pi - 0.49894448924039203
    phi_S = 2.23279697592
    theta_K = 0.8
    phi_K = 0.8
    D = 1.8
    LISA = True
    backint = True

    pars = {
        "backint": backint,
        "LISA": LISA,
        "length": length + 2,
        "dt": dt,
        "p": p,
        "T": T_fit,
        "f": 2.0e-3,
        "T_fit": T_fit,
        "mu": mu,
        "M": M,
        "s": s,
        "e": e,
        "iota": iota,
        "gamma": gamma,
        "psi": psi,
        "theta_S": theta_S,
        "phi_S": phi_S,
        "theta_K": theta_K,
        "phi_K": phi_K,
        "alpha": alph,
        "D": D,
    }

    key_order = [
        "cos_iota",
        "s",
        "p",
        "e",
        "ln_mT",
        "ln_mu",
        "gamma",
        "psi",
        "alph",
        "sin_theta_S",
        "phi_S",
        "sin_theta_K",
        "phi_K",
        "ln_distance",
    ]

    print("Running cpu waveform.")
    st = time.perf_counter()
    tc, hIc, hIIc, timing = cpu_AAK.wave(pars)
    # tc = np.arange(0.0, length * dt)
    # hIc = np.zeros_like(tc)
    # hIIc = np.zeros_like(tc)
    et = time.perf_counter()
    print("CPU Waveform complete: {} seconds".format(et - st))

    pad = 0
    hIc = np.pad(hIc, (0, pad), "constant")
    hIIc = np.pad(hIIc, (0, pad), "constant")
    hI_f, hII_f = np.fft.rfft(hIc), np.fft.rfft(hIIc)
    data_stream = {"channel1": hI_f, "channel2": hII_f}

    freqs = np.fft.rfftfreq(len(hIc), d=dt)
    freqs[0] = 1e-8
    deltaF = 1 / dt
    ASD = np.sqrt(LISA_Noise(freqs, dur=4))
    noise_channels = {"channel1": ASD ** 2, "channel2": ASD ** 2}

    kwargs = {"T_fit": 1.0, "LISA": True, "backint": True}

    like_class = pygpuAAK.pyGPUAAK(
        data_stream, noise_channels, length + pad, dt, init_dt, key_order, **kwargs
    )
    x = np.tile(
        np.array(
            [
                np.cos(iota),
                s,
                p,
                e,
                np.log(M),
                np.log(mu),
                gamma,
                psi,
                alph,
                np.sin(theta_S),
                phi_S,
                np.sin(theta_K),
                phi_K,
                np.log(D),
            ]
        ),
        (10, 1),
    ).T

    likelihood = like_class.getNLL(x)
    hI, hII = like_class.getNLL(x, return_waveform=True)
    import pdb

    pdb.set_trace()
    num = 1
    st = time.perf_counter()
    for i in range(num):
        likelihood = like_class.getNLL(x)
    et = time.perf_counter()
    print(
        "{} likelihood calculations with {} time points.".format(
            num * x.shape[1], length + pad
        )
    )
    print("Time per likelihood calculation:", (et - st) / (num * x.shape[1]))

    snr = like_class.getNLL(x, return_snr=True)
    fisher = like_class.get_Fisher(x.T[0])
    import pdb

    pdb.set_trace()
    print("SNR: ", snr, "Likelihood:", likelihood)


if __name__ == "__main__":
    test()
