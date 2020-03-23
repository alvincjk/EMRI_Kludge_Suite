import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
from pygpuAAK import pygpuAAK

import time


def test():

    iota = 2.142199999999910620e00
    s = 9.697000000000000e-01
    p = 13.0  # 1.081981106649140578e01
    e = 2.286566522026051151e-01
    T_fit = 1.0
    init_length = 100000
    length = 4194304
    init_dt = 100.0
    dt = 15.0
    M = 1.134944869275921956e06
    mu = 2.948999954778721033e01
    gamma = 5.659686260028827576e00
    psi = 2.391562081325919742e00
    alph = 1.175479042396851526e00
    theta_S = 1.071851837554504527e00
    phi_S = 2.232796975919999927e00
    theta_K = np.pi - 1.522100180545759907e00
    phi_K = 3.946698166040000011e00
    D = 1.0
    LISA = True
    backint = True

    """
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
    """

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
        "ln_D",
    ]

    """
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
    """

    kwargs = {"T_fit": 1.0, "LISA": True, "backint": True, "wd_dur": 0}

    injection_vals = [
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
    injection = {key: value for (key, value) in zip(key_order, injection_vals)}

    data_info = {"injection_params": injection}
    like_class = pygpuAAK.pyGPUAAK(data_info, length, dt, init_dt, key_order, **kwargs)
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

    e = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    x[3] = e

    likelihood = like_class.getNLL(x)

    fisher = like_class.get_Fisher(x.T[0])

    num = 100
    st = time.perf_counter()
    for i in range(num):
        likelihood = like_class.getNLL(x)
    et = time.perf_counter()
    print(
        "{} likelihood calculations with {} time points.".format(
            num * x.shape[1], length
        )
    )
    print("Time per likelihood calculation:", (et - st) / (num * x.shape[1]))

    snr = like_class.getNLL(x, return_snr=True)
    print("SNR: ", snr, "Likelihood:", likelihood)


if __name__ == "__main__":
    test()
