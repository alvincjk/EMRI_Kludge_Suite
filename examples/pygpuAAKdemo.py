import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
from pygpuAAK import pygpuAAK

import time


def test():
    iota = 0.5
    s = 0.5
    p = 7.0
    e = 0.5
    T_fit = 1.0
    init_length = 100000
    length = 4194304
    init_dt = 100.0
    dt = 15.0
    M = 1134944.869275098
    mu = 100.489999547765798
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
        "ln_distance",
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

    kwargs = {"T_fit": 1.0, "LISA": True, "backint": True}

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
                np.log(M) * 0.99999999,
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

    num = 1
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
    fisher = like_class.get_Fisher(x.T[0])
    print("SNR: ", snr, "Likelihood:", likelihood)


if __name__ == "__main__":
    test()
