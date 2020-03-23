import numpy as np

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
    if dur == 0 or dur is None:
        return 0.0
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
    S_n = 10.0 / (3.0 * L ** 2) * (
        P_OMS(f) + 2 * (1 + np.cos(f / f_star) ** 2) * P_acc(f) / (2 * np.pi * f) ** 4
    ) * (1 + 6 / 10 * (f / f_star) ** 2) + S_c(f, dur=dur)
    return S_n


def generate_noise_frequencies(Tobs, fs):
    df = 1.0 / Tobs
    number_of_samples = int(np.round(Tobs * fs))
    number_of_frequencies = int(np.round(number_of_samples / 2) + 1)

    noise_freqs = np.linspace(start=0, stop=fs / 2, num=number_of_frequencies)

    return noise_freqs


def generate_noise_single_channel(noise_func, noise_args, noise_kwargs, df, data_freqs):

    norm1 = 0.5 * (1.0 / df) ** 0.5
    re = np.random.normal(0, norm1, data_freqs.shape)
    im = np.random.normal(0, norm1, data_freqs.shape)
    htilde = re + 1j * im

    return np.sqrt(noise_func(data_freqs, *noise_args, **noise_kwargs)) * htilde
