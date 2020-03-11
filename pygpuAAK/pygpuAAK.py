import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
import tdi

import time

MTSUN = 1.989e30 * ct.G / ct.c ** 3


class Converter:
    def __init__(self, key_order, **kwargs):

        self.conversions = []
        if "ln_mT" in key_order:
            # replace ln_mT and mr with m1 and m2 respectively
            self.ind_ln_mT = key_order.index("ln_mT")
            self.conversions.append(self.ln_mT)

        if "ln_mu" in key_order:
            # replace ln_mT and mr with m1 and m2 respectively
            self.ind_ln_mu = key_order.index("ln_mu")
            self.conversions.append(self.ln_mu)

        if "cos_iota" in key_order:
            self.ind_iota = key_order.index("cos_iota")
            self.conversions.append(self.cos_iota)

        if "sin_theta_S" in key_order:
            self.ind_theta_S = key_order.index("sin_theta_S")
            self.conversions.append(self.sin_theta_S)

        if "sin_theta_K" in key_order:
            self.ind_theta_K = key_order.index("sin_theta_K")
            self.conversions.append(self.sin_theta_K)

        if "ln_distance" in key_order:
            self.ind_ln_distance = key_order.index("ln_distance")
            self.conversions.append(self.ln_distance)

    def ln_mT(self, x):
        x[self.ind_ln_mT] = np.exp(x[self.ind_ln_mT])
        return x

    def ln_mu(self, x):
        x[self.ind_ln_mu] = np.exp(x[self.ind_ln_mu])
        return x

    def ln_distance(self, x):
        x[self.ind_ln_distance] = np.exp(x[self.ind_ln_distance])  # Gpc
        return x

    def cos_iota(self, x):
        x[self.ind_iota] = np.arccos(x[self.ind_iota])
        return x

    def sin_theta_S(self, x):
        x[self.ind_theta_S] = np.arcsin(x[self.ind_theta_S])
        return x

    def sin_theta_K(self, x):
        x[self.ind_theta_K] = np.arcsin(x[self.ind_theta_K])
        return x

    def convert(self, x):
        for func in self.conversions:
            x = func(x)
        return x


class Recycler:
    def __init__(self, key_order, **kwargs):

        # Setup of recycler
        self.recycles = []

        if "phi_S" in key_order:
            # assumes beta is also there
            self.ind_phi_S = key_order.index("phi_S")
            self.recycles.append(self.phi_S)

        if "phi_K" in key_order:
            # assumes beta is also there
            self.ind_phi_K = key_order.index("phi_K")
            self.recycles.append(self.phi_K)

        if "gamma" in key_order:
            self.ind_gamma = key_order.index("gamma")
            self.recycles.append(self.gamma)

        if "psi" in key_order:
            self.ind_psi = key_order.index("psi")
            self.recycles.append(self.psi)

        if "alph" in key_order:
            self.ind_alph = key_order.index("alph")
            self.recycles.append(self.alph)

    def phi_S(self, x):
        x[self.ind_phi_S] = x[self.ind_phi_S] % (2 * np.pi)
        return x

    def phi_K(self, x):
        x[self.ind_phi_K] = x[self.ind_phi_K] % (2 * np.pi)
        return x

    def phiRef(self, x):
        x[self.ind_phiRef] = x[self.ind_phiRef] % (2 * np.pi)
        return x

    def gamma(self, x):
        x[self.ind_gamma] = x[self.ind_gamma] % (2 * np.pi)
        return x

    def psi(self, x):
        x[self.ind_psi] = x[self.ind_psi] % (2 * np.pi)
        return x

    def alph(self, x):
        x[self.ind_alph] = x[self.ind_alph] % (2 * np.pi)
        return x

    def recycle(self, x):
        for func in self.recycles:
            x = func(x)
        return x


class pyGPUAAK:
    def __init__(
        self, data_stream, noise_channels, length, dt, init_dt, key_order, **kwargs
    ):
        """
        data_stream (dict): keys of channel1, channel2
        """
        prop_defaults = {
            "T_fit": 1.0,
            "LISA": True,
            "backint": True,
            "eps": 1e-6,
            "nwalkers": 1,
            "ndevices": 1,
        }

        self.ndim = len(key_order)
        self.length = length
        self.init_dt = init_dt
        self.dt = dt

        self.converter = Converter(key_order)
        self.recycler = Recycler(key_order)

        self.total_time = self.length * self.dt
        self.init_length = int(np.ceil(self.total_time / self.init_dt))
        for prop, default in prop_defaults.items():
            setattr(self, prop, kwargs.get(prop, default))
            # TODO: check this
            kwargs[prop] = kwargs.get(prop, default)

        self.data_channel1 = data_stream["channel1"].copy()
        self.data_channel2 = data_stream["channel2"].copy()

        freqs = np.fft.rfftfreq(self.length, d=dt)
        self.fft_length = len(freqs)
        df = 1.0 / (self.length * self.dt)

        # whiten data stream
        self.data_stream = data_stream

        self.noise_channel1_inv = (
            1.0 / np.sqrt(noise_channels["channel1"]) * np.sqrt(dt / length)
        )  # dt dt df = dt
        self.noise_channel2_inv = (
            1.0 / np.sqrt(noise_channels["channel2"]) * np.sqrt(dt / length)
        )

        self.data_channel1 *= self.noise_channel1_inv
        self.data_channel2 *= self.noise_channel2_inv

        self.d_d = 4 * np.sum(
            [np.abs(self.data_channel1) ** 2, np.abs(self.data_channel2) ** 2]
        )

        self.generator = GPUAAK(
            self.T_fit,
            self.init_length,
            self.length,
            self.init_dt,
            self.dt,
            self.LISA,
            self.backint,
            self.data_channel1,
            self.data_channel2,
            self.noise_channel1_inv,
            self.noise_channel2_inv,
        )

    def NLL(
        self,
        iota,
        s,
        p,
        e,
        M,
        mu,
        gamma,
        psi,
        alph,
        theta_S,
        phi_S,
        theta_K,
        phi_K,
        D,
        return_snr=False,
        return_TDI=False,
        **kwargs,
    ):

        (iota, s, p, e, M, mu, gamma, psi, alph, theta_S, phi_S, theta_K, phi_K, D) = (
            np.atleast_1d(iota),
            np.atleast_1d(s),
            np.atleast_1d(p),
            np.atleast_1d(e),
            np.atleast_1d(M),
            np.atleast_1d(mu),
            np.atleast_1d(gamma),
            np.atleast_1d(psi),
            np.atleast_1d(alph),
            np.atleast_1d(theta_S),
            np.atleast_1d(phi_S),
            np.atleast_1d(theta_K),
            np.atleast_1d(phi_K),
            np.atleast_1d(D),
        )

        d_h, h_h = np.zeros_like(iota), np.zeros_like(iota)

        if return_TDI:
            hI = np.zeros((len(iota), self.fft_length), dtype=np.complex128)
            hII = np.zeros((len(iota), self.fft_length), dtype=np.complex128)

        for i in range(len(iota)):
            d_h[i], h_h[i] = self.generator.WaveformThroughLikelihood(
                iota[i],
                s[i],
                p[i],
                e[i],
                M[i],
                mu[i],
                gamma[i],
                psi[i],
                alph[i],
                theta_S[i],
                phi_S[i],
                theta_K[i],
                phi_K[i],
                D[i],
            )

            if return_TDI:
                _, hI[i], hII[i] = self.generator.GetWaveform(is_Fourier=True)

        if return_TDI:
            return hI, hII

        if return_snr:
            return np.sqrt(h_h)

        return 1 / 2 * (self.d_d + h_h - 2 * d_h)

    def getNLL(self, x, **kwargs):
        # changes parameters to in range in actual array (not copy)
        x = self.recycler.recycle(x)

        # converts parameters in copy, not original array
        x_in = self.converter.convert(x.copy())

        iota, s, p, e, M, mu, gamma, psi, alph, theta_S, phi_S, theta_K, phi_K, D = x_in

        return self.NLL(
            iota,
            s,
            p,
            e,
            M,
            mu,
            gamma,
            psi,
            alph,
            theta_S,
            phi_S,
            theta_K,
            phi_K,
            D,
            **kwargs,
        )

    def get_Fisher(self, x):
        Mij = np.zeros((self.ndim, self.ndim), dtype=x.dtype)
        x_in = np.tile(x, (2 * self.ndim, 1))

        for i in range(self.ndim):
            x_in[2 * i, i] += self.eps
            x_in[2 * i + 1, i] -= self.eps

        hI, hII = self.getNLL(x_in.T, return_TDI=True)

        for i in range(self.ndim):
            hIi_up, hIIi_up = hI[2 * i + 1], hII[2 * i + 1]
            hIi_down, hIIi_down = hI[2 * i], hII[2 * i]

            hi_hI = (hIi_up - hIi_down) / (2 * self.eps)
            hi_hII = (hIIi_up - hIIi_down) / (2 * self.eps)

            for j in range(i, self.ndim):
                hIj_up, hIIj_up = hI[2 * j + 1], hII[2 * j + 1]
                hIj_down, hIIj_down = hI[2 * j], hII[2 * j]

                hj_hI = (hIj_up - hIj_down) / (2 * self.eps)
                hj_hII = (hIIj_up - hIIj_down) / (2 * self.eps)

                inner_product = 4 * np.real(
                    (np.dot(hj_hI.conj(), hj_hI) + np.dot(hi_hII.conj(), hj_hII))
                )

                Mij[i][j] = inner_product
                Mij[j][i] = inner_product

        return Mij
