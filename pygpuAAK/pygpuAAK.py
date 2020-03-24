import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
from pygpuAAK.LISAnoise import (
    LISA_Noise,
    generate_noise_frequencies,
    generate_noise_single_channel,
)
from pygpuAAK.convert import Converter

import tdi

import time

MTSUN = 1.989e30 * ct.G / ct.c ** 3


class pyGPUAAK:
    def __init__(self, data_info, length, dt, init_dt, key_order, **kwargs):
        """
        data_info (dict): keys of channel1, channel2 or injection params
        """
        prop_defaults = {
            "T_fit": 1.0,
            "LISA": True,
            "backint": True,
            "eps": 1e-6,
            "nwalkers": 1,
            "ndevices": 1,
            "wd_dur": 1,
            "add_noise": None,
        }

        self.ndim = len(key_order)
        self.length = length
        self.init_dt = init_dt
        self.dt = dt

        self.converter = Converter("emri", key_order)
        self.key_order = key_order

        self.total_time = self.length * self.dt
        self.init_length = int(np.ceil(self.total_time / self.init_dt))
        for prop, default in prop_defaults.items():
            setattr(self, prop, kwargs.get(prop, default))
            # TODO: check this
            kwargs[prop] = kwargs.get(prop, default)

        self.generator = GPUAAK(
            self.T_fit,
            self.init_length,
            self.length,
            self.init_dt,
            self.dt,
            self.LISA,
            self.backint,
        )

        self.determine_freqs_noise(**kwargs)

        self.data_info = data_info

        self.create_input_data(**kwargs)

    def determine_freqs_noise(self, **kwargs):
        self.freqs = generate_noise_frequencies(self.length * self.dt, 1 / self.dt)
        self.fft_length = len(self.freqs)
        df = 1.0 / (self.length * self.dt)

        self.freqs[0] = self.freqs[1] / 10.0

        self.added_noise = [0.0, 0.0]
        if self.add_noise is not None:
            self.added_noise[0] = generate_noise_single_channel(
                LISA_Noise, [], {"dur": self.wd_dur}, df, self.freqs
            )

            self.added_noise[1] = generate_noise_single_channel(
                LISA_Noise, [], {"dur": self.wd_dur}, df, self.freqs
            )

    def create_input_data(self, **kwargs):

        self.noise_channel1_inv = (
            1.0
            / np.sqrt(LISA_Noise(self.freqs, dur=self.wd_dur))
            * np.sqrt(self.dt / self.length)
        )  # dt dt df = dt

        self.noise_channel2_inv = (
            1.0
            / np.sqrt(LISA_Noise(self.freqs, dur=self.wd_dur))
            * np.sqrt(self.dt / self.length)
        )  # dt dt df = dt

        if "channel1" in self.data_info:

            self.injection_channel1 = self.data_info["channel1"].copy()
            self.injection_channel2 = self.data_info["channel2"].copy()

        else:
            self.injection_params = self.data_info["injection_params"]

            data_channel_temp = np.zeros_like(
                self.noise_channel1_inv, dtype=np.complex128
            )

            noise_channel_inv_temp = np.ones_like(
                self.noise_channel1_inv, dtype=np.float64
            )

            self.generator.input_data(
                data_channel_temp,
                data_channel_temp,
                noise_channel_inv_temp,
                noise_channel_inv_temp,
            )

            injection_arr = np.tile(
                np.array([self.injection_params.get(key) for key in self.key_order]),
                (1, 1),
            )

            hI_f, hII_f = self.getNLL(injection_arr.T, return_TDI=True)

            self.injection_channel1_no_noise = hI_f
            self.injection_channel2_no_noise = hII_f
            self.injection_channel1 = hI_f[0] + self.added_noise[0]
            self.injection_channel2 = hII_f[0] + self.added_noise[1]

        self.data_channel1 = self.noise_channel1_inv * self.injection_channel1
        self.data_channel2 = self.noise_channel2_inv * self.injection_channel2

        self.generator.input_data(
            self.data_channel1,
            self.data_channel2,
            self.noise_channel1_inv,
            self.noise_channel2_inv,
        )

        self.d_d = 4 * np.sum(
            [np.abs(self.data_channel1) ** 2, np.abs(self.data_channel2) ** 2]
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
        return_waveform=False,
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
            factor = 1
            dtype = np.complex128
            is_Fourier = True

            hI = np.zeros((len(iota), self.fft_length), dtype=dtype)
            hII = np.zeros((len(iota), self.fft_length), dtype=dtype)

        if return_waveform:
            factor = 2
            dtype = np.float64
            is_Fourier = False

            hI = np.zeros((len(iota), self.length + 2), dtype=dtype)
            hII = np.zeros((len(iota), self.length + 2), dtype=dtype)

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
                return_waveform=return_waveform,
            )

            if return_TDI or return_waveform:
                _, temp1, temp2 = self.generator.GetWaveform(is_Fourier=is_Fourier)

                if is_Fourier:
                    hI[i].real, hI[i].imag = temp1[0::2], temp1[1::2]
                    hII[i].real, hII[i].imag = temp2[0::2], temp2[1::2]

                else:
                    hI[i] = temp1  # [:-2]
                    hII[i] = temp2  # [:-2]

        if return_TDI or return_waveform:
            return hI, hII

        if return_snr:
            return np.sqrt(d_h), np.sqrt(h_h)

        return 1 / 2 * (self.d_d + h_h - 2 * d_h)

    def getNLL(self, x, **kwargs):
        # changes parameters to in range in actual array (not copy)
        x = self.converter.recycle(x)

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
