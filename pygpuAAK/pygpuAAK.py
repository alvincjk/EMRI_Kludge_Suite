import numpy as np
from scipy import constants as ct

from AAKwrapper import AAKwrapper as cpu_AAK
from gpuAAK import GPUAAK
import tdi

import time

MTSUN = 1.989e30*ct.G/ct.c**3


class pyGPUAAK:
    def __init__(self, data_stream, noise_channels, length, dt, init_dt, **kwargs):
        """
        data_stream (dict): keys of channel1, channel2
        """
        prop_defaults = {
            'T_fit': 1.0,
            'LISA': True,
            'backint': True
        }

        self.length = length
        self.init_dt = init_dt
        self.dt = dt

        self.total_time = self.length*self.dt
        self.init_length = int(np.ceil(self.total_time/self.init_dt))
        for prop, default in prop_defaults.items():
            setattr(self, prop, kwargs.get(prop, default))
            #TODO: check this
            kwargs[prop] = kwargs.get(prop, default)

        self.data_channel1 = data_stream['channel1'].copy()
        self.data_channel2 = data_stream['channel2'].copy()

        freqs = np.fft.rfftfreq(self.length, d=dt)
        df = 1./(self.length*self.dt)

        # whiten data stream
        self.data_stream = data_stream

        self.noise_channel1_inv = 1./np.sqrt(noise_channels['channel1'])*np.sqrt(dt/length)  # dt dt df = dt
        self.noise_channel2_inv = 1./np.sqrt(noise_channels['channel2'])*np.sqrt(dt/length)

        self.data_channel1 *= self.noise_channel1_inv
        self.data_channel2 *= self.noise_channel2_inv

        self.d_d = 4*np.sum([np.abs(self.data_channel1)**2, np.abs(self.data_channel2)**2])

        self.generator = GPUAAK(self.T_fit,
                                self.init_length,
                                self.length,
                                self.init_dt,
                                self.dt,
                                self.LISA,
                                self.backint,
                                self.data_channel1,
                                self.data_channel2,
                                self.noise_channel1_inv,
                                self.noise_channel2_inv)

    def NLL(self, iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                        phi_S, theta_K, phi_K, D, return_snr=False, **kwargs):

        d_h, h_h = self.generator.WaveformThroughLikelihood(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                                                 phi_S, theta_K, phi_K, D)

        if return_snr:
            return np.sqrt(h_h)

        return self.d_d + h_h - 2*d_h
