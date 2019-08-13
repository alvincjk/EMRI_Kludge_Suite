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

        self.noise_channel1_inv = 1./np.sqrt(noise_channels['channel1'])*np.sqrt(df)
        self.noise_channel2_inv = 1./np.sqrt(noise_channels['channel2'])*np.sqrt(df)

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


def test():
    iota = 0.2
    s = 0.8
    p = 10.0
    e = 0.7
    T_fit = 1.0
    init_length = 5000
    length = 500000
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

    tc, hIc, hIIc, timing = cpu_AAK.wave(pars)
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

    like_class = pyGPUAAK(data_stream, noise_channels, length, dt, init_dt, **kwargs)
    num = 1000
    st = time.perf_counter()
    for i in range(num):
        check2 = like_class.NLL(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                                phi_S, theta_K, phi_K, D)
    et = time.perf_counter()
    print('Total time for {} likelihood calculations:'.format(num), (et - st), '\n', 'Time per likelihood calculation:', (et - st)/num, '\n')


if __name__ == "__main__":
    test()
