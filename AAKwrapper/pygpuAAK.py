import numpy as np
from scipy import constants as ct
from .utils.convert import Converter, Recycler

import tdi

from gpuAAK import GPUAAK

import time

MTSUN = 1.989e30*ct.G/ct.c**3


class pygpuAAK:
    def __init__(self, data_stream, noise_channels, dt, **kwargs):
        """
        data_stream (dict): keys X, Y, Z or A, E, T
        """
        prop_defaults = {
            'T_fit': 1.0,
            'LISA': True,
            'backint' True
        }

        self.dt = dt
        self.T_fit = T_fit
        self.length = len(data_stream['channel1'])
        for prop, default in prop_defaults.items():
            setattr(self, prop, kwargs.get(prop, default))
            #TODO: check this
            kwargs[prop] = kwargs.get(prop, default)

        self.data_channel1 = np.fft.rfft(data_stream['channel1'])
        self.data_channel2 = np.fft.rfft(data_stream['channel2'])

        freqs = np.fft.rfftfreq(len(data_stream['channel1']), d=dt)
        df = freqs[1] - freq[0]

        # whiten data stream
        self.data_stream = data_stream
        self.data_channel1 *= 1./np.sqrt(noise_channels['channel1'])*np.sqrt(df)
        self.data_channel2 *= 1./np.sqrt(noise_channels['channel2'])*np.sqrt(df)

        self.noise_channel1_inv = 1./np.sqrt(noise_channels['channel1'])*np.sqrt(df)
        self.noise_channel2_inv = 1./np.sqrt(noise_channels['channel2'])*np.sqrt(df)

        self.d_d = 4*np.sum([np.abs(self.data_channel1)**2, np.abs(self.data_channel2)**2])

        self.generator = GPUAAK(self.T_fit,
                                self.length,
                                self.dt,
                                self.LISA,
                                self.backint,
                                self.data_channel1,
                                self.data_channel2,
                                self.noise_channel1_inv,
                                self.noise_channel2_inv)

    def NLL(self, m1, m2, a1, a2, distance,
                 phiRef, inc, lam, beta,
                 psi, tRef_wave_frame, tRef_sampling_frame, freqs=None, return_amp_phase=False, return_TDI=False):

        Msec = (m1+m2)*MTSUN
        # merger frequency for 22 mode amplitude in phenomD
        merger_freq = 0.018/Msec
        fRef = 0.0

        if freqs is None:
            upper_freq = self.max_dimensionless_freq/Msec
            lower_freq = self.min_dimensionless_freq/Msec
            freqs = np.logspace(np.log10(lower_freq), np.log10(upper_freq), self.max_length_init)

        out = self.generator.WaveformThroughLikelihood(freqs,
                                              m1, m2,  # solar masses
                                              a1, a2,
                                              distance, phiRef, fRef,
                                              inc, lam, beta, psi,
                                              self.t0, tRef_wave_frame, tRef_sampling_frame, merger_freq,
                                              return_amp_phase=return_amp_phase,
                                              return_TDI=return_TDI)

        if return_amp_phase or return_TDI:
            return out

        d_h, h_h = out

        return self.d_d + h_h - 2*d_h


def create_data_set(l_vals,  m_vals, t0, waveform_params, data_freqs=None, TDItag='AET', num_data_points=int(2**19), num_generate_points=int(2**18), df=None, min_dimensionless_freq=1e-4, max_dimensionless_freq=1.0, add_noise=None, **kwargs):
    key_list = list(waveform_params.keys())
    converter = Converter(key_list, **kwargs)
    recycler = Recycler(key_list, **kwargs)

    vals = np.array([waveform_params[key] for key in key_list])

    tRef_sampling_frame = np.exp(vals[10])

    vals = recycler.recycle(vals)
    vals = converter.convert(vals)

    waveform_params = {key: vals[i] for i, key in enumerate(key_list)}

    waveform_params['tRef_sampling_frame'] = tRef_sampling_frame

    if 'ln_m1' in waveform_params:
        waveform_params['m1'] = waveform_params['ln_m1']
        waveform_params['m2'] = waveform_params['ln_m2']
    if 'ln_mT' in waveform_params:
        # has been converted
        waveform_params['m1'] = waveform_params['ln_mT']
        waveform_params['m2'] = waveform_params['mr']

    if 'chi_s' in waveform_params:
        waveform_params['a1'] = waveform_params['chi_s']
        waveform_params['a2'] = waveform_params['chi_a']

    waveform_params['distance'] = waveform_params['ln_distance']
    waveform_params['tRef_wave_frame'] = waveform_params['ln_tRef']
    waveform_params['fRef'] = 0.0

    m1 = waveform_params['m1']
    m2 = waveform_params['m2']
    Msec = (m1+m2)*MTSUN
    merger_freq = 0.018/Msec

    if data_freqs is None:
        if add_noise is not None:
            fs = add_noise['fs']
            t_obs_dur = waveform_params['t_obs_dur']
            num_data_points = int(t_obs_dur*ct.Julian_year*fs)
            noise_freqs = np.fft.rfftfreq(num_data_points, 1/fs)
            data_freqs = noise_freqs[noise_freqs >= add_noise['min_freq']]

        else:
            m1 = waveform_params['m1']
            m2 = waveform_params['m2']
            Msec = (m1+m2)*MTSUN
            upper_freq = max_dimensionless_freq/Msec
            lower_freq = min_dimensionless_freq/Msec
            merger_freq = 0.018/Msec
            if df is None:
                data_freqs = np.logspace(np.log10(lower_freq), np.log10(upper_freq), num_data_points)
            else:
                data_freqs = np.arange(fmin, fmax+df, df)

    if add_noise is not None:
        if TDItag == 'AET':
            # assumes gaussian noise
            noise_channel1 = np.abs(np.random.normal(0.0, np.sqrt(1./2.*tdi.noisepsd_AE(data_freqs, **kwargs['noise_kwargs']))))*np.exp(1j*np.random.uniform(0, 2*np.pi, size=data_freqs.shape))
            noise_channel2 = np.abs(np.random.normal(0.0, np.sqrt(1./2.*tdi.noisepsd_AE(data_freqs, **kwargs['noise_kwargs']))))*np.exp(1j*np.random.uniform(0, 2*np.pi, size=data_freqs.shape))
            noise_channel3 = np.abs(np.random.normal(0.0, np.sqrt(1./2.*tdi.noisepsd_T(data_freqs, **kwargs['noise_kwargs']))))*np.exp(1j*np.random.uniform(0, 2*np.pi, size=data_freqs.shape))

        else:
            # assumes gaussian noise
            noise_channel1 = np.abs(np.random.normal(0.0, np.sqrt(1./2.*tdi.noisepsd_XYZ(data_freqs, **kwargs['noise_kwargs']))))*np.exp(1j*np.random.uniform(0, 2*np.pi, size=data_freqs.shape))
            noise_channel2 = nnp.abs(np.random.normal(0.0, p.sqrt(1./2.*tdi.noisepsd_XYZ(data_freqs, **kwargs['noise_kwargs']))))*np.exp(1j*np.random.uniform(0, 2*np.pi, size=data_freqs.shape))
            noise_channel3 = np.abs(np.random.normal(0.0, np.sqrt(1./2.*tdi.noisepsd_XYZ(data_freqs, **kwargs['noise_kwargs']))))*np.exp(1j*np.random.uniform(0, 2*np.pi, size=data_freqs.shape))

    generate_freqs = np.logspace(np.log10(data_freqs.min()), np.log10(data_freqs.max()), num_generate_points)

    fake_data = np.zeros_like(data_freqs, dtype=np.complex128)
    fake_ASD = np.ones_like(data_freqs)

    if TDItag == 'AET':
        TDItag_in = 2

    elif TDItag == 'XYZ':
        TDItag_in = 1

    phenomHM = PhenomHM(len(generate_freqs), l_vals, m_vals, data_freqs, fake_data, fake_data, fake_data, fake_ASD, fake_ASD, fake_ASD, TDItag_in, waveform_params['t_obs_dur'])

    phenomHM.gen_amp_phase(generate_freqs, waveform_params['m1'],  # solar masses
                 waveform_params['m2'],  # solar masses
                 waveform_params['a1'],
                 waveform_params['a2'],
                 waveform_params['distance'],
                 waveform_params['phiRef'],
                 waveform_params['fRef'])

    phenomHM.LISAresponseFD(waveform_params['inc'], waveform_params['lam'], waveform_params['beta'], waveform_params['psi'], waveform_params['t0'], waveform_params['tRef_wave_frame'], waveform_params['tRef_sampling_frame'], merger_freq)
    phenomHM.setup_interp_wave()
    phenomHM.setup_interp_response()
    phenomHM.perform_interp()

    channel1, channel2, channel3 = phenomHM.GetTDI()

    if channel1.ndim > 1:
        channel1, channel2, channel3 = channel1.sum(axis=0), channel2.sum(axis=0), channel3.sum(axis=0)

    if add_noise is not None:
        channel1 = channel1 + noise_channel1
        channel2 = channel2 + noise_channel2
        channel3 = channel3 + noise_channel3

    data_stream = {TDItag[0]: channel1, TDItag[1]: channel2, TDItag[2]: channel3}
    return data_freqs, data_stream
