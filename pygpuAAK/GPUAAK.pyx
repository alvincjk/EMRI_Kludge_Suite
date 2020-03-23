import numpy as np
cimport numpy as np

assert sizeof(int) == sizeof(np.int32_t)

cdef extern from "AAK_manager.hh":
    cdef cppclass GPUAAKwrap "GPUAAK":
        GPUAAKwrap(double, int, int,
        double, double,
        bool,
        bool)

        void run_phase_trajectory(
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double);

        void input_data(np.complex128_t *hI_f_, np.complex128_t *hII_f_, np.float64_t *channel_ASDinv1_, np.float64_t *channel_ASDinv2_, int len)

        void gpu_gen_AAK(double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double,
                            double)

        void GetWaveform(np.float64_t*, np.float64_t*, np.float64_t*)
        void Likelihood(np.float64_t*)


cdef class GPUAAK:
    cdef int length
    cdef GPUAAKwrap* g

    def __cinit__(self, T_fit, init_length, length, init_dt, dt, LISA, backint):
        self.length = length
        self.g = new GPUAAKwrap(T_fit, init_length, length, init_dt, dt, LISA, backint)

    def input_data(self,
                   np.ndarray[ndim=1, dtype=np.complex128_t] hI_f_,
                   np.ndarray[ndim=1, dtype=np.complex128_t] hII_f_,
                   np.ndarray[ndim=1, dtype=np.float64_t] channel_ASDinv1_,
                   np.ndarray[ndim=1, dtype=np.float64_t] channel_ASDinv2_):

        length_in = len(hI_f_)
        self.g.input_data(&hI_f_[0], &hII_f_[0], &channel_ASDinv1_[0], &channel_ASDinv2_[0], length_in)

    def gpu_gen_AAK(self, iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                        phi_S, theta_K, phi_K, D):

        self.g.gpu_gen_AAK(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)

    def GetWaveform(self, is_Fourier=False):
        cdef np.ndarray[ndim=1, dtype=np.float64_t] t_ = np.zeros(self.length+2, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hI_ = np.zeros(self.length+2, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hII_ = np.zeros(self.length+2, dtype=np.float64)

        self.g.GetWaveform(&t_[0], &hI_[0], &hII_[0])

        return (t_, hI_, hII_)

    def Likelihood(self):
        cdef np.ndarray[ndim=1, dtype=np.float64_t] like_ = np.zeros((2,), dtype=np.float64)
        self.g.Likelihood(&like_[0])
        return like_

    def WaveformThroughLikelihood(self, iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                        phi_S, theta_K, phi_K, D, return_waveform=False):
        self.gpu_gen_AAK(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)
        if return_waveform == True:
            return (0.0, 0.0)
        return self.Likelihood()
