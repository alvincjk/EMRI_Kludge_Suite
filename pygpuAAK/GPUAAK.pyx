import numpy as np
cimport numpy as np

assert sizeof(int) == sizeof(np.int32_t)

cdef extern from "manager.hh":
    cdef cppclass GPUAAKwrap "GPUAAK":
        GPUAAKwrap(double, int, int,
        double, double,
        bool,
        bool,
        np.complex128_t*, np.complex128_t*, np.float64_t*, np.float64_t*)

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

    def __cinit__(self, T_fit, init_length, length, init_dt, dt, LISA, backint,
                  np.ndarray[ndim=1, dtype=np.complex128_t] data_channel1,
                  np.ndarray[ndim=1, dtype=np.complex128_t] data_channel2,
                  np.ndarray[ndim=1, dtype=np.float64_t] noise_channel1_inv,
                  np.ndarray[ndim=1, dtype=np.float64_t] noise_channel2_inv):
        self.length = length
        self.g = new GPUAAKwrap(T_fit, init_length, length, init_dt, dt, LISA, backint, &data_channel1[0], &data_channel2[0], &noise_channel1_inv[0], &noise_channel2_inv[0])

    def gpu_gen_AAK(self, iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                        phi_S, theta_K, phi_K, D):

        self.g.gpu_gen_AAK(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)

    def GetWaveform(self):
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
                        phi_S, theta_K, phi_K, D):
        self.gpu_gen_AAK(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)
        return self.Likelihood()
