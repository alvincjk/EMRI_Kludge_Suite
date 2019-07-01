import numpy as np
cimport numpy as np

assert sizeof(int) == sizeof(np.int32_t)

cdef extern from "include/manager.hh":
    cdef cppclass GPUAAKwrap "GPUAAK":
        GPUAAKwrap(double, int,
        double,
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

cdef class GPUAAK:
    cdef int length
    cdef GPUAAKwrap* g

    def __cinit__(self, T_fit, length, dt, LISA, backint):
        self.length = length
        self.g = new GPUAAKwrap(T_fit, length, dt, LISA, backint)

    def gpu_gen_AAK(self, iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                        phi_S, theta_K, phi_K, D):

        self.g.gpu_gen_AAK(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)

    def GetWaveform(self):
        cdef np.ndarray[ndim=1, dtype=np.float64_t] t_ = np.zeros(self.length, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hI_ = np.zeros(self.length, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hII_ = np.zeros(self.length, dtype=np.float64)

        self.g.GetWaveform(&t_[0], &hI_[0], &hII_[0])

        return (t_, hI_, hII_)
