import numpy as np
cimport numpy as np

assert sizeof(int) == sizeof(np.int32_t)

cdef extern from "include/manager.hh":
    cdef cppclass GPUAAKwrap "GPUAAK":
        GPUAAKwrap(double, int,
        double,
        bool,
        bool,
        int)

        void cpu_gen_AAK(double,
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

        void Get_Waveform(np.float64_t*, np.float64_t*, np.float64_t*)
        void gpu_Get_Waveform(np.float64_t*, np.float64_t*, np.float64_t*)

cdef class GPUAAK:
    cdef int length
    cdef GPUAAKwrap* g

    def __cinit__(self, T_fit, length, dt, LISA, backint, to_gpu):
        self.length = length
        self.g = new GPUAAKwrap(T_fit, length, dt, LISA, backint, to_gpu)

    def cpu_gen_AAK(self,iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                        phi_S, theta_K, phi_K, D):

        self.g.cpu_gen_AAK(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)

    def gpu_gen_AAK(self, iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                        phi_S, theta_K, phi_K, D):

        self.g.gpu_gen_AAK(iota, s, p, e, M, mu, gamma, psi, alph, theta_S,
                            phi_S, theta_K, phi_K, D)

    def Get_Waveform(self):
        cdef np.ndarray[ndim=1, dtype=np.float64_t] t_ = np.zeros(self.length, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hI_ = np.zeros(self.length, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hII_ = np.zeros(self.length, dtype=np.float64)

        self.g.Get_Waveform(&t_[0], &hI_[0], &hII_[0])

        return (t_, hI_, hII_)

    def gpu_Get_Waveform(self):
        cdef np.ndarray[ndim=1, dtype=np.float64_t] t_ = np.zeros(self.length, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hI_ = np.zeros(self.length, dtype=np.float64)
        cdef np.ndarray[ndim=1, dtype=np.float64_t] hII_ = np.zeros(self.length, dtype=np.float64)

        self.g.gpu_Get_Waveform(&t_[0], &hI_[0], &hII_[0])
        return (t_, hI_, hII_)
