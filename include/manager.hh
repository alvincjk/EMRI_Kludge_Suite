#ifndef __MANAGER_H__
#define __MANAGER_H__

#include "AAK.h"
#include "Globals.h"
#include <stdlib.h>
#include "cuComplex.h"
#include "cublas_v2.h"
#include <cufft.h>
#include <interpolate.hh>

typedef complex<double> cmplx;


class GPUAAK {
  // pointer to the GPU memory where the array is stored
  int length;
  double T_fit;
  double dt;
  bool LISA;
  bool backint;
  int to_gpu;

  size_t double_size;
  double *d_t;
  double *d_hI;
  double *d_hII;
  cmplx *data_channel1;
  cmplx *data_channel2;
  double *noise_channel1_inv;
  double *noise_channel2_inv;
  cuDoubleComplex *d_data_channel1;
  cuDoubleComplex *d_data_channel2;
  double *d_noise_channel1_inv;
  double *d_noise_channel2_inv;
  size_t double_plus_one_size;

  double *tvec;
  double *evec;
  double *vvec;
  double *Mvec;
  double *Svec;
  double *gimvec;
  double *Phivec;
  double *alpvec;
  double *nuvec;
  double *gimdotvec;

  double *d_tvec;
  InterpArrayContainer *d_trajectories;
  InterpArrayContainer *trajectories;
  size_t numBytes;
  InterpArrayContainer d_evec;
  InterpArrayContainer d_vvec;
  InterpArrayContainer d_Mvec;
  InterpArrayContainer d_Svec;
  InterpArrayContainer d_gimvec;
  InterpArrayContainer d_Phivec;
  InterpArrayContainer d_alpvec;
  InterpArrayContainer d_nuvec;
  InterpArrayContainer d_gimdotvec;

  Interpolate interp;

  double iota;
  double s;
  double p;
  double e;
  double M;
  double mu;
  double gamma;
  double psi;
  double alph;
  double theta_S;
  double phi_S;
  double theta_K;
  double phi_K;
  double D;

  int NUM_THREADS;
  int num_blocks;

  double zeta;
  int i_plunge;
  int i_buffer;

  double par[12];

  cufftHandle plan;
  int fft_length;

  cublasHandle_t handle;
  cublasStatus_t stat;

public:
  /* By using the swig default names INPLACE_ARRAY1, DIM1 in the header
     file (these aren't the names in the implementation file), we're giving
     swig the info it needs to cast to and from numpy arrays.

     If instead the constructor line said
       GPUAdder(int* myarray, int length);

     We would need a line like this in the swig.i file
       %apply (int* ARGOUT_ARRAY1, int DIM1) {(int* myarray, int length)}
   */

   GPUAAK (double T_fit_,
       int length_,
       double dt_,
       bool LISA_,
       bool backint_,
       cmplx *data_channel1_,
       cmplx *data_channel2_,
       double *noise_channel1_inv_,
       double *noise_channel2_inv_); // constructor (copies to GPU)

  ~GPUAAK(); // destructor

  void run_phase_trajectory(
      double iota_,
      double s_,
      double p_,
      double e_,
      double M_,
      double mu_,
      double gamma_,
      double psi_,
      double alph_,
      double theta_S_,
      double phi_S_,
      double theta_K_,
      double phi_K_,
      double D_);

  void gpu_gen_AAK(
        double iota_,
        double s_,
        double p_,
        double e_,
        double M_,
        double mu_,
        double gamma_,
        double psi_,
        double alph_,
        double theta_S_,
        double phi_S_,
        double theta_K,
        double phi_K,
        double D_);

    void Likelihood(double *like_out_);

  //gets results back from the gpu, putting them in the supplied memory location
  void GetWaveform (double *t_, double* hI_, double* hII_);

};

static char *_cudaGetErrorEnum(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}

#endif //__MANAGER_H__
