#include "AAK.h"
#include "Globals.h"
#include <stdlib.h>

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
  double *d_evec;
  double *d_vvec;
  double *d_Mvec;
  double *d_Svec;
  double *d_gimvec;
  double *d_Phivec;
  double *d_alpvec;
  double *d_nuvec;
  double *d_gimdotvec;

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

public:
  /* By using the swig default names INPLACE_ARRAY1, DIM1 in the header
     file (these aren't the names in the implementation file), we're giving
     swig the info it needs to cast to and from numpy arrays.

     If instead the constructor line said
       GPUAdder(int* myarray, int length);

     We would need a line like this in the swig.i file
       %apply (int* ARGOUT_ARRAY1, int DIM1) {(int* myarray, int length)}
   */

  GPUAAK(double T_fit_,
      int length_,
      double dt_,
      bool LISA_,
      bool backint_); // constructor (copies to GPU)

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


  //gets results back from the gpu, putting them in the supplied memory location
  void GetWaveform (double *t_, double* hI_, double* hII_);

};
