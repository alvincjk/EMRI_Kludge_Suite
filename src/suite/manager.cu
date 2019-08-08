/*
This is the central piece of code. This file implements a class
(interface in gpuadder.hh) that takes data in on the cpu side, copies
it to the gpu, and exposes functions (increment and retreive) that let
you perform actions with the GPU

This class will get translated into python via swig
*/

#include <kernel.cu>
#include <manager.hh>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include "cuComplex.h"
#include "cublas_v2.h"
#include <cufft.h>
#include <complex.h>

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAK.h"
#include "interpolate.cu"

using namespace std;

#define BATCH 1

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


GPUAAK::GPUAAK (double T_fit_,
    int init_length_,
    int length_,
    double init_dt_,
    double dt_,
    bool LISA_,
    bool backint_,
    cmplx *data_channel1_,
    cmplx *data_channel2_,
    double *noise_channel1_inv_,
    double *noise_channel2_inv_){

    T_fit = T_fit_;
    init_length = init_length_;
    length = length_;
    init_dt = init_dt_;
    dt = dt_;
    LISA = LISA_;
    backint = backint_;
    data_channel1 = data_channel1_;
    data_channel2 = data_channel2_;
    noise_channel1_inv = noise_channel1_inv_;
    noise_channel2_inv = noise_channel2_inv_;

    to_gpu = 1;

     fft_length = ((int) (length/2)) + 1;

    cudaError_t err;

    // DECLARE ALL THE  NECESSARY STRUCTS

    tvec = new double[init_length+1];
    evec = new double[init_length+1];
    vvec = new double[init_length+1];
    Mvec = new double[init_length+1];
    Svec = new double[init_length+1];
    gimvec = new double[init_length+1];
    Phivec = new double[init_length+1];
    alpvec = new double[init_length+1];
    nuvec = new double[init_length+1];
    gimdotvec = new double[init_length+1];



        size_t numBytes_ = 0;
        trajectories = createInterpArrayContainer(&numBytes_, 9, init_length+1);
        numBytes = numBytes_;
        d_trajectories = createInterpArrayContainer_gpu(numBytes);

        d_evec = trajectories[0];
        d_vvec = trajectories[1];
        d_Mvec = trajectories[2];
        d_Svec = trajectories[3];
        d_gimvec = trajectories[4];
        d_Phivec = trajectories[5];
        d_alpvec = trajectories[6];
        d_nuvec = trajectories[7];
        d_gimdotvec = trajectories[8];

      double_size = length*sizeof(double);
      gpuErrchk(cudaMalloc(&d_t, (length+2)*sizeof(double)));

      gpuErrchk(cudaMalloc(&d_hI, (length+2)*sizeof(double)));
      gpuErrchk(cudaMalloc(&d_hII, (length+2)*sizeof(double)));

      gpuErrchk(cudaMalloc(&d_data_channel1, fft_length*sizeof(cuDoubleComplex)));
      gpuErrchk(cudaMalloc(&d_data_channel2, fft_length*sizeof(cuDoubleComplex)));

      gpuErrchk(cudaMemcpy(d_data_channel1, data_channel1, fft_length*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));
      gpuErrchk(cudaMemcpy(d_data_channel2, data_channel2, fft_length*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice));

      gpuErrchk(cudaMalloc(&d_noise_channel1_inv, fft_length*sizeof(double)));
      gpuErrchk(cudaMalloc(&d_noise_channel2_inv, fft_length*sizeof(double)));

      gpuErrchk(cudaMemcpy(d_noise_channel1_inv, noise_channel1_inv, fft_length*sizeof(double), cudaMemcpyHostToDevice));
      gpuErrchk(cudaMemcpy(d_noise_channel2_inv, noise_channel2_inv, fft_length*sizeof(double), cudaMemcpyHostToDevice));

      double_plus_one_size = (length+1)*sizeof(double);  // TODO reduce size properly
      gpuErrchk(cudaMalloc(&d_tvec, (length+1)*sizeof(double)));


      NUM_THREADS = 256;
      num_blocks = std::ceil((init_length + 1 + NUM_THREADS -1)/NUM_THREADS);
      num_blocks_wave = std::ceil((length + 1 + NUM_THREADS -1)/NUM_THREADS);

     // cufftHandle plan_;
      //plan = plan_;

      //cufftComplex *data;
      if (cufftPlan1d(&plan, length, CUFFT_D2Z, BATCH) != CUFFT_SUCCESS){
        	fprintf(stderr, "CUFFT error: Plan creation failed");
        	return;	}

    stat = cublasCreate(&handle);
  if (stat != CUBLAS_STATUS_SUCCESS) {
          printf ("CUBLAS initialization failed\n");
          exit(0);
      }

    interp.alloc_arrays(init_length + 1, 9);

}

void GPUAAK::gpu_gen_AAK(
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
    double D_){

    GPUAAK::run_phase_trajectory(
        iota_,
        s_,
        p_,
        e_,
        M_,
        mu_,
        gamma_,
        psi_,
        alph_,
        theta_S_,
        phi_S_,
        theta_K_,
        phi_K_,
        D_);

    float milliseconds = 0;


    // Initialize inputs
    // ----- number of modes summed -----
    int nmodes=(int)(30*par[3]);
    if (par[3]<0.135) nmodes=4;
    // ----------

    zeta=par[0]/D/Gpc; // M/D

    cudaError_t err;
    gpuErrchk(cudaMemcpy(d_tvec, tvec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_evec.array, evec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_vvec.array, vvec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_Mvec.array, Mvec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_Svec.array, Svec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_gimvec.array, gimvec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_Phivec.array, Phivec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_alpvec.array, alpvec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_nuvec.array, nuvec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_gimdotvec.array, gimdotvec, (init_length+1)*sizeof(double), cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(d_trajectories, trajectories, numBytes, cudaMemcpyHostToDevice));

    interp.setup(d_trajectories, (init_length + 1), 9);

    /* main: evaluate model at given frequencies */
    kernel_create_waveform<<<num_blocks_wave, NUM_THREADS>>>(d_t, d_hI, d_hII, d_tvec, d_evec, d_vvec, d_Mvec, d_Svec, d_gimvec, d_Phivec, d_alpvec, d_nuvec, d_gimdotvec, iota, theta_S, phi_S, theta_K, phi_K, LISA, init_length, length, nmodes, i_plunge, i_buffer, zeta, M, init_dt, dt);  //iota = lam

     cudaDeviceSynchronize();
     gpuErrchk(cudaGetLastError());


         /*double *hI = new double[length+2];
     cudaMemcpy(hI, d_data_channel1, (length+2)*sizeof(double), cudaMemcpyDeviceToHost);
     for (int i=0; i<200; i+=1){
         //if (i == fft_length-1) hI[2*i + 1] = 0.0;
         printf("%d after: , %e + %e j, %e + %e j\n", i, hI[2*i], hI[2*i + 1], data_channel1[i].real(), data_channel1[i].imag());
     }
     delete[] hI;//*/
}


void GPUAAK::run_phase_trajectory(
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
    double D_){

    iota = iota_;
    s = s_;
    p = p_;
    e = e_;
    M = M_;
    mu = mu_;
    gamma = gamma_;
    psi = psi_;
    alph = alph_;
    theta_S = theta_S_;
    phi_S = phi_S_;
    theta_K = theta_K_;
    phi_K = phi_K_;
    D = D_;

    clock_t ticks=clock();

    GKTrajFast gktraj3(cos(iota),s);
    gktraj3.p=p;
    gktraj3.ecc=e;
    int maxsteps=100;
    int steps=0;
    double dt_fit=min(T_fit,init_length*init_dt/SOLARMASSINSEC/M/M*mu)/(maxsteps-1);
    TrajData *traj3;
    traj3=(TrajData*)malloc((size_t)((maxsteps+1)*sizeof(TrajData)));
    gktraj3.Eccentric(dt_fit,traj3,maxsteps,steps);
    double Omega_t[3],ang[3],map_t[3],e_traj[steps],v_map[steps],M_map[steps],s_map[steps],dt_map[steps];
    double Phi;
    for(int i=1;i<=steps;i++){
      IEKG geodesic_t(traj3[i].p,traj3[i].ecc,traj3[i].cosiota,s);
      geodesic_t.Frequencies(Omega_t);
      if(i==1){
        ParAng(ang,e,iota,gamma,psi,theta_S,phi_S,theta_K,phi_K,alph,geodesic_t.zedminus);
        Phi=ang[0]; // initial mean anomaly
      }
      ParMap(map_t,Omega_t,traj3[i].p,M,s,traj3[i].ecc,iota);
      e_traj[i-1]=traj3[i].ecc;
      v_map[i-1]=map_t[0]; // mapped initial velocity in c
      M_map[i-1]=map_t[1]; // mapped BH mass in solar masses
      s_map[i-1]=map_t[2]; // mapped spin parameter a/M = S/M^2
      dt_map[i-1]=traj3[i].t*SOLARMASSINSEC*M*M/mu;
    }

    //GenWave(t,hI,hII,AAK.dt,AAK.length,e_traj,v_map,AAK.M,M_map,AAK.mu,AAK.s,s_map,AAK.D,AAK.iota,AAK.gamma,Phi,AAK.theta_S,AAK.phi_S,AAK.alpha,AAK.theta_K,AAK.phi_K,dt_map,steps,AAK.backint,AAK.LISA,false);

    par[0]=mu*SOLARMASSINSEC;
    par[1]=M_map[0]*SOLARMASSINSEC;
    par[2]=s_map[0];
    par[3]=e_traj[0];
    par[4]=iota;  // TODO: check this
    par[5]=gamma;
    par[6]=Phi;
    par[7]=theta_S;
    par[8]=phi_S;
    par[9]=theta_K;
    par[10]=phi_K;
    par[11]=alph;


    PNevolution(tvec,evec,vvec,Mvec,Svec,gimvec,Phivec,alpvec,nuvec,gimdotvec,init_dt,init_length,par,e_traj,v_map,M,M_map,s,s_map,dt_map,steps,&i_plunge,&i_buffer,backint);

}

void GPUAAK::Likelihood (double *like_out_){

    //cudaMemcpy(hI, d_hI, (length+2)*sizeof(double), cudaMemcpyDeviceToHost);

    if (cufftExecD2Z(plan, d_hI, (cufftDoubleComplex*)d_hI) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: ExecC2C Forward failed");
    return;}
    cudaDeviceSynchronize();
    gpuErrchk(cudaGetLastError());

    if (cufftExecD2Z(plan, d_hII, (cufftDoubleComplex*)d_hII) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: ExecC2C Forward failed");
    return;}

    cudaDeviceSynchronize();
    gpuErrchk(cudaGetLastError());

    likelihood_prep<<<num_blocks, NUM_THREADS>>>((cuDoubleComplex*)d_hI, (cuDoubleComplex*)d_hII, d_noise_channel1_inv, d_noise_channel2_inv, fft_length);
    cudaDeviceSynchronize();
    gpuErrchk(cudaGetLastError());


    //printf("checkcheckcheck\n");

     double d_h = 0.0;
     double h_h = 0.0;
     char * status;
     double res;
     cuDoubleComplex result;

         stat = cublasZdotc(handle, fft_length,
                 (cuDoubleComplex*)d_hI, 1,
                 (cuDoubleComplex*)d_data_channel1, 1,
                 &result);
         status = _cudaGetErrorEnum(stat);
          cudaDeviceSynchronize();

          if (stat != CUBLAS_STATUS_SUCCESS) {
                  exit(0);
              }
         d_h += cuCreal(result);
         //printf("channel1 d_h: %e\n", cuCreal(result));

         stat = cublasZdotc(handle, fft_length,
                 (cuDoubleComplex*)d_hII, 1,
                 (cuDoubleComplex*)d_data_channel2, 1,
                 &result);
         status = _cudaGetErrorEnum(stat);
          cudaDeviceSynchronize();

          if (stat != CUBLAS_STATUS_SUCCESS) {
                  exit(0);
              }
         d_h += cuCreal(result);
         //printf("channel2 d_h: %e\n", cuCreal(result));

        stat = cublasZdotc(handle, fft_length,
                     (cuDoubleComplex*)d_hI, 1,
                     (cuDoubleComplex*)d_hI, 1,
                     &result);
             status = _cudaGetErrorEnum(stat);
              cudaDeviceSynchronize();

              if (stat != CUBLAS_STATUS_SUCCESS) {
                      exit(0);
                  }
             h_h += cuCreal(result);
             //printf("channel1 h_h: %e\n", cuCreal(result));

             stat = cublasZdotc(handle, fft_length,
                     (cuDoubleComplex*)d_hII, 1,
                     (cuDoubleComplex*)d_hII, 1,
                     &result);
             status = _cudaGetErrorEnum(stat);
              cudaDeviceSynchronize();

              if (stat != CUBLAS_STATUS_SUCCESS) {
                      exit(0);
                  }
             h_h += cuCreal(result);
             //printf("channel2 h_h: %e\n", cuCreal(result));

    //printf("dh: %e, hh: %e\n", d_h, h_h);
     like_out_[0] = 4*d_h;
     like_out_[1] = 4*h_h;

}

void GPUAAK::GetWaveform (double *t_, double* hI_, double* hII_) {
 gpuErrchk(cudaMemcpy(t_, d_t, (length+2)*sizeof(double), cudaMemcpyDeviceToHost));
 gpuErrchk(cudaMemcpy(hI_, d_hI, (length+2)*sizeof(double), cudaMemcpyDeviceToHost));
 gpuErrchk(cudaMemcpy(hII_, d_hII, (length+2)*sizeof(double), cudaMemcpyDeviceToHost));
}//*/

GPUAAK::~GPUAAK() {
  delete[] tvec;
  delete[] evec;
  delete[] vvec;
  delete[] Mvec;
  delete[] Svec;
  delete[] gimvec;
  delete[] Phivec;
  delete[] alpvec;
  delete[] nuvec;
  delete[] gimdotvec;

  cudaFree(d_t);
  cudaFree(d_hI);
  cudaFree(d_hII);
  cudaFree(d_tvec);
  destroyInterpArrayContainer(d_trajectories, trajectories, 9);
  cudaFree(d_data_channel1);
  cudaFree(d_data_channel2);
  cudaFree(d_noise_channel1_inv);
  cudaFree(d_noise_channel2_inv);

  cufftDestroy(plan);

}
