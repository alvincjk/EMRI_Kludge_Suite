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
#include "WaveformContainer.h"

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAK.h"

using namespace std;

GPUAAK::GPUAAK (double T_fit_,
    int length_,
    double dt_,
    bool LISA_,
    bool backint_){

    T_fit = T_fit_;
    length = length_;
    dt = dt_;
    LISA = LISA_;
    backint = backint_;

    to_gpu = 1;

    cudaError_t err;

    // DECLARE ALL THE  NECESSARY STRUCTS
    transfer_info = new WaveformContainer;

    tvec = new double[length+1];
    evec = new double[length+1];
    vvec = new double[length+1];
    Mvec = new double[length+1];
    Svec = new double[length+1];
    gimvec = new double[length+1];
    Phivec = new double[length+1];
    alpvec = new double[length+1];
    nuvec = new double[length+1];
    gimdotvec = new double[length+1];

      double_size = length*sizeof(double);
      err = cudaMalloc(&d_t, length*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_hI, length*sizeof(double));
      assert(err == 0);
      err = cudaMalloc(&d_hII, length*sizeof(double));
      assert(err == 0);

      double_plus_one_size = (length+1)*sizeof(double);  // TODO reduce size properly
      err = cudaMalloc(&d_tvec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_evec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_vvec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_Mvec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_Svec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_gimvec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_Phivec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_alpvec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_nuvec, (length+1)*sizeof(double));
      assert(err == 0);

      err = cudaMalloc(&d_gimdotvec, (length+1)*sizeof(double));
      assert(err == 0);

      NUM_THREADS = 256;
      num_blocks = std::ceil((length + NUM_THREADS -1)/NUM_THREADS);
      printf("blocks %d\n", num_blocks);

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


    // Initialize inputs
    // ----- number of modes summed -----
    int nmodes=(int)(30*par[3]);
    if (par[3]<0.135) nmodes=4;
    // ----------

    zeta=par[0]/dist/Gpc; // M/D

    cudaError_t err;
    err = cudaMemcpy(d_tvec, tvec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_evec, evec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_vvec, vvec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_Mvec, Mvec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_Svec, Svec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_gimvec, gimvec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_Phivec, Phivec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_alpvec, alpvec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_nuvec, nuvec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    err = cudaMemcpy(d_gimdotvec, gimdotvec, (length+1)*sizeof(double), cudaMemcpyHostToDevice);
    assert(err == 0);

    /* main: evaluate model at given frequencies */

    kernel_create_waveform<<<num_blocks, NUM_THREADS>>>(d_t, d_hI, d_hII, d_tvec, d_evec, d_vvec, d_Mvec, d_Svec, d_gimvec, d_Phivec, d_alpvec, d_nuvec, d_gimdotvec, iota, theta_S, phi_S, theta_K, phi_K, LISA, length, nmodes, i_plunge, i_buffer, zeta, M, dt);  //iota = lam

     cudaDeviceSynchronize();
     err = cudaGetLastError();
     assert(err == 0);
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
    double dt_fit=min(T_fit,length*dt/SOLARMASSINSEC/M/M*mu)/(maxsteps-1);
    TrajData *traj3;
    traj3=(TrajData*)malloc((size_t)((maxsteps+1)*sizeof(TrajData)));
    gktraj3.Eccentric(dt_fit,traj3,maxsteps,steps);
    double Omega_t[3],ang[3],map_t[3],e_traj[steps],v_map[steps],M_map[steps],s_map[steps],dt_map[steps];
    double Phi;
    for(int i=1;i<=steps;i++){
      IEKG geodesic_t(traj3[i].p,traj3[i].ecc,traj3[i].cosiota,s);
      geodesic_t.Frequencies(Omega_t);
      if(i==1){
        ParAng(ang,e,iota,gamma,psi,theta_S,phi_S,theta_K,phi_K,alpha,geodesic_t.zedminus);
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

    double par[12];
    par[0]=mu*SOLARMASSINSEC;
    par[1]=M_map[0]*SOLARMASSINSEC;
    par[2]=S_map[0];
    par[3]=e_traj[0];
    par[4]=inc;
    par[5]=gim0;
    par[6]=Phi0;
    par[7]=qS;
    par[8]=phiS;
    par[9]=qK;
    par[10]=phiK;
    par[11]=alph;


    PNevolution(tvec,evec,vvec,Mvec,Svec,gimvec,Phivec,alpvec,nuvec,gimdotvec,timestep,vlength,par,e_traj,v_map,M,M_map,s,S_map,dt_map,steps,&i_plunge,&i_buffer,backint);

}

void GPUAAK::GetWaveform (double *t_, double* hI_, double* hII_) {
assert ((to_gpu == 0) || (to_gpu == 2));
 memcpy(t_, t, length*sizeof(double));
 memcpy(hI_, hI, length*sizeof(double));
 memcpy(hII_, hII, length*sizeof(double));
}

GPUAAK::~GPUAAK() {
  delete[] evec;
  delete[] vvec;
  delete[] Mvec;
  delete[] Svec;
  delete[] gimvec;
  delete[] Phivec;
  delete[] alpvec;
  delete[] nuvec;
  delete[] gimdotvec;
  delete[] transfer_info;

  cudaFree(d_t);
  cudaFree(d_hI);
  cudaFree(d_hII);
  cudaFree(d_evec);
  cudaFree(d_vvec);
  cudaFree(d_Mvec);
  cudaFree(d_Svec);
  cudaFree(d_gimvec);
  cudaFree(d_Phivec);
  cudaFree(d_alpvec);
  cudaFree(d_nuvec);
  cudaFree(d_gimdotvec);

}
