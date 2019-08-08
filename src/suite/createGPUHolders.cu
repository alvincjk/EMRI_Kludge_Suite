#include <assert.h>
#include <iostream>
#include "interpolate.hh"
#include <complex>
#include "cuComplex.h"


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


TrajectoryContainer * gpu_create_container(int max_length){
        TrajectoryContainer * trajectories;

        double *tvec, *evec, *Mvec, *Svec, *gimvec, *Phivec, *alpvec, *nuvec, *gimdotvec;
        double *evec_coeff_1, *Mvec_coeff_1, *Svec_coeff_1, *gimvec_coeff_1, *Phivec_coeff_1, *alpvec_coeff_1, *nuvec_coeff_1, *gimdotvec_coeff_1;
        double *evec_coeff_2, *Mvec_coeff_2, *Svec_coeff_2, *gimvec_coeff_2, *Phivec_coeff_2, *alpvec_coeff_2, *nuvec_coeff_2, *gimdotvec_coeff_2;
        double *evec_coeff_3, *Mvec_coeff_3, *Svec_coeff_3, *gimvec_coeff_3, *Phivec_coeff_3, *alpvec_coeff_3, *nuvec_coeff_3, *gimdotvec_coeff_3;

        gpuErrchk(cudaMalloc(&trajectories, sizeof(TrajectoryContainer)));

            // waveform
            gpuErrchk(cudaMalloc(&tvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&evec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec, max_length*sizeof(double)));

            gpuErrchk(cudaMalloc(&evec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec_coeff_1, (max_length-1)*sizeof(double)));

            gpuErrchk(cudaMalloc(&evec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec_coeff_2, (max_length-1)*sizeof(double)));

            gpuErrchk(cudaMalloc(&evec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec_coeff_3, (max_length-1)*sizeof(double)));


            gpuErrchk(cudaMemcpy(&(trajectories.tvec), &(tvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.evec), &(evec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec), &(Mvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec), &(Svec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec), &(gimvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec), &(Phivec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec), &(alpvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec), &(nuvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec), &(gimdotvec), sizeof(double *), cudaMemcpyHostToDevice));

            gpuErrchk(cudaMemcpy(&(trajectories.evec_coeff_1), &(evec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec_coeff_1), &(Mvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec_coeff_1), &(Svec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec_coeff_1), &(gimvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec_coeff_1), &(Phivec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec_coeff_1), &(alpvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec_coeff_1), &(nuvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec_coeff_1), &(gimdotvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));

            gpuErrchk(cudaMemcpy(&(trajectories.evec_coeff_2), &(evec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec_coeff_2), &(Mvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec_coeff_2), &(Svec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec_coeff_2), &(gimvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec_coeff_2), &(Phivec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec_coeff_2), &(alpvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec_coeff_2), &(nuvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec_coeff_2), &(gimdotvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));

            gpuErrchk(cudaMemcpy(&(trajectories.evec_coeff_3), &(evec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec_coeff_3), &(Mvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec_coeff_3), &(Svec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec_coeff_3), &(gimvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec_coeff_3), &(Phivec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec_coeff_3), &(alpvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec_coeff_3), &(nuvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec_coeff_3), &(gimdotvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));


        return trajectories;
}

void gpu_destroy_container(TrajectoryContainer * trajectories){
        gpuErrchk(cudaFree(trajectories.tvec));
        gpuErrchk(cudaFree(trajectories.evec));
        gpuErrchk(cudaFree(trajectories.Mvec));
        gpuErrchk(cudaFree(trajectories.Svec));
        gpuErrchk(cudaFree(trajectories.gimvec));
        gpuErrchk(cudaFree(trajectories.Phivec));
        gpuErrchk(cudaFree(trajectories.alpvec));
        gpuErrchk(cudaFree(trajectories.nuvec));
        gpuErrchk(cudaFree(trajectories.gimdotvec));

        gpuErrchk(cudaFree(trajectories.evec_coeff_1));
        gpuErrchk(cudaFree(trajectories.Mvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.Svec_coeff_1));
        gpuErrchk(cudaFree(trajectories.gimvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.Phivec_coeff_1));
        gpuErrchk(cudaFree(trajectories.alpvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.nuvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.gimdotvec_coeff_1));

        gpuErrchk(cudaFree(trajectories.evec_coeff_2));
        gpuErrchk(cudaFree(trajectories.Mvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.Svec_coeff_2));
        gpuErrchk(cudaFree(trajectories.gimvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.Phivec_coeff_2));
        gpuErrchk(cudaFree(trajectories.alpvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.nuvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.gimdotvec_coeff_2));

        gpuErrchk(cudaFree(trajectories.evec_coeff_3));
        gpuErrchk(cudaFree(trajectories.Mvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.Svec_coeff_3));
        gpuErrchk(cudaFree(trajectories.gimvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.Phivec_coeff_3));
        gpuErrchk(cudaFree(trajectories.alpvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.nuvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.gimdotvec_coeff_3));


    gpuErrchk(cudaFree(trajectories));
}



TrajectoryContainer * cpu_create_container(int max_length){
        TrajectoryContainer * trajectories;

        double *tvec, *evec, *Mvec, *Svec, *gimvec, *Phivec, *alpvec, *nuvec, *gimdotvec;
        double *evec_coeff_1, *Mvec_coeff_1, *Svec_coeff_1, *gimvec_coeff_1, *Phivec_coeff_1, *alpvec_coeff_1, *nuvec_coeff_1, *gimdotvec_coeff_1;
        double *evec_coeff_2, *Mvec_coeff_2, *Svec_coeff_2, *gimvec_coeff_2, *Phivec_coeff_2, *alpvec_coeff_2, *nuvec_coeff_2, *gimdotvec_coeff_2;
        double *evec_coeff_3, *Mvec_coeff_3, *Svec_coeff_3, *gimvec_coeff_3, *Phivec_coeff_3, *alpvec_coeff_3, *nuvec_coeff_3, *gimdotvec_coeff_3;

        gpuErrchk(cudaMalloc(&trajectories, sizeof(TrajectoryContainer)));

            // waveform
            gpuErrchk(cudaMalloc(&tvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&evec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec, max_length*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec, max_length*sizeof(double)));

            gpuErrchk(cudaMalloc(&evec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec_coeff_1, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec_coeff_1, (max_length-1)*sizeof(double)));

            gpuErrchk(cudaMalloc(&evec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec_coeff_2, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec_coeff_2, (max_length-1)*sizeof(double)));

            gpuErrchk(cudaMalloc(&evec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Mvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Svec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&Phivec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&alpvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&nuvec_coeff_3, (max_length-1)*sizeof(double)));
            gpuErrchk(cudaMalloc(&gimdotvec_coeff_3, (max_length-1)*sizeof(double)));


            gpuErrchk(cudaMemcpy(&(trajectories.tvec), &(tvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.evec), &(evec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec), &(Mvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec), &(Svec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec), &(gimvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec), &(Phivec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec), &(alpvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec), &(nuvec), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec), &(gimdotvec), sizeof(double *), cudaMemcpyHostToDevice));

            gpuErrchk(cudaMemcpy(&(trajectories.evec_coeff_1), &(evec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec_coeff_1), &(Mvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec_coeff_1), &(Svec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec_coeff_1), &(gimvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec_coeff_1), &(Phivec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec_coeff_1), &(alpvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec_coeff_1), &(nuvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec_coeff_1), &(gimdotvec_coeff_1), sizeof(double *), cudaMemcpyHostToDevice));

            gpuErrchk(cudaMemcpy(&(trajectories.evec_coeff_2), &(evec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec_coeff_2), &(Mvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec_coeff_2), &(Svec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec_coeff_2), &(gimvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec_coeff_2), &(Phivec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec_coeff_2), &(alpvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec_coeff_2), &(nuvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec_coeff_2), &(gimdotvec_coeff_2), sizeof(double *), cudaMemcpyHostToDevice));

            gpuErrchk(cudaMemcpy(&(trajectories.evec_coeff_3), &(evec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Mvec_coeff_3), &(Mvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Svec_coeff_3), &(Svec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimvec_coeff_3), &(gimvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.Phivec_coeff_3), &(Phivec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.alpvec_coeff_3), &(alpvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.nuvec_coeff_3), &(nuvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(&(trajectories.gimdotvec_coeff_3), &(gimdotvec_coeff_3), sizeof(double *), cudaMemcpyHostToDevice));


        return trajectories;
}

void cpu_destroy_container(TrajectoryContainer * trajectories){
        gpuErrchk(cudaFree(trajectories.tvec));
        gpuErrchk(cudaFree(trajectories.evec));
        gpuErrchk(cudaFree(trajectories.Mvec));
        gpuErrchk(cudaFree(trajectories.Svec));
        gpuErrchk(cudaFree(trajectories.gimvec));
        gpuErrchk(cudaFree(trajectories.Phivec));
        gpuErrchk(cudaFree(trajectories.alpvec));
        gpuErrchk(cudaFree(trajectories.nuvec));
        gpuErrchk(cudaFree(trajectories.gimdotvec));

        gpuErrchk(cudaFree(trajectories.evec_coeff_1));
        gpuErrchk(cudaFree(trajectories.Mvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.Svec_coeff_1));
        gpuErrchk(cudaFree(trajectories.gimvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.Phivec_coeff_1));
        gpuErrchk(cudaFree(trajectories.alpvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.nuvec_coeff_1));
        gpuErrchk(cudaFree(trajectories.gimdotvec_coeff_1));

        gpuErrchk(cudaFree(trajectories.evec_coeff_2));
        gpuErrchk(cudaFree(trajectories.Mvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.Svec_coeff_2));
        gpuErrchk(cudaFree(trajectories.gimvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.Phivec_coeff_2));
        gpuErrchk(cudaFree(trajectories.alpvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.nuvec_coeff_2));
        gpuErrchk(cudaFree(trajectories.gimdotvec_coeff_2));

        gpuErrchk(cudaFree(trajectories.evec_coeff_3));
        gpuErrchk(cudaFree(trajectories.Mvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.Svec_coeff_3));
        gpuErrchk(cudaFree(trajectories.gimvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.Phivec_coeff_3));
        gpuErrchk(cudaFree(trajectories.alpvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.nuvec_coeff_3));
        gpuErrchk(cudaFree(trajectories.gimdotvec_coeff_3));


    gpuErrchk(cudaFree(trajectories));
}
