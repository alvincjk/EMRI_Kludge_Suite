#include "stdio.h"
#include <assert.h>
#include <cusparse_v2.h>
#include "interpolate.hh"

#define gpuErrchk_here(ans) { gpuAssert_here((ans), __FILE__, __LINE__); }
inline void gpuAssert_here(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define ERR_NE(X,Y) do { if ((X) != (Y)) { \
                             fprintf(stderr,"Error in %s at %s:%d\n",__func__,__FILE__,__LINE__); \
                             exit(-1);}} while(0)

#define CUDA_CALL(X) ERR_NE((X),cudaSuccess)
#define CUSPARSE_CALL(X) ERR_NE((X),CUSPARSE_STATUS_SUCCESS)
using namespace std;


InterpArrayContainer * createInterpArrayContainer(size_t *numBytes, int num_arr, int num_points){

    InterpArrayContainer *cpu_array_container;
    size_t InterpArrayContainer_size = sizeof(InterpArrayContainer);
    *numBytes = num_arr*InterpArrayContainer_size;

    cpu_array_container = (InterpArrayContainer*)malloc(*numBytes);

    for (int i=0; i<num_arr; i++){

        gpuErrchk_here(cudaMalloc( (void**)&(cpu_array_container[i].array), num_points*sizeof(double) ));
        gpuErrchk_here(cudaMalloc( (void**)&(cpu_array_container[i].coeff_1), (num_points-1)*sizeof(double) ));
        gpuErrchk_here(cudaMalloc( (void**)&(cpu_array_container[i].coeff_2), (num_points-1)*sizeof(double) ));
        gpuErrchk_here(cudaMalloc( (void**)&(cpu_array_container[i].coeff_3), (num_points-1)*sizeof(double) ));
    }
    return cpu_array_container;

    //cudaMalloc((void**)&gpu_array_container, *numBytes);
}

InterpArrayContainer * createInterpArrayContainer_gpu(size_t numBytes, InterpArrayContainer *cpu_array_container){
    InterpArrayContainer *gpu_array_container;
    cudaMalloc((void**)&gpu_array_container, numBytes);
    cudaMemcpy(gpu_array_container, cpu_array_container, numBytes, cudaMemcpyHostToDevice);
    return gpu_array_container;
}

void destroyInterpArrayContainer(InterpArrayContainer * gpu_array_container, InterpArrayContainer *cpu_array_container, int num_arr){
    for (int i=0; i<num_arr; i++){
        gpuErrchk_here(cudaFree(cpu_array_container[i].array));
        gpuErrchk_here(cudaFree(cpu_array_container[i].coeff_1));
        gpuErrchk_here(cudaFree(cpu_array_container[i].coeff_2));
        gpuErrchk_here(cudaFree(cpu_array_container[i].coeff_3));
    }

    gpuErrchk_here(cudaFree(gpu_array_container));
    free(cpu_array_container);
}



Interpolate::Interpolate(){
    int pass = 0;
}

__host__
void Interpolate::alloc_arrays(int max_length_init, int num_arr){
    gpuErrchk_here(cudaMalloc(&d_B, max_length_init*num_arr*sizeof(double)));
    gpuErrchk_here(cudaMalloc(&d_dl, max_length_init*num_arr*sizeof(double)));
    gpuErrchk_here(cudaMalloc(&d_d, max_length_init*num_arr*sizeof(double)));
    gpuErrchk_here(cudaMalloc(&d_du, max_length_init*num_arr*sizeof(double)));
}

__device__ void prep_splines(int i, int length, double *b, double *ud, double *diag, double *ld, double *x, double *y){
  double dx1, dx2, d, slope1, slope2;
  if (i == length - 1){
    dx1 = x[length - 2] - x[length - 3];
    dx2 = x[length - 1] - x[length - 2];
    d = x[length - 1] - x[length - 3];

    slope1 = (y[length - 2] - y[length - 3])/dx1;
    slope2 = (y[length - 1] - y[length - 2])/dx2;

    b[length - 1] = ((dx2*dx2*slope1 +
                             (2*d + dx2)*dx1*slope2) / d);
    diag[length - 1] = dx1;
    ld[length - 1] = d;
    ud[length - 1] = 0.0;

  } else if (i == 0){
      dx1 = x[1] - x[0];
      dx2 = x[2] - x[1];
      d = x[2] - x[0];

      //amp
      slope1 = (y[1] - y[0])/dx1;
      slope2 = (y[2] - y[1])/dx2;

      b[0] = ((dx1 + 2*d) * dx2 * slope1 +
                          dx1*dx1 * slope2) / d;
      diag[0] = dx2;
      ud[0] = d;
      ld[0] = 0.0;

  } else{
    dx1 = x[i] - x[i-1];
    dx2 = x[i+1] - x[i];

    //amp
    slope1 = (y[i] - y[i-1])/dx1;
    slope2 = (y[i+1] - y[i])/dx2;

    b[i] = 3.0* (dx2*slope1 + dx1*slope2);
    diag[i] = 2*(dx1 + dx2);
    ud[i] = dx1;
    ld[i] = dx2;
  }
}

/*
fill the B array on the GPU for response transfer functions.
*/

__device__
void fill_B(InterpArrayContainer *array_container, double *B, double *tvals, double *upper_diag, double *diag, double *lower_diag, int length, int num_splines, int spline_i, int i){
    int num_pars = 8;
    int lead_ind;

    // phaseRdelay
    lead_ind = spline_i*length;
    prep_splines(i, length, &B[lead_ind], &upper_diag[lead_ind], &diag[lead_ind], &lower_diag[lead_ind], tvals, array_container[spline_i].array);
}

__global__
void fill_B_wrap(InterpArrayContainer *array_container, double *tvals, double *B, double *upper_diag, double *diag, double *lower_diag, int length, int num_splines){
    for (int spline_i = blockIdx.y * blockDim.y + threadIdx.y;
         spline_i < num_splines;
         spline_i += blockDim.y * gridDim.y){

       for (int i = blockIdx.x * blockDim.x + threadIdx.x;
            i < length;
            i += blockDim.x * gridDim.x){

              fill_B(array_container, B, tvals, upper_diag, diag, lower_diag, length, num_splines, spline_i, i);

}
}
}


__device__
void fill_coefficients(int i, int length, double *dydx, double dx, double *y, double *coeff1, double *coeff2, double *coeff3){
  double slope, t, dydx_i;

  slope = (y[i+1] - y[i])/dx;

  dydx_i = dydx[i];

  t = (dydx_i + dydx[i+1] - 2*slope)/dx;

  coeff1[i] = dydx_i;
  coeff2[i] = (slope - dydx_i) / dx - t;
  coeff3[i] = t/dx;
}

/*
find spline constants based on matrix solution for response transfer functions.
*/
__device__
void set_spline_constants(InterpArrayContainer *array_container, double *B, int length, int num_splines, int spline_i, int i, double dt){

    int lead_ind;

    // phaseRdelay
    lead_ind = spline_i*length;
    fill_coefficients(i, length, &B[lead_ind], dt, array_container[spline_i].array, array_container[spline_i].coeff_1, array_container[spline_i].coeff_2, array_container[spline_i].coeff_3);

}

__global__
void set_spline_constants_wrap(InterpArrayContainer *array_container, double *B, int length, int num_splines, double *tvals){
    int num_pars = 8;
    int spline_index;
    double dt;

    for (int spline_i = blockIdx.y * blockDim.y + threadIdx.y;
         spline_i < num_splines;
         spline_i += blockDim.y * gridDim.y){

       for (int i = blockIdx.x * blockDim.x + threadIdx.x;
            i < length-1;
            i += blockDim.x * gridDim.x){

              dt = tvals[i + 1] - tvals[i];
              set_spline_constants(array_container, B, length, num_splines, spline_i, i, dt);
  }
  }
}

void fit_constants_serial_wrap(int m, int n, double *a, double *b, double *c, double *d_in){

  void *pBuffer;
  cusparseStatus_t stat;
  cusparseHandle_t handle;

  size_t bufferSizeInBytes;

  CUSPARSE_CALL(cusparseCreate(&handle));
  CUSPARSE_CALL( cusparseDgtsv2StridedBatch_bufferSizeExt(handle, m, a, b, c, d_in, n, m, &bufferSizeInBytes));
  gpuErrchk_here(cudaMalloc(&pBuffer, bufferSizeInBytes));

    CUSPARSE_CALL(cusparseDgtsv2StridedBatch(handle,
                                              m,
                                              a, // dl
                                              b, //diag
                                              c, // du
                                              d_in,
                                              n,
                                              m,
                                              pBuffer));


CUSPARSE_CALL(cusparseDestroy(handle));
gpuErrchk_here(cudaFree(pBuffer));
}


void Interpolate::setup(InterpArrayContainer *array_container, double *d_tvec, int m_, int n_){

    m = m_;
    n = n_;

    int NUM_THREADS = 256;

    int num_blocks = std::ceil((m + NUM_THREADS -1)/NUM_THREADS);

    dim3 interpGrid(num_blocks, n);

    fill_B_wrap<<<interpGrid, NUM_THREADS>>>(array_container, d_tvec, d_B, d_du, d_d, d_dl, m, n);
    cudaDeviceSynchronize();
    gpuErrchk_here(cudaGetLastError());

    fit_constants_serial_wrap(m, n, d_dl, d_d, d_du, d_B);

    set_spline_constants_wrap<<<interpGrid, NUM_THREADS>>>(array_container, d_B, m, n, d_tvec);
    cudaDeviceSynchronize();
    gpuErrchk_here(cudaGetLastError());
}

__host__ Interpolate::~Interpolate(){
    cudaFree(d_dl);
    cudaFree(d_du);
    cudaFree(d_d);
    cudaFree(d_B);
    //delete[] d;
    //delete[] dl;
    //delete[] du;
}
