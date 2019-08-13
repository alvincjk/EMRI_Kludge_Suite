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

InterpArrayContainer * createInterpArrayContainer_gpu(size_t numBytes){
    InterpArrayContainer *gpu_array_container;
    cudaMalloc((void**)&gpu_array_container, numBytes);
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


__global__
void fill_B(InterpArrayContainer *arr_container, double *B, int length_per_arr, int num_arr){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    if (j >= num_arr) return;
    if (i >= length_per_arr) return;

            if (i == length_per_arr - 1){
                B[j*length_per_arr + i] = 3.0* (arr_container[j].array[i] - arr_container[j].array[(i-1)]);

            } else if (i == 0){
                B[j*length_per_arr + i] = 3.0* (arr_container[j].array[1] - arr_container[j].array[0]);

            } else{
                B[j*length_per_arr + i] = 3.0* (arr_container[j].array[(i+1)] - arr_container[j].array[(i-1)]);
            }
        /*# if __CUDA_ARCH__>=200
        if ((i < 100) && (j ==8))
            printf("%d %d, %.18e, %.18e, %.18e, %.18e\n", i, j, B[j*length_per_arr + i], arr_container[j].array[i+1], arr_container[j].array[i], arr_container[j].array[i-1]);
        #endif //*/
}

__global__
void set_spline_constants(InterpArrayContainer *arr_container, double *B, int length_per_arr, int num_arr){
    double D_i, D_ip1, y_i, y_ip1;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    if (j >= num_arr) return;
    if (i >= length_per_arr) return;

            D_i = B[(j*length_per_arr) + i];
            D_ip1 = B[(j*length_per_arr) + i + 1];
            y_i = arr_container[j].array[i];
            y_ip1 = arr_container[j].array[i+1];
            arr_container[j].coeff_1[i] = D_i;
            arr_container[j].coeff_2[i] = 3.0 * (y_ip1 - y_i) - 2.0*D_i - D_ip1;
            arr_container[j].coeff_3[i] = 2.0 * (y_i - y_ip1) + D_i + D_ip1;

            /*# if __CUDA_ARCH__>=200
            if ((i % 2000 == 0) && (j == 3))
                printf("%d, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e, %.18e \n", i, B[(j*length_per_arr) + i], B[(j*length_per_arr) + i+1], arr_container[j].array[i], arr_container[j].array[i], arr_container[j].array[i+1], arr_container[j].coeff_1[i], arr_container[j].coeff_2[i], arr_container[j].coeff_3[i]);
            #endif //*/
}



Interpolate::Interpolate(){
    int pass = 0;
}

__host__
void Interpolate::alloc_arrays(int max_length_init, int num_arr){
    gpuErrchk_here(cudaMalloc(&d_B, max_length_init*num_arr*sizeof(double)));
    gpuErrchk_here(cudaMalloc(&d_dl, max_length_init*sizeof(double)));
    gpuErrchk_here(cudaMalloc(&d_d, max_length_init*sizeof(double)));
    gpuErrchk_here(cudaMalloc(&d_du, max_length_init*sizeof(double)));
}

__global__
void setup_d_vals(double *dl, double *d, double *du, int current_length){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= current_length) return;
    if (i == 0){
        dl[0] = 0.0;
        d[0] = 2.0;
        du[0] = 1.0;
    } else if (i == current_length - 1){
        dl[current_length-1] = 1.0;
        d[current_length-1] = 2.0;
        du[current_length-1] = 0.0;
    } else{
        dl[i] = 1.0;
        d[i] = 4.0;
        du[i] = 1.0;
    }
    # /*if __CUDA_ARCH__>=200
    if ((i == 0) || (i == current_length-1) || (i == 10))
        printf("%d, %e, %e, %e \n", i, dl[i], d[i], du[i]);
    #endif //*/
}

/*__device__
void find_index_and_xout(int *index, double *x_out, double dx, double x_new, double *x_old, int length){
    *index = floor((x_new - x_old[0])/dx);
    if (*index >= length - 1) *index = length - 2;
    *x_out = (x_new - x_old[*index])/(x_old[*index+1] - x_old[*index]);
}

__device__
double interpolate_array(InterpArrayContainer array_container, double x, int index){
    double coeff_0 = array_container.array[index];
    double coeff_1 = array_container.coeff_1[index];
    double coeff_2 = array_container.coeff_2[index];
    double coeff_3 = array_container.coeff_3[index];
    double x2 = x*x;
    double x3 = x*x2;
    double return_val = coeff_0 + coeff_1*x + coeff_2*x2 + coeff_3*x3;
    return return_val;
}//*/

void Interpolate::setup(InterpArrayContainer *array_container, int m_, int n_){
    m = m_;
    n = n_;


    int NUM_THREADS = 256;

    int num_blocks = std::ceil((m + NUM_THREADS -1)/NUM_THREADS);

    dim3 interpGrid(num_blocks, n);

    setup_d_vals<<<num_blocks, NUM_THREADS>>>(d_dl, d_d, d_du, m);
    cudaDeviceSynchronize();
    gpuErrchk_here(cudaGetLastError());

    fill_B<<<interpGrid, NUM_THREADS>>>(array_container, d_B, m, n);
    cudaDeviceSynchronize();
    gpuErrchk_here(cudaGetLastError());

    /*cudaMemcpy(checker, d_B, m*n*sizeof(double), cudaMemcpyDeviceToHost);
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j+=100){
            printf("%d %d, %e\n", i, j, checker[i*m + j]);
        }
    }//*/

    double *checker = new double[m*n];

    CUSPARSE_CALL( cusparseCreate(&handle) );
    cusparseStatus_t status = cusparseDgtsv_nopivot(handle, m, n, d_dl, d_d, d_du, d_B, m);
    if (status !=  CUSPARSE_STATUS_SUCCESS) assert(0);
    cusparseDestroy(handle);

    set_spline_constants<<<interpGrid, NUM_THREADS>>>(array_container, d_B, m, n);
    cudaDeviceSynchronize();
    gpuErrchk_here(cudaGetLastError());

    delete[] checker;
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
