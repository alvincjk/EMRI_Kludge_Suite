#ifndef __INTERPOLATE_H_
#define __INTERPOLATE_H_
#include <cusparse_v2.h>


typedef struct tagInterpArrayContainer{
    double *array;
    double *coeff_1;
    double *coeff_2;
    double *coeff_3;
} InterpArrayContainer;

class Interpolate{
    double *w;
    double *D;

    double *dl;
    double *d;
    double *du;
    double *d_dl;
    double *d_d;
    double *d_du;
    cusparseHandle_t  handle;
    cudaError_t err;
    int m;
    int n;
    double *d_B;

public:
    // FOR NOW WE ASSUME dLOGX is evenly spaced // TODO: allocate at the beginning
    Interpolate();

    __host__ void alloc_arrays(int max_length_init, int num_arr);
    __host__ void setup(InterpArrayContainer *array_container, double *d_tvec, int m_, int n_);

    __host__ ~Interpolate(); //destructor
};

void destroyInterpArrayContainer(InterpArrayContainer *array_container, int num_arr);

InterpArrayContainer * createInterpArrayContainer(int num_arr, int num_points);

__global__
void set_spline_constants(InterpArrayContainer *arr_container, double *B, int length_per_arr, int num_arr);

__global__
void fill_B(InterpArrayContainer *arr_container, double *B, int length_per_arr, int num_arr);

/*extern __device__
double interpolate_array(InterpArrayContainer array_container, double x, int index);

extern __device__
void find_index_and_xout(int *index, double *x_out, double dx, double x_new, double *x_old, int length);//*/

#endif //__INTERPOLATE_H_
