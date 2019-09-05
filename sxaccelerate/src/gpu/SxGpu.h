#ifndef _SX_GPU_H_
#define _SX_GPU_H_

void sx_gpu_gemmm3 (ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *A,
              const SxComplex16 *B,
              const double *C,
              SxComplex16 *res, ssize_t rldc, ssize_t rldbc);

void sx_gpu_gemm3m (ssize_t Na, ssize_t Nb, ssize_t Nc, ssize_t Nd,
              const SxComplex16 *B,
              const double *C,
              const SxComplex16 *X,
              SxComplex16 *res, ssize_t xldc, ssize_t xldbc);

void sx_gpu_gemmm_free ();

void sx_gpu_set_device (int id);

#endif
