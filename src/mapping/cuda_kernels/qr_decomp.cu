#include "ginkgo/ginkgo.hpp"
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include "cuda_runtime.h"
#include <cuda.h>
#include "device_launch_parameters.h"
#include <iostream>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include "utils/Event.hpp"
#include "utils/EventUtils.hpp"

using GinkgoMatrix = gko::matrix::Dense<>;

// Handles for low-level CUDA libraries
cusolverDnHandle_t solverHandle;
cublasHandle_t cublasHandle;
cublasStatus_t cublasStatus = CUBLAS_STATUS_SUCCESS;
cusolverStatus_t cusolverStatus = CUSOLVER_STATUS_SUCCESS;
cudaError_t cudaErrorCode = cudaSuccess;

// Important variables which track the state of the solver routines
double *dTau = nullptr;
double *dWork = nullptr;
int *devInfo = nullptr;

void initCuSolver(){
    // Allocating important CUDA variables
    cudaMalloc((void **)&dWork, sizeof(double));
    cudaMalloc((void **)&devInfo, sizeof(int));
    cudaMalloc((void **)&dTau, sizeof(double));

}

void deInitCuSolver(){
    // Freeing CUDA variables
    cudaFree(dTau);
    cudaFree(dWork);
    cudaFree(devInfo);
}

void computeQR(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoMatrix *R)
{
    cusolverDnCreate(&solverHandle);
    cublasCreate(&cublasHandle);

    // NOTE: It's important to transpose since cuSolver assumes column-major memory layout
    // Making a copy since every value will be overridden
    auto A_T = gko::share(GinkgoMatrix::create(exec, gko::dim<2>(A_Q->get_size()[1], A_Q->get_size()[0])));
    A_Q->transpose(gko::lend(A_T));

    // Setting dimensions for solver
    const unsigned int M = A_T->get_size()[1];
    const unsigned int N = A_T->get_size()[0];

    const int lda = max(1, M);
    const int k = min(M, N);

    int lwork_geqrf = 0;
    int lwork_orgqr = 0;
    int lwork = 0;

    precice::utils::Event calculateQRDecompEvent{"calculateQRDecomp"};

    // Query working space of geqrf and orgqr
    cusolverStatus = cusolverDnDgeqrf_bufferSize(solverHandle, M, N, A_T->get_values(), lda, &lwork_geqrf);
    assert(cusolverStatus == cusolverStatus_SUCCESS);
    cusolverStatus = cusolverDnDorgqr_bufferSize(solverHandle, M, N, k, A_T->get_values(), lda, dTau, &lwork_orgqr);
    assert(cusolverStatus == cusolverStatus_SUCCESS);
    lwork = (lwork_geqrf > lwork_orgqr) ? lwork_geqrf : lwork_orgqr;
    cudaErrorCode = cudaMalloc((void **)&dWork, sizeof(double) * lwork);
    assert(cudaSuccess == cudaErrorCode);

    // Compute QR factorization
    cusolverStatus = cusolverDnDgeqrf(solverHandle, M, N, A_T->get_values(), lda, dTau, dWork, lwork, devInfo);
    cudaErrorCode = cudaDeviceSynchronize();
    assert(cusolverStatus_SUCCESS == cusolverStatus);
    assert(cudaSuccess == cudaErrorCode);

    // Copy A_T to R s.t. the upper triangle corresponds to R
    A_T->transpose(gko::lend(R));

    // Compute Q
    cusolverStatus = cusolverDnDorgqr(solverHandle, M, N, k, A_T->get_values(), lda, dTau, dWork, lwork, devInfo);
    cudaErrorCode = cudaDeviceSynchronize();
    assert(cusolverStatus_SUCCESS == cusolverStatus);
    assert(cudaSuccess == cudaErrorCode);

    A_T->transpose(gko::lend(A_Q));

    cudaDeviceSynchronize();

    calculateQRDecompEvent.stop();

    return;
}
