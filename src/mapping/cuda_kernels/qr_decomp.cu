#ifdef GINKGO_BUILD_CUDA
#include <ginkgo/ginkgo.hpp>
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "device_launch_parameters.h"
#include <cuda_runtime.h>
#include "utils/Event.hpp"
#include "utils/EventUtils.hpp"

using GinkgoMatrix = gko::matrix::Dense<>;

// Handles for low-level CUDA libraries
cusolverDnHandle_t solverHandle;
cusolverStatus_t cusolverStatus = CUSOLVER_STATUS_SUCCESS;
cudaError_t cudaErrorCode = cudaSuccess;

// Important variables which track the state of the solver routines
double *dTau = nullptr;
void *dWork = nullptr;
void *hWork = nullptr;
int *devInfo = nullptr;

int cudaBackupDeviceId = 0;

void initQRSolver(const int deviceId=0){
    cudaGetDevice(&cudaBackupDeviceId);
    cudaSetDevice(deviceId);

    // Allocating important CUDA variables
    cudaMalloc((void **)&dWork, sizeof(double));
    cudaMalloc((void **)&devInfo, sizeof(int));

    cusolverDnCreate(&solverHandle);

}

void deInitQRSolver(){
    // Freeing CUDA variables
    cudaFree(dTau);
    cudaFree(dWork);
    cudaFree(devInfo);

    if (nullptr != solverHandle){
        cusolverDnDestroy(solverHandle);
    }
}

void computeQR(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoMatrix *R)
{
    cusolverDnCreate(&solverHandle);

    // NOTE: It's important to transpose since cuSolver assumes column-major memory layout
    // Making a copy since every value will be overridden
    auto A_T = gko::share(GinkgoMatrix::create(exec, gko::dim<2>(A_Q->get_size()[1], A_Q->get_size()[0])));
    A_Q->transpose(gko::lend(A_T));

    // Setting dimensions for solver
    const unsigned int M = A_T->get_size()[1];
    const unsigned int N = A_T->get_size()[0];

    const int lda = max(1, M);
    const int k = min(M, N);

    size_t dLwork_geqrf = 0;
    size_t dLwork_orgqr = 0;
    size_t dLwork = 0;

    size_t hLwork_geqrf = 0;
    size_t hLwork = 0;

    cudaMalloc((void **)&dTau, sizeof(double) * M);

    precice::utils::Event calculateQRDecompEvent{"calculateQRDecomp"};

    // Query working space of geqrf and orgqr
    cusolverStatus = cusolverDnXgeqrf_bufferSize(solverHandle, nullptr, M, N, CUDA_R_64F, A_T->get_values(), lda, CUDA_R_64F, dTau, CUDA_R_64F, &dLwork_geqrf, &hLwork_geqrf);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
    cusolverStatus = cusolverDnDorgqr_bufferSize(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (int*) &dLwork_orgqr);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
    dLwork = (dLwork_geqrf > dLwork_orgqr) ? dLwork_geqrf : dLwork_orgqr;
    cudaErrorCode = cudaMalloc((void **)&dWork, sizeof(double) * dLwork);
    assert(cudaSuccess == cudaErrorCode);

    // Compute QR factorization
    cusolverStatus = cusolverDnXgeqrf(solverHandle, nullptr, M, N, CUDA_R_64F, A_T->get_values(), lda, CUDA_R_64F, dTau, CUDA_R_64F, dWork, dLwork, hWork, hLwork, devInfo);
    cudaErrorCode = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolverStatus);
    assert(cudaSuccess == cudaErrorCode);

    // Copy A_T to R s.t. the upper triangle corresponds to R
    A_T->transpose(gko::lend(R));

    // Compute Q
    cusolverStatus = cusolverDnDorgqr(solverHandle, M, N, k, A_T->get_values(), lda, dTau, (double*) dWork, dLwork, devInfo);
    cudaErrorCode = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolverStatus);
    assert(cudaSuccess == cudaErrorCode);

    A_T->transpose(gko::lend(A_Q));

    cudaDeviceSynchronize();

    calculateQRDecompEvent.stop();

    cudaSetDevice(cudaBackupDeviceId); // Switch back to the GPU used for all coupled solvers

    return;
}
#endif
