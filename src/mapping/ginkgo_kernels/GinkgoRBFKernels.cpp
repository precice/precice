/*******************************<GINKGO LICENSE>******************************
Copyright (c) 2017-2021, the Ginkgo authors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************<GINKGO LICENSE>*******************************/

#include "mapping/impl/DeviceBasisFunctions.cuh"

#include <functional>
#include <ginkgo/ginkgo.hpp>
#include <stdio.h>

#include <ginkgo/kernels/kernel_launch.hpp>

namespace GKO_DEVICE_NAMESPACE {

using namespace gko::kernels::GKO_DEVICE_NAMESPACE;
using vec = gko::matrix::Dense<double>;
using mat = gko::matrix::Dense<double>;

template <typename ValueType, typename EvalFunctionType>
void create_rbf_system_matrix(std::shared_ptr<const DefaultExecutor> exec,
                              const std::size_t n1, const std::size_t n2, const std::size_t dataDimensionality, const std::array<bool, 3> activeAxis, ValueType *mtx, ValueType *supportPoints,
                              ValueType *targetPoints, EvalFunctionType f, const std::array<ValueType, 3> rbf_params, const bool addPolynomial, const unsigned int extraDims = 0)
{
  run_kernel(
      exec,
      GKO_KERNEL(auto i, auto j, auto N, auto dataDimensionality, auto activeAxis, auto mtx, auto supportPoints, auto targetPoints, auto f, auto rbf_params, auto addPolynomial, auto extraDims) {
        const unsigned int rowLength          = N + extraDims;
        const unsigned int supportPointOffset = dataDimensionality * j; // Point of current column
        const unsigned int evalPointOffset    = dataDimensionality * i; // Point of current row
        double             dist               = 0;

        // Make each entry zero if polynomial is on since not every entry will be adjusted below
        if (addPolynomial) {
          mtx[i * rowLength + j] = 0;
        }

#ifdef __NVCC__

        // Use this as readonly shared buffer
        __shared__ double prefetchedEvalPoint[3];
        double            y;

        // Check if current block is at the end of a matrix row and induces a line break
        // However, if a block is larger than an entire row, we need to disable it since we now can't make sure
        // it does not still span across two lines.
        // We have to use uint64_t because these matrices become so large that using this
        // if-condition overflows with int32
        uint32_t blockID        = blockIdx.x;
        uint64_t leftThreadNum  = blockID * blockDim.x;
        uint64_t rightThreadNum = fma(blockID + 1, blockID, -1); // static_cast<uint64_t>(blockIdx.x + 1) * blockDim.x - 1;
        if (blockDim.x < rowLength && leftThreadNum % rowLength > rightThreadNum % rowLength) {

          // Since this block spans across two lines, we have to use global memory and cannot use prefetched memory
          for (size_t k = 0; k < dataDimensionality; ++k) {
            y    = supportPoints[supportPointOffset + k] - targetPoints[evalPointOffset + k];
            dist = fma(y, y, dist);
          }
        } else {

          // If this block is indeed only in one row, we can make thread 0 in each block responsible for prefetching values into shared memory
          if (0 == threadIdx.x) {
            prefetchedEvalPoint[0] = targetPoints[evalPointOffset];
            prefetchedEvalPoint[1] = targetPoints[evalPointOffset + 1];
            prefetchedEvalPoint[2] = targetPoints[evalPointOffset + 2];
          }
          // Let all threads in a block wait until memory is prefetched
          __syncthreads();

          for (size_t k = 0; k < dataDimensionality; ++k) {
            y    = supportPoints[supportPointOffset + k] - prefetchedEvalPoint[k];
            dist = fma(y, y, dist); //(supportPoints[supportPointOffset + k] - prefetchedEvalPoint[k]) * (supportPoints[supportPointOffset + k] - prefetchedEvalPoint[k]) * static_cast<int>(activeAxis.at(k));
          }
        }

        dist = sqrt(dist);

#else
        // Loop over each dimension and calculate euclidian distance
        for (size_t k = 0; k < dataDimensionality; ++k) {
          dist += std::pow(supportPoints[supportPointOffset + k] - targetPoints[evalPointOffset + k], 2) * static_cast<int>(activeAxis.at(k));
        }

        dist = std::sqrt(dist);
#endif

        mtx[i * rowLength + j] = f(dist, rbf_params);
      },
      gko::dim<2>{n1, n2}, n2, dataDimensionality, activeAxis, mtx, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
}

// Here, we need to instantiate all possible variants for each basis function

template void create_rbf_system_matrix<double, precice::mapping::ThinPlateSplinesFunctor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                          double *, double *, double *, precice::mapping::ThinPlateSplinesFunctor, const std::array<double, 3>,
                                                                                          const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::MultiQuadraticsFunctor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                         double *, double *, double *, precice::mapping::MultiQuadraticsFunctor, const std::array<double, 3>,
                                                                                         const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::InverseMultiquadricsFunctor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                              double *, double *, double *, precice::mapping::InverseMultiquadricsFunctor, const std::array<double, 3>,
                                                                                              const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::VolumeSplinesFunctor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                       double *, double *, double *, precice::mapping::VolumeSplinesFunctor, const std::array<double, 3>,
                                                                                       const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::GaussianFunctor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                  double *, double *, double *, precice::mapping::GaussianFunctor, const std::array<double, 3>,
                                                                                  const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::CompactThinPlateSplinesC2Functor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                                   double *, double *, double *, precice::mapping::CompactThinPlateSplinesC2Functor, const std::array<double, 3>,
                                                                                                   const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC0Functor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                             double *, double *, double *, precice::mapping::CompactPolynomialC0Functor, const std::array<double, 3>,
                                                                                             const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC2Functor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                             double *, double *, double *, precice::mapping::CompactPolynomialC2Functor, const std::array<double, 3>,
                                                                                             const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC4Functor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                             double *, double *, double *, precice::mapping::CompactPolynomialC4Functor, const std::array<double, 3>,
                                                                                             const bool, const unsigned int);

template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC6Functor>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
                                                                                             double *, double *, double *, precice::mapping::CompactPolynomialC6Functor, const std::array<double, 3>,
                                                                                             const bool, const unsigned int);

template <typename ValueType>
void fill_polynomial_matrix(std::shared_ptr<const DefaultExecutor> exec,
                            const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *x, const unsigned int dims = 4)
{
  run_kernel(
      exec,
      GKO_KERNEL(auto i, auto j, auto N1, auto N2, auto mtx, auto x, auto dims) {
        const unsigned int supportPointOffset = 3 * i;
        if (j < dims - 1) {
          mtx[i * dims + j] = x[supportPointOffset + j];
        } else {
          mtx[i * dims + j] = 1;
        }
      },
      gko::dim<2>{n1, n2}, n1, n2, mtx, x, dims);
}

template void fill_polynomial_matrix<double>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, double *, double *, const unsigned int);

template <typename ValueType>
void extract_upper_triangular(std::shared_ptr<const DefaultExecutor> exec, ValueType *src, ValueType *dest, const std::size_t i, const std::size_t j, const std::size_t N)
{
  run_kernel(
      exec,
      GKO_KERNEL(auto i, auto j, auto src, auto dest, auto N) {
        dest[i * N + j] = src[i * N + j] * (int) (j >= i);
      },
      gko::dim<2>{i, j}, src, dest, N);
}

template void extract_upper_triangular<double>(std::shared_ptr<const DefaultExecutor>, double *, double *, const std::size_t, const std::size_t, const std::size_t);

} // namespace GKO_DEVICE_NAMESPACE
