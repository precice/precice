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

#include <ginkgo/ginkgo.hpp>
#include <stdio.h>

#include <ginkgo/kernels/kernel_launch.hpp>

namespace GKO_DEVICE_NAMESPACE {

using namespace gko::kernels::GKO_DEVICE_NAMESPACE;
using vec = gko::matrix::Dense<double>;
using mat = gko::matrix::Dense<double>;

template <typename ValueType>
void create_rbf_system_matrix(std::shared_ptr<const DefaultExecutor> exec,
                              const std::size_t n1, const std::size_t n2, const ValueType rbf_shape,
                              const ValueType support_radius, ValueType *mtx, ValueType *support_points,
                              ValueType *target_points, const bool add_polynomial, const unsigned int extra_dims = 0)
{
  run_kernel(
      exec,
      GKO_KERNEL(auto i, auto j, auto N, auto mtx, auto support_points, auto target_points, auto support_radius, auto rbf_shape, auto add_polynomial, auto extra_dims) {
        const unsigned int dim                  = 3;
        const unsigned int row_length           = N + extra_dims;
        const unsigned int center_coords_offset = 3 * j;
        const unsigned int eval_coords_offset   = 3 * i;
        float              dist                 = 0;

        // Make each entry zero if polynomial is on since not every entry will be adjusted below
        if (add_polynomial) {
          mtx[i * row_length + j] = 0;
        }

        for (size_t k = 0; k < dim; ++k) {
          dist += (support_points[center_coords_offset + k] - target_points[eval_coords_offset + k]) * (support_points[center_coords_offset + k] - target_points[eval_coords_offset + k]);
        }

        dist         = sqrt(dist);
        float result = 0.f;
        if (dist <= support_radius) {
          result = exp(-(rbf_shape * dist) * (rbf_shape * dist));
        }
        mtx[i * row_length + j] = result;

        // Give first column responsibility for adding additional polynomial part (TODO: Optimize for GPU)
        if (add_polynomial) {
          if (j == 0) {
            for (std::size_t k = 0; k < extra_dims - 1; ++k) {
              mtx[i * row_length + N + k]              = target_points[eval_coords_offset + k];
              mtx[N * row_length + k * row_length + i] = target_points[eval_coords_offset + k];
            }
            mtx[i * row_length + N + extra_dims - 1] = 1;
            mtx[(row_length - 1) * row_length + i]   = 1;
          }

          // Fill lower right block of 0's (PAY ATTENTION TO BLOCKS!)
          if (i == 0 && j < extra_dims) {
            for (std::size_t k = 0; k < extra_dims; ++k) {
              mtx[N * row_length + j * row_length + N + k] = 0;
            }
          }
        }
      },
      gko::dim<2>{n1, n2}, n2, mtx, support_points, target_points, support_radius, rbf_shape, add_polynomial, extra_dims);
}

template void create_rbf_system_matrix<double>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t,
                                               const double, const double, double *, double *, double *, const bool, const unsigned int);

template <typename ValueType>
void fill_polynomial_matrix(std::shared_ptr<const DefaultExecutor> exec,
                            const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *x, const unsigned int dims = 4)
{
  run_kernel(
      exec,
      GKO_KERNEL(auto i, auto j, auto N1, auto N2, auto mtx, auto x, auto dims) {
        const unsigned int center_coords_offset = 3 * i;
        if (j < dims - 1) {
          mtx[i * dims + j] = x[center_coords_offset + j];
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
