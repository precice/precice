#include "mapping/GinkgoKernels.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "math/math.hpp"

#include <functional>
#include <ginkgo/ginkgo.hpp>
#include <ginkgo/extensions.hpp>

using precice::mapping::RadialBasisParameters;
using precice::math::pow_int;

template <typename T, typename MemorySpace>
using native_type = gko::ext::kokkos::native_type<T, MemorySpace>;

using gko::ext::kokkos::map_data;

namespace precice::mapping::kernel {

template <typename EvalFunctionType>
void create_rbf_system_matrix(std::shared_ptr<const gko::Executor> exec,
                              gko::ptr_param<GinkgoMatrix> mtx, const std::array<bool, 3> activeAxis,
                              gko::ptr_param<GinkgoMatrix> supportPoints, gko::ptr_param<GinkgoMatrix> targetPoints,
                              EvalFunctionType &&f, ::precice::mapping::RadialBasisParameters rbf_params, bool addPolynomial, unsigned int extraDims)
{
  exec->run(
      gko::ext::kokkos::parallel_for(
          "create_rbf_system_matrix",
          gko::ext::kokkos::make_policy_top<Kokkos::MDRangePolicy, Kokkos::Rank<2>>(Kokkos::Array<int64_t, 2>{0, 0},
                                                                                    Kokkos::Array<int64_t, 2>{static_cast<int64_t>(mtx->get_size()[0]),
                                                                                                              static_cast<int64_t>(mtx->get_size()[1])}),
          [=] GKO_KOKKOS_FN(const int &i, const int &j, auto k_mtx, auto k_supportPoints, auto k_targetPoints) {
            double dist = 0;

            // Make each entry zero if polynomial is on since not every entry will be adjusted below
            if (addPolynomial) {
              k_mtx(i, j) = 0;
            }

            // Loop over each dimension and calculate euclidean distance
            for (size_t k = 0; k < k_supportPoints.values.extent(1); ++k) {
              dist += pow_int<2>(k_supportPoints(j, k) - k_targetPoints(i, k)) * static_cast<int>(activeAxis[k]);
            }

            dist = Kokkos::sqrt(dist);

            k_mtx(i, j) = f(dist, rbf_params);
          },
          mtx.get(), supportPoints.get(), targetPoints.get()));
}

#define INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(_function_type)                                                                                     \
  template void create_rbf_system_matrix<_function_type>(std::shared_ptr<const gko::Executor> exec,                                             \
                                                         gko::ptr_param<GinkgoMatrix> mtx, const std::array<bool, 3> activeAxis,                \
                                                         gko::ptr_param<GinkgoMatrix> supportPoints, gko::ptr_param<GinkgoMatrix> targetPoints, \
                                                         _function_type &&f, RadialBasisParameters rbf_params, bool addPolynomial, unsigned int extraDims)

INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(ThinPlateSplines);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(Multiquadrics);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(InverseMultiquadrics);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(VolumeSplines);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(Gaussian);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactThinPlateSplinesC2);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC0);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC2);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC4);
//INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC6);

} // namespace precice::mapping::kernel

// template <typename ValueType, typename EvalFunctionType>
// void create_rbf_system_matrix(std::shared_ptr<const DefaultExecutor> exec,
//                               const std::size_t n1, const std::size_t n2, const std::size_t dataDimensionality, const std::array<bool, 3> activeAxis, ValueType *mtx, ValueType *supportPoints,
//                               ValueType *targetPoints, EvalFunctionType f, const RadialBasisParameters rbf_params, const std::size_t inputRowLength, const std::size_t outputRowLength, const bool addPolynomial, const unsigned int extraDims = 0)
//{
//   run_kernel(
//       exec,
//       GKO_KERNEL(auto i, auto j, auto N, auto dataDimensionality, auto activeAxis, auto mtx, auto supportPoints, auto targetPoints, auto f, auto rbf_params, auto inputRowLength, auto outputRowLength, auto addPolynomial, auto extraDims) {
//         const unsigned int rowLength = N + extraDims;
//         double             dist      = 0;
//
//         // Make each entry zero if polynomial is on since not every entry will be adjusted below
//         if (addPolynomial) {
//           mtx[i * rowLength + j] = 0;
//         }
//
// #if defined(__NVCC__) || defined(__HIPCC__)
//
//         double y;
//         for (size_t k = 0; k < dataDimensionality; ++k) {
//           y    = supportPoints[k * inputRowLength + j] - targetPoints[k * outputRowLength + i];
//           dist = fma(y, y, dist);
//         }
//
//         dist = sqrt(dist);
//
// #else
//         const unsigned int supportPointOffset = dataDimensionality * j; // Point of current column
//         const unsigned int targetPointOffset  = dataDimensionality * i; // Point of current row
//         // Loop over each dimension and calculate euclidean distance
//         for (size_t k = 0; k < dataDimensionality; ++k) {
//           dist += pow_int<2>(supportPoints[supportPointOffset + k] - targetPoints[targetPointOffset + k]) * static_cast<int>(activeAxis.at(k));
//         }
//
//         dist                                  = std::sqrt(dist);
// #endif
//
//         mtx[i * rowLength + j] = f(dist, rbf_params);
//       },
//       gko::dim<2>{n1, n2}, n2, dataDimensionality, activeAxis, mtx, supportPoints, targetPoints, f, rbf_params, inputRowLength, outputRowLength, addPolynomial, extraDims);
// }
//
//// Here, we need to instantiate all possible variants for each basis function
//
// template void create_rbf_system_matrix<double, precice::mapping::ThinPlateSplines>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                   double *, double *, double *, precice::mapping::ThinPlateSplines, const RadialBasisParameters,
//                                                                                   const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::Multiquadrics>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                double *, double *, double *, precice::mapping::Multiquadrics, const RadialBasisParameters,
//                                                                                const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::InverseMultiquadrics>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                       double *, double *, double *, precice::mapping::InverseMultiquadrics, const RadialBasisParameters,
//                                                                                       const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::VolumeSplines>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                double *, double *, double *, precice::mapping::VolumeSplines, const RadialBasisParameters,
//                                                                                const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::Gaussian>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                           double *, double *, double *, precice::mapping::Gaussian, const RadialBasisParameters,
//                                                                           const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::CompactThinPlateSplinesC2>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                            double *, double *, double *, precice::mapping::CompactThinPlateSplinesC2, const RadialBasisParameters,
//                                                                                            const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC0>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                      double *, double *, double *, precice::mapping::CompactPolynomialC0, const RadialBasisParameters,
//                                                                                      const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC2>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                      double *, double *, double *, precice::mapping::CompactPolynomialC2, const RadialBasisParameters,
//                                                                                      const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC4>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                      double *, double *, double *, precice::mapping::CompactPolynomialC4, const RadialBasisParameters,
//                                                                                      const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template void create_rbf_system_matrix<double, precice::mapping::CompactPolynomialC6>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, const std::size_t, const std::array<bool, 3>,
//                                                                                      double *, double *, double *, precice::mapping::CompactPolynomialC6, const RadialBasisParameters,
//                                                                                      const std::size_t, const std::size_t, const bool, const unsigned int);
//
// template <typename ValueType>
// void fill_polynomial_matrix(std::shared_ptr<const DefaultExecutor> exec,
//                            const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *x, const std::size_t supportPointsRowLength, const unsigned int dims = 4)
//{
//  run_kernel(
//      exec,
//      GKO_KERNEL(auto i, auto j, auto N1, auto N2, auto mtx, auto x, auto supportPointsRowLength, auto dims) {
// #if defined(__NVCC__) || defined(__HIPCC__)
//        if (j < dims - 1) {
//          mtx[i * dims + j] = x[j * supportPointsRowLength + i];
//        } else {
//          mtx[i * dims + j] = 1;
//        }
// #else
//        const unsigned int supportPointOffset = (dims - 1) * i;
//        if (j < dims - 1) {
//          mtx[i * dims + j] = x[supportPointOffset + j];
//        } else {
//          mtx[i * dims + j] = 1;
//        }
// #endif
//      },
//      gko::dim<2>{n1, n2}, n1, n2, mtx, x, supportPointsRowLength, dims);
//}
//
// template void fill_polynomial_matrix<double>(std::shared_ptr<const DefaultExecutor>, const std::size_t, const std::size_t, double *, double *, const std::size_t, const unsigned int);
