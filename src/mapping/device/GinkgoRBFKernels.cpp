#include "mapping/GinkgoKernels.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "math/math.hpp"

#include <functional>
#include <ginkgo/extensions.hpp>
#include <ginkgo/ginkgo.hpp>

using precice::mapping::RadialBasisParameters;
using precice::math::pow_int;

template <typename T, typename MemorySpace>
using native_type = gko::ext::kokkos::native_type<T, MemorySpace>;

using gko::ext::kokkos::map_data;

namespace precice::mapping::kernel {

auto create_row_major_acc()
{
  return [](const int &i, const int &j, auto view) {
    return view(i, j);
  };
}

auto create_col_major_acc()
{
  return [](const int &i, const int &j, auto view) {
    return view(j, i);
  };
};

template <typename EvalFunctionType>
void create_rbf_system_matrix(std::shared_ptr<const gko::Executor> exec,
                              gko::ptr_param<GinkgoMatrix> mtx, const std::array<bool, 3> activeAxis,
                              gko::ptr_param<GinkgoMatrix> supportPoints, gko::ptr_param<GinkgoMatrix> targetPoints,
                              EvalFunctionType f, ::precice::mapping::RadialBasisParameters rbf_params, bool addPolynomial, unsigned int extraDims)
{
  auto run_kernel = [=](auto acc) {
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
              for (size_t k = 0; k < activeAxis.size(); ++k) {
                if (activeAxis[k]) {
                  dist += pow_int<2>(acc(j, k, k_supportPoints) - acc(i, k, k_targetPoints));
                }
              }

              dist = Kokkos::sqrt(dist);

              k_mtx(i, j) = f(dist, rbf_params);
            },
            mtx.get(), supportPoints.get(), targetPoints.get()));
  };

  if (dynamic_cast<const gko::ReferenceExecutor *>(exec.get()) ||
      dynamic_cast<const gko::OmpExecutor *>(exec.get())) {
    run_kernel(create_row_major_acc());
  } else {
    run_kernel(create_col_major_acc());
  }
}

#define INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(_function_type)                                                                                     \
  template void create_rbf_system_matrix<_function_type>(std::shared_ptr<const gko::Executor> exec,                                             \
                                                         gko::ptr_param<GinkgoMatrix> mtx, const std::array<bool, 3> activeAxis,                \
                                                         gko::ptr_param<GinkgoMatrix> supportPoints, gko::ptr_param<GinkgoMatrix> targetPoints, \
                                                         _function_type f, RadialBasisParameters rbf_params, bool addPolynomial, unsigned int extraDims)

INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(ThinPlateSplines);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(Multiquadrics);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(InverseMultiquadrics);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(VolumeSplines);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(Gaussian);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactThinPlateSplinesC2);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC0);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC2);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC4);
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC6);

void fill_polynomial_matrix(std::shared_ptr<const gko::Executor> exec,
                            gko::ptr_param<GinkgoMatrix> mtx, gko::ptr_param<const GinkgoMatrix> x, const unsigned int dims)
{
  auto run_kernel = [=](auto acc) {
    exec->run(
        gko::ext::kokkos::parallel_for(
            "create_rbf_system_matrix",
            gko::ext::kokkos::make_policy_top<Kokkos::MDRangePolicy, Kokkos::Rank<2>>(Kokkos::Array<int64_t, 2>{0, 0},
                                                                                      Kokkos::Array<int64_t, 2>{static_cast<int64_t>(mtx->get_size()[0]),
                                                                                                                static_cast<int64_t>(mtx->get_size()[1])}),
            [dims, acc] GKO_KOKKOS_FN(const int &i, const int &j, auto k_mtx, auto k_x) {
              if (j < dims - 1) {
                k_mtx(i, j) = acc(i, j, k_x);
              } else {
                k_mtx(i, j) = 1;
              }
            },
            mtx.get(), x.get()));
  };

  if (dynamic_cast<const gko::ReferenceExecutor *>(exec.get()) ||
      dynamic_cast<const gko::OmpExecutor *>(exec.get())) {
    run_kernel(create_row_major_acc());
  } else {
    run_kernel(create_col_major_acc());
  }
}

} // namespace precice::mapping::kernel
