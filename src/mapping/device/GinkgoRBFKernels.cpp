#include "mapping/device/GinkgoRBFKernels.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "math/math.hpp"

#include <functional>
#include <ginkgo/extensions/kokkos.hpp>
#include <ginkgo/ginkgo.hpp>

using precice::mapping::RadialBasisParameters;
using precice::math::pow_int;

using gko::ext::kokkos::map_data;

namespace precice::mapping {

std::shared_ptr<gko::Executor> create_device_executor(const std::string &execName, bool enableUnifiedMemory)
{

#ifdef KOKKOS_ENABLE_SERIAL
  if (execName == "reference-executor") {
    return gko::ext::kokkos::create_executor(Kokkos::Serial{});
  }
#endif
  if (enableUnifiedMemory) {
    return gko::ext::kokkos::create_executor(Kokkos::DefaultExecutionSpace{},
                                             Kokkos::SharedSpace{});
  } else {
    return gko::ext::kokkos::create_executor(Kokkos::DefaultExecutionSpace{});
  }
}

namespace kernel {

template <typename MemorySpace, typename EvalFunctionType>
void create_rbf_system_matrix_impl(std::shared_ptr<const gko::Executor>      exec,
                                   gko::ptr_param<GinkgoMatrix>              mtx,
                                   const std::array<bool, 3>                 activeAxis,
                                   gko::ptr_param<GinkgoMatrix>              supportPoints,
                                   gko::ptr_param<GinkgoMatrix>              targetPoints,
                                   EvalFunctionType                          f,
                                   ::precice::mapping::RadialBasisParameters rbf_params,
                                   bool                                      addPolynomial,
                                   unsigned int                              extraDims)
{
  auto k_mtx           = map_data<MemorySpace>(mtx.get());
  auto k_supportPoints = map_data<MemorySpace>(supportPoints.get());
  auto k_targetPoints  = map_data<MemorySpace>(targetPoints.get());

  // the outcome of this if-condition is known at compile time. However, Cuda
  // has issues with constexpr
  if (std::is_same_v<MemorySpace, Kokkos::HostSpace>) {
    // Row-major access
    Kokkos::parallel_for(
        "create_rbf_system_matrix_row_major",
        Kokkos::MDRangePolicy<typename MemorySpace::execution_space, Kokkos::Rank<2>>{{0, 0}, {mtx->get_size()[0], mtx->get_size()[1]}},
        KOKKOS_LAMBDA(const int &i, const int &j) {
          double dist = 0;
          if (addPolynomial) {
            k_mtx(i, j) = 0; // Zero the matrix entry if polynomial terms are added
          }
          // We need to use a pointer here, because the bound checking of std::array
          // contins some host-only code, which yields errors when compiling in
          // debug mode
          const bool *deviceActiveAxis = activeAxis.data();

          // Compute Euclidean distance using row-major indexing
          for (size_t k = 0; k < activeAxis.size(); ++k) {
            if (deviceActiveAxis[k]) {
              double diff = k_supportPoints(j, k) - k_targetPoints(i, k);
              dist += diff * diff;
            }
          }
          dist        = Kokkos::sqrt(dist);
          k_mtx(i, j) = f(dist, rbf_params); // Evaluate the RBF function
        });
  } else {
    // Column-major access
    Kokkos::parallel_for(
        "create_rbf_system_matrix_col_major",
        Kokkos::MDRangePolicy<typename MemorySpace::execution_space, Kokkos::Rank<2>>{{0, 0}, {mtx->get_size()[0], mtx->get_size()[1]}},
        KOKKOS_LAMBDA(const int &i, const int &j) {
          double dist = 0;
          if (addPolynomial) {
            k_mtx(i, j) = 0; // Zero the matrix entry if polynomial terms are added
          }
          const bool *deviceActiveAxis = activeAxis.data();

          // Compute Euclidean distance using column-major indexing
          for (size_t k = 0; k < activeAxis.size(); ++k) {
            if (deviceActiveAxis[k]) {
              double diff = k_supportPoints(k, j) - k_targetPoints(k, i);
              dist += diff * diff;
            }
          }
          dist        = Kokkos::sqrt(dist);
          k_mtx(i, j) = f(dist, rbf_params); // Evaluate the RBF function
        });
  }
}

template <typename EvalFunctionType>
void create_rbf_system_matrix(std::shared_ptr<const gko::Executor>      exec,
                              gko::ptr_param<GinkgoMatrix>              mtx,
                              const std::array<bool, 3>                 activeAxis,
                              gko::ptr_param<GinkgoMatrix>              supportPoints,
                              gko::ptr_param<GinkgoMatrix>              targetPoints,
                              EvalFunctionType                          f,
                              ::precice::mapping::RadialBasisParameters rbf_params,
                              bool                                      addPolynomial,
                              unsigned int                              extraDims)
{
  if (std::dynamic_pointer_cast<const gko::ReferenceExecutor>(exec)) {
    create_rbf_system_matrix_impl<Kokkos::HostSpace>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
  } else {
    create_rbf_system_matrix_impl<Kokkos::DefaultExecutionSpace::memory_space>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
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
INSTATIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC8);

template <typename MemorySpace>
void fill_polynomial_matrix_impl(std::shared_ptr<const gko::Executor> exec,
                                 gko::ptr_param<GinkgoMatrix>         mtx,
                                 gko::ptr_param<const GinkgoMatrix>   x,
                                 const unsigned int                   dims)
{

  auto k_mtx = map_data<MemorySpace>(mtx.get());
  auto k_x   = map_data<MemorySpace>(x.get());

  if (std::is_same_v<MemorySpace, Kokkos::HostSpace>) {
    Kokkos::parallel_for(
        "fill_polynomial_matrix_row_major",
        Kokkos::MDRangePolicy<typename MemorySpace::execution_space, Kokkos::Rank<2>>{{0, 0}, {mtx->get_size()[0], mtx->get_size()[1]}},
        KOKKOS_LAMBDA(const int &i, const int &j) {
          double value;
          if (j < dims - 1) {
            value = k_x(i, j); // Row-major access
          } else {
            value = 1;
          }
          k_mtx(i, j) = value;
        });
  } else {
    Kokkos::parallel_for(
        "fill_polynomial_matrix_col_major",
        Kokkos::MDRangePolicy<typename MemorySpace::execution_space, Kokkos::Rank<2>>{{0, 0}, {mtx->get_size()[0], mtx->get_size()[1]}},
        KOKKOS_LAMBDA(const int &i, const int &j) {
          double value;
          if (j < dims - 1) {
            value = k_x(j, i); // Column-major access
          } else {
            value = 1;
          }
          k_mtx(i, j) = value;
        });
  }
}

void fill_polynomial_matrix(std::shared_ptr<const gko::Executor> exec,
                            gko::ptr_param<GinkgoMatrix>         mtx,
                            gko::ptr_param<const GinkgoMatrix>   x,
                            const unsigned int                   dims)
{
  if (std::dynamic_pointer_cast<const gko::ReferenceExecutor>(exec)) {
    fill_polynomial_matrix_impl<Kokkos::HostSpace>(exec, mtx, x, dims);
  } else {
    fill_polynomial_matrix_impl<Kokkos::DefaultExecutionSpace::memory_space>(exec, mtx, x, dims);
  }
}
} // namespace kernel
} // namespace precice::mapping
