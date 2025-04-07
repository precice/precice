#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
// #include <KokkosBatched_Trsv_TeamVector_Internal.hpp>
// #include <KokkosBatched_Trsv_TeamVector_Impl.hpp>
#include <KokkosBatched_Util.hpp>

// #include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBlas2_gemv.hpp>
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
#ifdef PRECICE_WITH_OPENMP
  if (execName == "omp-executor") {
    return gko::ext::kokkos::create_executor(Kokkos::OpenMP{});
  }
#endif
#ifdef PRECICE_WITH_CUDA
  if (execName == "cuda-executor") {
    if (enableUnifiedMemory) {
      return gko::ext::kokkos::create_executor(Kokkos::Cuda{}, Kokkos::CudaUVMSpace{});
    } else {
      return gko::ext::kokkos::create_executor(Kokkos::Cuda{}, Kokkos::CudaSpace{});
    }
  }
#endif
#ifdef PRECICE_WITH_HIP
  if (execName == "hip-executor") {
    if (enableUnifiedMemory) {
      return gko::ext::kokkos::create_executor(Kokkos::HIP{}, Kokkos::HIPManagedSpace{});
    } else {
      return gko::ext::kokkos::create_executor(Kokkos::HIP{}, Kokkos::HIPSpace{});
    }
  }
#endif
  PRECICE_UNREACHABLE("Executor {} unknown to preCICE", execName);
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
          // contains some host-only code, which yields errors when compiling in
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
                              bool                                      unifiedMemory,
                              gko::ptr_param<GinkgoMatrix>              mtx,
                              const std::array<bool, 3>                 activeAxis,
                              gko::ptr_param<GinkgoMatrix>              supportPoints,
                              gko::ptr_param<GinkgoMatrix>              targetPoints,
                              EvalFunctionType                          f,
                              ::precice::mapping::RadialBasisParameters rbf_params,
                              bool                                      addPolynomial,
                              unsigned int                              extraDims)
{
#ifdef KOKKOS_ENABLE_SERIAL
  if (std::dynamic_pointer_cast<const gko::ReferenceExecutor>(exec)) {
    create_rbf_system_matrix_impl<Kokkos::Serial::memory_space>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
    return;
  }
#endif
#ifdef PRECICE_WITH_OPENMP
  if (std::dynamic_pointer_cast<const gko::OmpExecutor>(exec)) {
    create_rbf_system_matrix_impl<Kokkos::OpenMP::memory_space>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
    return;
  }
#endif
#ifdef PRECICE_WITH_CUDA
  if (std::dynamic_pointer_cast<const gko::CudaExecutor>(exec)) {
    if (unifiedMemory) {
      create_rbf_system_matrix_impl<Kokkos::CudaUVMSpace>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
    } else {
      create_rbf_system_matrix_impl<Kokkos::CudaSpace>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
    }
    return;
  }
#endif
#ifdef PRECICE_WITH_HIP
  if (std::dynamic_pointer_cast<const gko::HipExecutor>(exec)) {
    if (unifiedMemory) {
      create_rbf_system_matrix_impl<Kokkos::HIPManagedSpace>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
    } else {
      create_rbf_system_matrix_impl<Kokkos::HIPSpace>(exec, mtx, activeAxis, supportPoints, targetPoints, f, rbf_params, addPolynomial, extraDims);
    }
    return;
  }
#endif
  PRECICE_UNREACHABLE("Executor unknown to preCICE");
}

#define PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(_function_type)                                                                            \
  template void create_rbf_system_matrix<_function_type>(std::shared_ptr<const gko::Executor> exec,                                             \
                                                         bool                                 unifiedMemory,                                    \
                                                         gko::ptr_param<GinkgoMatrix> mtx, const std::array<bool, 3> activeAxis,                \
                                                         gko::ptr_param<GinkgoMatrix> supportPoints, gko::ptr_param<GinkgoMatrix> targetPoints, \
                                                         _function_type f, RadialBasisParameters rbf_params, bool addPolynomial, unsigned int extraDims)

PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(ThinPlateSplines);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(Multiquadrics);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(InverseMultiquadrics);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(VolumeSplines);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(Gaussian);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(CompactThinPlateSplinesC2);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC0);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC2);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC4);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC6);
PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX(CompactPolynomialC8);
#undef PRECICE_INSTANTIATE_CREATE_RBF_SYSTEM_MATRIX

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
                            bool                                 unifiedMemory,
                            gko::ptr_param<GinkgoMatrix>         mtx,
                            gko::ptr_param<const GinkgoMatrix>   x,
                            const unsigned int                   dims)
{
#ifdef KOKKOS_ENABLE_SERIAL
  if (std::dynamic_pointer_cast<const gko::ReferenceExecutor>(exec)) {
    fill_polynomial_matrix_impl<Kokkos::Serial::memory_space>(exec, mtx, x, dims);
    return;
  }
#endif
#ifdef PRECICE_WITH_OPENMP
  if (std::dynamic_pointer_cast<const gko::OmpExecutor>(exec)) {
    fill_polynomial_matrix_impl<Kokkos::OpenMP::memory_space>(exec, mtx, x, dims);
    return;
  }
#endif
#ifdef PRECICE_WITH_CUDA
  if (auto p = std::dynamic_pointer_cast<const gko::CudaExecutor>(exec); p) {
    if (unifiedMemory) {
      fill_polynomial_matrix_impl<Kokkos::CudaUVMSpace>(exec, mtx, x, dims);
    } else {
      fill_polynomial_matrix_impl<Kokkos::CudaSpace>(exec, mtx, x, dims);
    }
    return;
  }
#endif
#ifdef PRECICE_WITH_HIP
  if (std::dynamic_pointer_cast<const gko::HipExecutor>(exec)) {
    if (unifiedMemory) {
      fill_polynomial_matrix_impl<Kokkos::HIPManagedSpace>(exec, mtx, x, dims);
    } else {
      fill_polynomial_matrix_impl<Kokkos::HIPSpace>(exec, mtx, x, dims);
    }
    return;
  }
#endif
  PRECICE_UNREACHABLE("Executor unknown to preCICE");
}

void compute_offsets(const Kokkos::View<int *> src1, const Kokkos::View<int *> src2,
                     Kokkos::View<std::size_t *> dst, int N)
{
  PRECICE_ASSERT(src1.extent(0) == src2.extent(0));
  PRECICE_ASSERT(src2.extent(0) == dst.extent(0));
  Kokkos::parallel_scan("compute_offsets", N, KOKKOS_LAMBDA(const int i, size_t &update, const bool final) {
    // Number of rows for local system i
    int nrows = src1(i+1) - src1(i);
    // Number of columns for local system i
    int ncols = src2(i+1) - src2(i);

    // Number of entries in the i-th local matrix
    size_t localSize = size_t(nrows) * size_t(ncols);

    // Add to running sum
    update += localSize;

    // 'final == true' indicates we should write to matrixOffsets
    if (final) {
      // matrixOffsets(i+1) = partial sum up to i
      dst(i+1) = update;
    } });
}

template <typename EvalFunctionType, typename MemorySpace>
void do_batched_assembly(
    int                                                              N,   // Number of local systems
    int                                                              dim, // Dimension of points
    EvalFunctionType                                                 f,
    ::precice::mapping::RadialBasisParameters                        rbf_params,
    const Kokkos::View<int *, MemorySpace>                          &inOffsets, // vertex offsets (length N+1)
    const Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> &inCoords,  // meshes
    const Kokkos::View<int *, MemorySpace>                          &targetOffsets,
    const Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> &targetCoords,
    const Kokkos::View<size_t *, MemorySpace>                       &matrixOffsets,
    Kokkos::View<double *, MemorySpace>                              matrices) // 1D view of batched matrices
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMember = typename TeamPolicy::member_type;

  // We launch one team per local system
  Kokkos::parallel_for("do_batched_assembly", TeamPolicy(N, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamMember &team) {
    const int batch = team.league_rank();
    // Ranges
    const int inBegin     = inOffsets(batch);
    const int inEnd       = inOffsets(batch + 1);
    const int targetBegin = targetOffsets(batch);
    const int targetEnd   = targetOffsets(batch + 1);

    // For our batched matrix, this results in
    const int nrows = targetEnd - targetBegin;
    const int ncols = inEnd - inBegin;

    // The matrix offset
    const size_t matrixBegin = matrixOffsets(batch);
    // const size_t matrixEnd   = matrixOffsets(batch + 1);

    // Create an unmanaged 2D subview pointing into matrices
    // This constructor: View(pointer, layout)
    Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        localMatrix(&matrices(matrixBegin), nrows, ncols);

    // Now fill localMatrix(r,c). We'll do a standard 2D nested parallel loop
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nrows),
        [=](int r) {
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(team, ncols),
              [=](int c) {
                // global indices in the original support/target arrays
                int targetIdx = targetBegin + r;
                int inIdx     = inBegin + c;

                // 1) Compute Euclidean distance
                double dist = 0;
                // supportIdx, targetIdx,
                //     inCoords, targetCoords, dim;

                for (int d = 0; d < dim; ++d) {
                  double diff = inCoords(inIdx, d) - targetCoords(targetIdx, d);
                  dist += diff * diff;
                }
                dist = Kokkos::sqrt(dist);

                // 2) Evaluate your RBF or similar function
                double val = f(dist, rbf_params);

                // 3) Store into localMatrix (2D)
                localMatrix(r, c) = val;
              }); // ThreadVectorRange
        });       // TeamThreadRange
  });
}

#define PRECICE_INSTANTIATE(_function_type)                                                             \
  template void do_batched_assembly<_function_type, Kokkos::DefaultExecutionSpace>(                     \
      int                                                                                N,             \
      int                                                                                dim,           \
      _function_type                                                                     f,             \
      ::precice::mapping::RadialBasisParameters                                          rbf_params,    \
      const Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          &inOffsets,     \
      const Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> &inCoords,      \
      const Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          &targetOffsets, \
      const Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> &targetCoords,  \
      const Kokkos::View<size_t *, Kokkos::DefaultExecutionSpace>                       &matrixOffsets, \
      Kokkos::View<double *, Kokkos::DefaultExecutionSpace>                              matrices)

PRECICE_INSTANTIATE(ThinPlateSplines);
PRECICE_INSTANTIATE(Multiquadrics);
PRECICE_INSTANTIATE(InverseMultiquadrics);
PRECICE_INSTANTIATE(VolumeSplines);
PRECICE_INSTANTIATE(Gaussian);
PRECICE_INSTANTIATE(CompactThinPlateSplinesC2);
PRECICE_INSTANTIATE(CompactPolynomialC0);
PRECICE_INSTANTIATE(CompactPolynomialC2);
PRECICE_INSTANTIATE(CompactPolynomialC4);
PRECICE_INSTANTIATE(CompactPolynomialC6);
PRECICE_INSTANTIATE(CompactPolynomialC8);
#undef PRECICE_INSTANTIATE

template <typename MemorySpace>
void do_batched_lu(
    int                                        N,
    const Kokkos::View<size_t *, MemorySpace> &matrixOffsets,
    Kokkos::View<double *, MemorySpace>        matrices)
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;

  Kokkos::parallel_for("do_batched_lu", TeamPolicy(N, Kokkos::AUTO), KOKKOS_LAMBDA(const MemberType &team) {
    const int i = team.league_rank();
    size_t start = matrixOffsets(i);
    size_t end = matrixOffsets(i + 1);
    size_t n = static_cast<size_t>(Kokkos::sqrt(end - start));

    Kokkos::View<double**, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
              A(&matrices(start), n, n);

   KokkosBatched::TeamLU<MemberType,KokkosBatched::Algo::LU::Blocked>::invoke(team,A); });
}

template <typename MemorySpace>
void do_batched_solve(
    int                                        N,
    const Kokkos::View<int *, MemorySpace>    &rhsOffsets,
    Kokkos::View<double *, MemorySpace>        rhs,
    const Kokkos::View<size_t *, MemorySpace> &matrixOffsets,
    const Kokkos::View<double *, MemorySpace> &matrices,
    const Kokkos::View<size_t *, MemorySpace> &evalOffsets,
    const Kokkos::View<double *, MemorySpace> &evalMat,
    const Kokkos::View<int *, MemorySpace>    &outOffsets,
    Kokkos::View<double *, MemorySpace>        out)
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;
  Kokkos::parallel_for("do_batched_solve", TeamPolicy(N, Kokkos::AUTO), KOKKOS_LAMBDA(const MemberType &team) {
    const int i = team.league_rank();

    size_t start = matrixOffsets(i);

    int bStart = rhsOffsets(i);
    int bEnd   = rhsOffsets(i + 1);

    auto n = bEnd - bStart;
    // The lu inplace lu decomposition computed with Kokkosbatched
    Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        A(&matrices(start), n, n);

    // The RHS
    Kokkos::View<double *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        b(&rhs(bStart), n);

    // Forward substitution: solve L * y = b
    // TODO: Check again how we can use TeamVector instead
    // Seems to be available as Unblocked version only
    KokkosBatched::Trsv<
        MemberType,
        KokkosBatched::Uplo::Lower,
        KokkosBatched::Trans::NoTranspose,
        KokkosBatched::Diag::Unit,
        KokkosBatched::Mode::Team,
        KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, A, b);

    team.team_barrier();

    // Backward substitution: solve U * x = y
    KokkosBatched::Trsv<
        MemberType,
        KokkosBatched::Uplo::Upper,
        KokkosBatched::Trans::NoTranspose,
        KokkosBatched::Diag::NonUnit,
        KokkosBatched::Mode::Team,
        KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, A, b);

    team.team_barrier();

    // Next we need the evaluation matrix
    size_t startEval = evalOffsets(i);
    auto   startOut  = outOffsets(i);
    auto   m         = (outOffsets(i + 1) - startOut);

    Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        eval(&evalMat(startEval), m, n);

    Kokkos::View<double *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        result(&out(startOut), m);

KokkosBlas::Experimental::Gemv<
        KokkosBlas::Mode::Team,
        KokkosBlas::Algo::Gemv::Blocked>::invoke(team,'N', 1.0, eval, b, 0.0, result); });
}

} // namespace kernel
} // namespace precice::mapping
