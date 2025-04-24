#include <KokkosBatched_LU_Decl.hpp>
// #include <KokkosBatched_LU_Team_Impl.hpp>
// #include <KokkosBatched_Trsv_TeamVector_Internal.hpp>
// #include <KokkosBatched_Trsv_TeamVector_Impl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>
#include <KokkosBatched_Util.hpp>
// #include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_ApplyPivot_Decl.hpp>
#include <KokkosBatched_ApplyQ_Decl.hpp>
#include <KokkosBatched_QR_WithColumnPivoting_Decl.hpp>
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

bool compute_weights(const std::size_t                                                           nCenters,
                     const std::size_t                                                           nWeights,
                     const std::size_t                                                           nMeshVertices,
                     const int                                                                   dim,
                     Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          offsets,
                     Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> centers,
                     Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          globalIDs,
                     Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> mesh,
                     const CompactPolynomialC2                                                  &w,
                     Kokkos::View<double *, Kokkos::DefaultExecutionSpace>                       normalizedWeights)
{
  using ExecSpace  = typename Kokkos::DefaultExecutionSpace;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMember = typename TeamPolicy::member_type;

  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> weightSum("weightSum", nMeshVertices);
  Kokkos::deep_copy(weightSum, 0.0);
  Kokkos::fence();

  const auto rbf_params = w.getFunctionParameters();

  // We launch one team per local system
  Kokkos::parallel_for("compute_weights", TeamPolicy(nCenters, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamMember &team) {
    const int batch = team.league_rank();
    const int begin = offsets(batch);
    const int end   = offsets(batch + 1);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, begin, end),
        [&](int i) {
          auto   globalID = globalIDs(i);
          double dist     = 0.0;
          for (int d = 0; d < dim; ++d) {
            const double diff = mesh(globalID, d) - centers(batch, d);
            dist += diff * diff;
          }
          dist       = Kokkos::sqrt(dist);
          double val = w(dist, rbf_params);
          // is NUMERICAL_ZERO_DIFFERENCE_DEVICE
          double res           = std::max(val, 1.0e-14);
          normalizedWeights(i) = res;

          Kokkos::atomic_add(&weightSum(globalID), res);
        }); // TeamThreadRange
  });

  // Check for output mesh vertices which are unassigned
  // This check is a pure sanity check
  bool hasZero = false;
  Kokkos::parallel_reduce(
      "check_zero",
      nMeshVertices,
      KOKKOS_LAMBDA(std::size_t i, bool &local) {
        if (weightSum(i) == 0.0)
          local = true;
      },
      Kokkos::LOr<bool>(hasZero));

  if (hasZero) {
    return false;
  }

  // Now scale back the sum
  Kokkos::parallel_for(
      "scale_weights",
      Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, nWeights),
      KOKKOS_LAMBDA(const std::size_t i) {
        const int id = globalIDs(i);
        normalizedWeights(i) /= weightSum(id);
      });

  return true;
}

// TODO: Check if we can specify something meaningful in the launch policy
template <typename MemorySpace>
void do_batched_qr(std::size_t                                               nCluster,
                   int                                                       dim,
                   int                                                       maxClusterSize,
                   Kokkos::View<int *, MemorySpace>                          offsets,
                   Kokkos::View<int *, MemorySpace>                          globalIDs,
                   Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> mesh,
                   Kokkos::View<double *, MemorySpace>                       qrMatrix,
                   Kokkos::View<double *, MemorySpace>                       qrTau,
                   Kokkos::View<int *, MemorySpace>                          qrP)
{
  using TeamPolicy  = Kokkos::TeamPolicy<MemorySpace>;
  using MemberType  = typename TeamPolicy::member_type;
  using ScratchView = Kokkos::View<double *, typename MemorySpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // First, we fill the entire matrix with ones such that we don't need to fill the 1 separately later
  Kokkos::deep_copy(qrMatrix, 1.0);
  Kokkos::fence();

  // Required as workspace for the pivoted QR, see code comment
  // workspace (norm and householder application, 2 * max(m,n) is needed)
  auto scratchSize = ScratchView::shmem_size(2 * maxClusterSize);
  Kokkos::parallel_for("do_batched_qr", TeamPolicy(nCluster, Kokkos::AUTO).set_scratch_size(
                                            /* level = */ 0, Kokkos::PerTeam(scratchSize)),
                       KOKKOS_LAMBDA(const MemberType &team) {
                         // Step 1: define some pointers we need
                         const int batch = team.league_rank();

                         // For the batch
                         const int  begin              = offsets(batch);
                         const auto verticesPerCluster = offsets(batch + 1) - begin; // the local cluster size
                         const int  matrixCols         = dim + 1;                    // our polyParams
                         const int  matrixBegin        = begin * matrixCols;

                         // Step 2: fill the polynomial matrix
                         Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                             qr(&qrMatrix(matrixBegin), verticesPerCluster, matrixCols);

                         Kokkos::parallel_for(
                             Kokkos::TeamThreadRange(team, verticesPerCluster),
                             [&](int i) {
                               auto globalID = globalIDs(i + begin);
                               // the 1 is already set in the last column
                               for (int d = 0; d < dim; ++d) {
                                 qr(i, d) = mesh(globalID, d);
                               }
                             });

                         // Step 3: Compute the QR decomposition
                         const int tauBegin = batch * matrixCols;
                         // the 1 here is for the local rank and has nothing to do with the permutation itself
                         const int PBegin = batch * (matrixCols + 1);

                         Kokkos::View<double *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                             tau(&qrTau(tauBegin), matrixCols);

                         Kokkos::View<int *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                             P(&qrP(PBegin), matrixCols);

                         // Essentially P.end() - 1 =  matrixCols + 1 - 1
                         int &rank = qrP(PBegin + matrixCols);

                         // The scratch memory (shared memory for the device)
                         ScratchView work(team.team_scratch(0), 2 * verticesPerCluster);

                         //  auto b   = Kokkos::subview(work, std::pair<int, int>(0, n));
                         //  auto res = Kokkos::subview(work, std::pair<int, int>(n, n + m));
                         KokkosBatched::TeamVectorQR_WithColumnPivoting<MemberType,
                                                                        KokkosBatched::Algo::QR::Unblocked>::invoke(team, qr, tau,
                                                                                                                    P, work,
                                                                                                                    rank);
                         // We have to define our own criterion for the rank, as the one provided is not stable enough
                         // |pivot|⩽threshold×|maxpivot|
                         // A pivot will be considered nonzero if its absolute value is strictly greater than |pivot|⩽threshold×|maxpivot| where maxpivot is the biggest pivot.
                         double threshold = 1e-6;
                         if (team.team_rank() == 0) {
                           const double maxp = Kokkos::abs(qr(0, 0)); // largest pivot
                           int          r    = 0;
                           for (int i = 0; i < matrixCols; ++i) {
                             if (Kokkos::abs(qr(i, i)) > (threshold * maxp)) {
                               ++r;
                             }
                           }
                           rank = Kokkos::min(r, rank);
                         }
                         // parallel_for
                       });
}

template <typename MemorySpace>
void do_qr_solve(std::size_t                                               nCluster,
                 int                                                       dim,
                 int                                                       maxInClusterSize,
                 Kokkos::View<int *, MemorySpace>                          inOffsets,
                 Kokkos::View<int *, MemorySpace>                          globalInIDs,
                 Kokkos::View<double *, MemorySpace>                       inData,
                 Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> inMesh,
                 Kokkos::View<double *, MemorySpace>                       qrMatrix,
                 Kokkos::View<double *, MemorySpace>                       qrTau,
                 Kokkos::View<int *, MemorySpace>                          qrP,
                 const Kokkos::View<double *, MemorySpace>                 weights,
                 Kokkos::View<int *, MemorySpace>                          outOffsets,
                 Kokkos::View<int *, MemorySpace>                          globalOutIDs,
                 Kokkos::View<double *, MemorySpace>                       outData,
                 Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> outMesh)
{
  using TeamPolicy  = Kokkos::TeamPolicy<MemorySpace>;
  using MemberType  = typename TeamPolicy::member_type;
  using ScratchView = Kokkos::View<double *[2], typename MemorySpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // Used for the inData solution vector and the QR solving
  auto scratchSize = ScratchView::shmem_size(2 * maxInClusterSize);
  Kokkos::parallel_for("do_qr_solve", TeamPolicy(nCluster, Kokkos::AUTO).set_scratch_size(
                                          /* level = */ 0, Kokkos::PerTeam(scratchSize)),
                       KOKKOS_LAMBDA(const MemberType &team) {
                         // Step 1: Some pointers we need
                         const int batch = team.league_rank();
                         // Ranges
                         const int inBegin  = inOffsets(batch);
                         const int inEnd    = inOffsets(batch + 1);
                         const int inSize   = inEnd - inBegin;
                         const int outBegin = outOffsets(batch);
                         const int outEnd   = outOffsets(batch + 1);
                         const int outSize  = outEnd - outBegin;

                         // Step 2: Collect the inData
                         ScratchView tmp(team.team_scratch(0), inSize, 2);
                         auto        in   = Kokkos::subview(tmp, Kokkos::ALL, 0);
                         auto        work = Kokkos::subview(tmp, Kokkos::ALL, 1); // work size is just a guess at the moment

                         Kokkos::parallel_for(
                             Kokkos::TeamThreadRange(team, inSize),
                             [&](int i) {
                               auto globalID = globalInIDs(i + inBegin);
                               in(i)         = inData(globalID);
                             });

                         team.team_barrier();

                         // Step 3: Gather the QR data structures
                         const int matrixCols  = dim + 1;
                         const int matrixBegin = inBegin * matrixCols;
                         const int tauBegin    = batch * matrixCols;
                         const int PBegin      = batch * (matrixCols + 1);
                         const int rank        = qrP(PBegin + matrixCols);

                         Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                             qr(&qrMatrix(matrixBegin), inSize, matrixCols);

                         Kokkos::View<double *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                             tau(&qrTau(tauBegin), matrixCols);

                         Kokkos::View<int *, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                             P(&qrP(PBegin), matrixCols);

                         // Step 4: Solve the linear least-square system A x = b, in our case Q R P^T x = in

                         // Step 4a: Apply Q on the left of in, i.e., y = Q^T * in
                         KokkosBatched::TeamVectorApplyQ<MemberType,
                                                         KokkosBatched::Side::Left,
                                                         KokkosBatched::Trans::Transpose,
                                                         KokkosBatched::Algo::ApplyQ::Unblocked>::invoke(team, qr, tau, in, work);
                         team.team_barrier();

                         auto in_r = Kokkos::subview(in, std::pair<int, int>(0, rank));
                         auto R    = Kokkos::subview(qr, std::pair<int, int>(0, rank), std::pair<int, int>(0, rank));

                         // Step 4b: Solve triangular solve R z = y
                         KokkosBatched::Trsv<
                             MemberType,
                             KokkosBatched::Uplo::Upper,
                             KokkosBatched::Trans::NoTranspose,
                             KokkosBatched::Diag::NonUnit,
                             KokkosBatched::Mode::Team,
                             KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, R, in_r);
                         team.team_barrier();

                         auto res = Kokkos::subview(in, std::pair<int, int>(0, matrixCols));

                         // Steo 4c: zero out the entries which are not within the rank region
                         Kokkos::parallel_for(
                             Kokkos::TeamThreadRange(team, rank, matrixCols), [&](int i) { res(i) = 0; });

                         team.team_barrier();

                         // Step 4d: Apply pivoting x = P z
                         KokkosBatched::TeamVectorApplyPivot<MemberType,
                                                             KokkosBatched::Side::Left,
                                                             KokkosBatched::Direct::Backward>::invoke(team, P, res);
                         team.team_barrier();

                         // Step 5: Subtract polynomial portion from the input data: in -= Q * p
                         // threading over inSize
                         Kokkos::parallel_for(
                             Kokkos::TeamThreadRange(team, inBegin, inEnd),
                             [&](int i) {
                               auto globalID = globalInIDs(i);

                               // The "1"/constant term is the last value in the result
                               double tmp = res(matrixCols - 1);
                               // ... and the linear polynomial
                               for (int d = 0; d < dim; ++d)
                                 tmp += inMesh(globalID, d) * res(d);

                               Kokkos::atomic_sub(&inData(globalID), tmp);
                             });
                         // no barrier needed here

                         // Step 6: Add polynomial portion to the output data: out += V * p
                         // threading over outSize
                         Kokkos::parallel_for(
                             Kokkos::TeamThreadRange(team, outBegin, outEnd),
                             [&](int i) {
                               auto globalID = globalOutIDs(i);

                               // The "1"/constant term is the last value in the result
                               double tmp = res(matrixCols - 1);
                               // ... and the linear polynomial
                               for (int d = 0; d < dim; ++d)
                                 tmp += outMesh(globalID, d) * res(d);

                               // and the weight as usual
                               double w = weights(i);
                               Kokkos::atomic_add(&outData(globalID), tmp * w);
                             });
                         // end parallel_for
                       });
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
    const Kokkos::View<int *, MemorySpace>                          &globalInIDs,
    const Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> &inCoords, // meshes
    const Kokkos::View<int *, MemorySpace>                          &targetOffsets,
    const Kokkos::View<int *, MemorySpace>                          &globalTargetIDs,
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

                auto globalIn     = globalInIDs(inIdx);
                auto globalTarget = globalTargetIDs(targetIdx);
                // 1) Compute Euclidean distance
                double dist = 0;
                for (int d = 0; d < dim; ++d) {
                  double diff = inCoords(globalIn, d) - targetCoords(globalTarget, d);
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

#define PRECICE_INSTANTIATE(_function_type)                                                               \
  template void do_batched_assembly<_function_type, Kokkos::DefaultExecutionSpace>(                       \
      int                                                                                N,               \
      int                                                                                dim,             \
      _function_type                                                                     f,               \
      ::precice::mapping::RadialBasisParameters                                          rbf_params,      \
      const Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          &inOffsets,       \
      const Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          &globalInIDs,     \
      const Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> &inCoords,        \
      const Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          &targetOffsets,   \
      const Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          &globalTargetIDs, \
      const Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> &targetCoords,    \
      const Kokkos::View<size_t *, Kokkos::DefaultExecutionSpace>                       &matrixOffsets,   \
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
    const int i     = team.league_rank();
    size_t    start = matrixOffsets(i);
    size_t    end   = matrixOffsets(i + 1);
    size_t    n     = static_cast<size_t>(Kokkos::sqrt(end - start));

    Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        A(&matrices(start), n, n);

    KokkosBatched::TeamLU<MemberType, KokkosBatched::Algo::LU::Blocked>::invoke(team, A);
    // Parallel end
  });
}

template <typename MemorySpace>
void do_batched_solve(
    std::size_t                                N,
    std::size_t                                maxClusterSize,
    const Kokkos::View<int *, MemorySpace>    &rhsOffsets,
    const Kokkos::View<int *, MemorySpace>    &globalRhsIDs,
    Kokkos::View<double *, MemorySpace>        rhs,
    const Kokkos::View<size_t *, MemorySpace> &matrixOffsets,
    const Kokkos::View<double *, MemorySpace> &matrices,
    const Kokkos::View<double *, MemorySpace> &normalizedWeights,
    const Kokkos::View<size_t *, MemorySpace> &evalOffsets,
    const Kokkos::View<double *, MemorySpace> &evalMat,
    const Kokkos::View<int *, MemorySpace>    &outOffsets,
    const Kokkos::View<int *, MemorySpace>    &globalOutIDs,
    Kokkos::View<double *, MemorySpace>        out)
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;

  using ScratchView = Kokkos::View<double *,
                                   Kokkos::DefaultExecutionSpace::scratch_memory_space,
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // Once for indata and once for outdata
  auto       scratchSize = ScratchView::shmem_size(2 * maxClusterSize);
  TeamPolicy policy(N, Kokkos::AUTO);
  policy.set_scratch_size(
      /* level = */ 0, Kokkos::PerTeam(scratchSize));

  Kokkos::parallel_for("do_batched_solve", policy, KOKKOS_LAMBDA(const MemberType &team) {
    const int batch = team.league_rank();

    size_t start = matrixOffsets(batch);

    // TODO: We could potentially remove the rhsOffsets here and use a sqrt instead
    int  bStart = rhsOffsets(batch);
    int  bEnd   = rhsOffsets(batch + 1);
    auto n      = bEnd - bStart;

    auto startOut = outOffsets(batch);
    auto m        = (outOffsets(batch + 1) - startOut);

    // The lu inplace lu decomposition computed with Kokkosbatched
    Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        A(&matrices(start), n, n);

    // The scratch memory (shared memory for the device)
    ScratchView work(team.team_scratch(0), n + m);
    auto        b   = Kokkos::subview(work, std::pair<int, int>(0, n));
    auto        res = Kokkos::subview(work, std::pair<int, int>(n, n + m));

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, n),
        [&](int i) {
          auto globalID = globalRhsIDs(i + bStart);
          b(i)          = rhs(globalID);
        });

    // Zero out the result (more of a safety feature)
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, m),
        [&](int i) { res(i) = 0; });

    team.team_barrier();

    // Forward substitution: solve L * y = b
    // TODO: Check again how we can use TeamVector instead
    // Seems to be available as Unblocked version only

    // There is also a convenience routine for LU
    // KokkosBatched::SolveLU<MemberType,
    //                        KokkosBatched::Trans::NoTranspose,
    //                        KokkosBatched::Mode::Team,
    //                        KokkosBatched::Algo::SolveLU::Blocked>(team, 1.0, A, b);
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
    size_t startEval = evalOffsets(batch);

    Kokkos::View<double **, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        eval(&evalMat(startEval), m, n);

    // res := 1.0 * eval * b + 0.0 * res
    KokkosBlas::Experimental::Gemv<
        KokkosBlas::Mode::Team,
        KokkosBlas::Algo::Gemv::Blocked>::invoke(team, 'N', 1.0, eval, b, 0.0, res);

    team.team_barrier();

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, m),
        [&](int i) {
          auto w        = normalizedWeights(i + startOut);
          auto globalID = globalOutIDs(i + startOut);
          Kokkos::atomic_add(&out(globalID), res(i) * w);
        }); // TeamThreadRange
    // End Team parallel loop
  });
}

} // namespace kernel
} // namespace precice::mapping
