#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_Util.hpp>

#include <KokkosBatched_ApplyPivot_Decl.hpp>
#include <KokkosBatched_ApplyQ_Decl.hpp>
#include <KokkosBatched_QR_WithColumnPivoting_Decl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBlas2_gemv.hpp>

#include "mapping/impl/BasisFunctions.hpp"
#include "math/math.hpp"

#include <functional>

using precice::mapping::RadialBasisParameters;
using precice::math::pow_int;

namespace precice::mapping::kernel {

namespace impl {
// Helper to compute the closest power of 2^x for the given n
// used only below to determine a reasonable team size
inline int nearestPowerOfTwo(int n)
{
  if (n < 1)
    return 1;
  // find highest power-of-two <= n
  int pow2 = 1;
  while (pow2 <= n) {
    pow2 = pow2 * 2;
  }
  // now that we have the value larger, we look for the pow2 below n
  int lower = pow2 / 2;
  // decide whether n is closer to lower or upper (pow2)
  if (n - lower < pow2 - n) {
    return lower;
  } else {
    return pow2;
  }
}

// Try to find a good team size based on the cluster size
// avgWork is in our case simply the average cluster size, as it determines
// the local matrix sizes
template <typename ExecSpace, typename FunctorType, typename Policy>
auto findTeamSize(int avgWork, const FunctorType &functor, const Policy &policy)
{
  // If using OpenMP, Kokkos::AUTO works best
  if constexpr (std::is_same_v<ExecSpace, Kokkos::HostSpace::execution_space>) {
    return Kokkos::AUTO;
  } else {
    // Compute the closest power of two for the work we have
    int teamSize = nearestPowerOfTwo(avgWork);
    // Ensure minimum of one warp for the GPU to avoid partial-warp inefficiency
    int warpSize = Kokkos::TeamPolicy<ExecSpace>::vector_length_max();
    if (teamSize < warpSize) {
      teamSize = warpSize;
    }
    // Query Kokkos for the max recommended team size for this functor
    // (takes also memory constraints into account)
    int maxRecommended = policy.team_size_recommended(functor, Kokkos::ParallelForTag());
    if (teamSize > maxRecommended) {
      teamSize = maxRecommended;
    }
    return teamSize;
  }
}

// small helper function to make the compiler handle variables in lambdas which are only conditionally used
template <typename... Args>
KOKKOS_INLINE_FUNCTION constexpr void capture_conditional_variables(const Args &...) {}
} // namespace impl

// For within the kernels
// The batch matrix has a default layout, which is consistent across kernels
template <typename MemorySpace = ExecutionSpace>
using BatchMatrix = Kokkos::View<double **, MemorySpace, UnmanagedMemory>;
template <typename T = double *, typename MemorySpace = ExecutionSpace>
using BatchVector = Kokkos::View<T, MemorySpace, UnmanagedMemory>;

template <typename MemorySpace>
bool compute_weights(const int                     nCenters,
                     const offset_1d_type          nWeights,
                     const int                     nMeshVertices,
                     const int                     dim,
                     VectorOffsetView<MemorySpace> offsets,
                     MeshView<MemorySpace>         centers,
                     GlobalIDView<MemorySpace>     globalIDs,
                     MeshView<MemorySpace>         mesh,
                     const CompactPolynomialC2    &w,
                     VectorView<MemorySpace>       normalizedWeights)
{
  using TeamPolicy = Kokkos::TeamPolicy<MemorySpace>;
  using TeamMember = typename TeamPolicy::member_type;

  VectorView<MemorySpace> weightSum("weightSum", nMeshVertices);
  Kokkos::deep_copy(weightSum, 0.0);
  Kokkos::fence();

  const auto rbf_params = w.getFunctionParameters();

  // We launch one team per local system
  Kokkos::parallel_for("compute_weights", TeamPolicy(nCenters, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamMember &team) {
    const int  batch = team.league_rank();
    const auto begin = offsets(batch);
    const auto end   = offsets(batch + 1);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, begin, end),
        [&](offset_1d_type i) {
          auto   globalID = globalIDs(i);
          double dist     = 0.0;
          for (int d = 0; d < dim; ++d) {
            double diff = mesh(globalID, d) - centers(batch, d);
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
      KOKKOS_LAMBDA(int i, bool &local) {
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
      Kokkos::RangePolicy<MemorySpace>(0, nWeights),
      KOKKOS_LAMBDA(const int i) {
        const int id = globalIDs(i);
        normalizedWeights(i) /= weightSum(id);
      });

  return true;
}

// TODO: Check if we can specify something meaningful in the launch policy
template <typename MemorySpace>
void do_batched_qr(int                           nCluster,
                   int                           dim,
                   int                           avgClusterSize,
                   int                           maxClusterSize,
                   VectorOffsetView<MemorySpace> offsets,
                   GlobalIDView<MemorySpace>     globalIDs,
                   MeshView<MemorySpace>         mesh,
                   VectorView<MemorySpace>       qrMatrix,
                   VectorView<MemorySpace>       qrTau,
                   PivotView<MemorySpace>        qrP)
{
  using TeamPolicy  = Kokkos::TeamPolicy<MemorySpace>;
  using MemberType  = typename TeamPolicy::member_type;
  using ScratchView = Kokkos::View<double *, typename MemorySpace::scratch_memory_space, UnmanagedMemory>;

  // First, we fill the entire matrix with ones such that we don't need to fill the 1 separately later
  Kokkos::deep_copy(qrMatrix, 1.0);
  Kokkos::fence();

  auto kernel = KOKKOS_LAMBDA(const MemberType &team)
  {
    // Step 1: define some pointers we need
    const int batch = team.league_rank();

    // For the batch
    const auto           begin              = offsets(batch);
    const int            verticesPerCluster = offsets(batch + 1) - begin; // the local cluster size
    const int            matrixCols         = dim + 1;                    // our polyParams
    const offset_1d_type matrixBegin        = begin * matrixCols;

    // Step 2: fill the polynomial matrix
    BatchMatrix<MemorySpace> qr(&qrMatrix(matrixBegin), verticesPerCluster, matrixCols);

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
    const offset_1d_type tauBegin = batch * matrixCols;
    // the 1 here is for the local rank and has nothing to do with the permutation itself
    const offset_1d_type PBegin = batch * (matrixCols + 1);

    BatchVector<double *, MemorySpace> tau(&qrTau(tauBegin), matrixCols);
    BatchVector<int *, MemorySpace>    P(&qrP(PBegin), matrixCols);

    // Essentially P.end() - 1 =  matrixCols + 1 - 1
    int &rank = qrP(PBegin + matrixCols);

    // The scratch memory (shared memory for the device)
    ScratchView work(team.team_scratch(0), 2 * verticesPerCluster);

    //  auto b   = Kokkos::subview(work, std::pair<int, int>(0, n));
    //  auto res = Kokkos::subview(work, std::pair<int, int>(n, n + m));
    // A Serial QR_WithColumnPivoting would probably be better suited here, but this is as of now not available
    // The workload to put a Team on the whole matrix is most likely not a good balance
    KokkosBatched::TeamVectorQR_WithColumnPivoting<MemberType,
                                                   KokkosBatched::Algo::QR::Unblocked>::invoke(team, qr, tau,
                                                                                               P, work,
                                                                                               rank);
    // We have to define our own criterion for the rank, as the one provided is not stable enough
    // A pivot will be considered nonzero if its absolute value is strictly greater than |pivot|⩽threshold×|maxpivot| where maxpivot is the biggest pivot.
    // We use 1e-6 as threshold, which is much stricter than the Kokkos internal criterion
    // The kokkos internal criterion can be found in KokkosBatched_QR_WithColumnPivoting_TeamVector_Internal.hpp
    // if diagonal value is smaller than threshold(10 * max_diag * ats::epsilon())
    // note that the pivoted algorithm aborts once the rank deficiency is detected, thus, values larger than rank just contain rubbish
    double threshold = 1e-5;
    if (team.team_rank() == 0) {
      const double maxp = Kokkos::abs(qr(0, 0)); // largest pivot
      int          r    = 0;
      for (int i = 0; i < rank; ++i) {
        r += static_cast<int>(Kokkos::abs(qr(i, i)) > (threshold * maxp));
      }
      rank = r;
    }
    // parallel_for
  };

  // Required as workspace for the pivoted QR, see code comment
  // workspace (norm and householder application, 2 * max(m,n) is needed)

  // This configuration gave by far the best performance, using Kokkos::AUTO results in 70ms, whereas
  // this configuration gave about 20ms on the A100
  auto scratchSize = ScratchView::shmem_size(2 * maxClusterSize);
  // Rather than selecting the clustersize, we the time size is here more proportional to 4, since it is the maximum rank.
  auto teamSize = impl::findTeamSize<typename MemorySpace::execution_space>(4, kernel, TeamPolicy(nCluster, Kokkos::AUTO).set_scratch_size(0, Kokkos::PerTeam(scratchSize)));
  // const int VL = TeamPolicy::vector_length_max();
  // The inner loop uses vector parallelism, so we should try to configure accordingly
  std::unique_ptr<TeamPolicy> policy;
  if constexpr (std::is_same_v<decltype(teamSize), Kokkos::AUTO_t>) {
    policy = std::make_unique<TeamPolicy>(nCluster, Kokkos::AUTO);
  } else {
    policy = std::make_unique<TeamPolicy>(nCluster, 4, teamSize / 4);
  }
  policy->set_scratch_size(/* level = */ 0, Kokkos::PerTeam(scratchSize));
  Kokkos::parallel_for("do_batched_qr", *policy, kernel);
}

template <typename MemorySpace>
void compute_offsets(const VectorOffsetView<MemorySpace> src1, const VectorOffsetView<MemorySpace> src2,
                     MatrixOffsetView<MemorySpace> dst, int nCluster)
{
  PRECICE_ASSERT(src1.extent(0) == src2.extent(0));
  PRECICE_ASSERT(src2.extent(0) == dst.extent(0));
  Kokkos::parallel_scan("compute_offsets", nCluster, KOKKOS_LAMBDA(const int i, offset_2d_type &update, const bool final) {
    // Number of rows for local system i
    int nrows = src1(i + 1) - src1(i);
    // Number of columns for local system i
    int ncols = src2(i + 1) - src2(i);

    // Number of entries in the i-th local matrix
    int localSize = nrows * ncols;

    // Add to running sum
    update += static_cast<offset_2d_type>(localSize);

    // 'final == true' indicates we should write to matrixOffsets
    if (final) {
      // matrixOffsets(i+1) = partial sum up to i
      dst(i + 1) = update;
    }
    // end parallel_for
  });
}

template <typename EvalFunctionType, typename MemorySpace>
void do_input_assembly(
    int                                  nCluster, // Number of local systems
    int                                  dim,      // Dimension of points
    int                                  avgClusterSize,
    int                                  maxInClusterSize,
    EvalFunctionType                     f,
    const VectorOffsetView<MemorySpace> &inOffsets, // vertex offsets (length N+1)
    const GlobalIDView<MemorySpace>     &globalInIDs,
    const MeshView<MemorySpace>         &inCoords, // meshes
    const MatrixOffsetView<MemorySpace> &matrixOffsets,
    VectorView<MemorySpace>              matrices) // 1D view of batched matrices
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;

  using ScratchSpace  = typename MemorySpace::scratch_memory_space;
  using ScratchView1d = Kokkos::View<double *, ScratchSpace, UnmanagedMemory>;
  using ScratchMatrix = Kokkos::View<double **, Kokkos::LayoutRight, ScratchSpace, UnmanagedMemory>;

  const auto rbf_params = f.getFunctionParameters();
  auto       kernel     = KOKKOS_LAMBDA(const MemberType &team)
  {
    const int batch = team.league_rank();
    // Ranges
    const auto inBegin = inOffsets(batch);
    const auto inEnd   = inOffsets(batch + 1);
    const int  n       = inEnd - inBegin;

    ScratchMatrix mesh(team.team_scratch(0), n, dim);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, n), [&](int i) {
          auto globalID = globalInIDs(i + inBegin);
          for (int d = 0; d < dim; ++d) {
            mesh(i, d) = inCoords(globalID, d);
          }
        });

    // The matrix offset
    const auto matrixBegin = matrixOffsets(batch);

    team.team_barrier();

    // Create an unmanaged 2D subview pointing into matrices
    // This constructor: View(pointer, layout)
    BatchMatrix<MemorySpace> localMatrix(&matrices(matrixBegin), n, n);

    // Now fill localMatrix(r,c). We'll do a standard 2D nested parallel loop
    Kokkos::parallel_for(
        Kokkos::TeamThreadMDRange(team, n, n),
        [=](int r, int c) {
          // global indices in the original support/target arrays
          // 1) Compute Euclidean distance
          double dist = 0;
          for (int d = 0; d < dim; ++d) {
            double diff = mesh(r, d) - mesh(c, d);
            dist += diff * diff;
          }
          dist = Kokkos::sqrt(dist);

          // 2) Evaluate your RBF or similar function
          double val = f(dist, rbf_params);

          // 3) Store into localMatrix (2D)
          localMatrix(c, r) = val;
        }); // TeamThreadRange
  };

  // We put the solution and the in data values into shared memory
  auto inBytes = ScratchView1d::shmem_size(maxInClusterSize);

  // We could also select here manually 64 for the average of 50 for example
  auto teamSize = impl::findTeamSize<ExecSpace>(avgClusterSize, kernel, TeamPolicy(nCluster, Kokkos::AUTO).set_scratch_size(/* level = */ 0, Kokkos::PerTeam(dim * inBytes)));
  Kokkos::parallel_for("do_input_assembly", TeamPolicy(nCluster, teamSize).set_scratch_size(/* level = */ 0, Kokkos::PerTeam(dim * inBytes)), kernel);
}

// TODO: Using Kokkos::LayoutRight for the Coords performs a bit better for the assembly,
// but it might deteriorate performance related to the polynomials. Especially for the gemv
// for the polynomial contributions etc
// For the full GPU porting, the Layout can only be LayoutRight, as we don't access the coordinates
// coalesced, but rather first pick what we need
template <typename EvalFunctionType, typename MemorySpace>
void do_batched_assembly(
    int                                  nCluster, // Number of local systems
    int                                  dim,      // Dimension of points
    int                                  avgClusterSize,
    EvalFunctionType                     f,
    const VectorOffsetView<MemorySpace> &inOffsets, // vertex offsets (length N+1)
    const GlobalIDView<MemorySpace>     &globalInIDs,
    const MeshView<MemorySpace>         &inCoords, // meshes
    const VectorOffsetView<MemorySpace> &targetOffsets,
    const GlobalIDView<MemorySpace>     &globalTargetIDs,
    const MeshView<MemorySpace>         &targetCoords,
    const MatrixOffsetView<MemorySpace> &matrixOffsets,
    VectorView<MemorySpace>              matrices) // 1D view of batched matrices
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;

  const auto rbf_params = f.getFunctionParameters();
  auto       kernel     = KOKKOS_LAMBDA(const MemberType &team)
  {
    const int batch = team.league_rank();
    // Ranges
    const auto inBegin     = inOffsets(batch);
    const auto inEnd       = inOffsets(batch + 1);
    const auto targetBegin = targetOffsets(batch);
    const auto targetEnd   = targetOffsets(batch + 1);

    // For our batched matrix, this results in
    const int nrows = targetEnd - targetBegin;
    const int ncols = inEnd - inBegin;

    // The matrix offset
    const auto matrixBegin = matrixOffsets(batch);
    // const size_t matrixEnd   = matrixOffsets(batch + 1);

    // Create an unmanaged 2D subview pointing into matrices
    // This constructor: View(pointer, layout)
    BatchMatrix<MemorySpace> localMatrix(&matrices(matrixBegin), nrows, ncols);

    // Now fill localMatrix(r,c). We'll do a standard 2D nested parallel loop
    Kokkos::parallel_for(
        Kokkos::TeamThreadMDRange(team, nrows, ncols),
        [=](int r, int c) {
          // global indices in the original support/target arrays
          offset_1d_type targetIdx = targetBegin + r;
          offset_1d_type inIdx     = inBegin + c;

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
  };

  auto teamSize = impl::findTeamSize<ExecSpace>(avgClusterSize, kernel, TeamPolicy(nCluster, Kokkos::AUTO));
  Kokkos::parallel_for("do_batched_assembly", TeamPolicy(nCluster, teamSize), kernel);
}

template <typename MemorySpace>
void do_batched_lu(
    int                                  nCluster,
    int                                  avgClusterSize,
    const MatrixOffsetView<MemorySpace> &matrixOffsets,
    VectorView<MemorySpace>              matrices)
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;

  auto kernel = KOKKOS_LAMBDA(const MemberType &team)
  {
    const int i     = team.league_rank();
    auto      start = matrixOffsets(i);
    auto      end   = matrixOffsets(i + 1);
    int       n     = static_cast<int>(Kokkos::sqrt(end - start));

    BatchMatrix<MemorySpace> A(&matrices(start), n, n);

    KokkosBatched::TeamLU<MemberType, KokkosBatched::Algo::LU::Blocked>::invoke(team, A);
    // Parallel end
  };

  // Using Kokkos::AUTO resulted in the best performance for all cases
  // auto teamSize = impl::findTeamSize<ExecSpace>(avgClusterSize, kernel, TeamPolicy(nCluster, Kokkos::AUTO));
  Kokkos::parallel_for("do_batched_lu", TeamPolicy(nCluster, Kokkos::AUTO), kernel);
}

/// polynomial: bool whether to evaluate the polynomial or not
/// evaluation_op_available: bool whether to compute the evaluation on the fly (merged)
/// or whether these entries are already precomputed, i.e., evalMat may be a dummy argument
/// in case evaluation_op_available = false
template <bool polynomial, bool evaluation_op_available, typename EvalFunctionType, typename MemorySpace>
void do_batched_solve(
    int                                  nCluster,
    int                                  dim,
    int                                  avgInClusterSize,
    int                                  maxInClusterSize,
    int                                  maxOutClusterSize,
    EvalFunctionType                     f,
    const VectorOffsetView<MemorySpace> &rhsOffsets,
    const GlobalIDView<MemorySpace>     &globalRhsIDs,
    VectorView<MemorySpace>              rhs,
    const MatrixOffsetView<MemorySpace> &matrixOffsets,
    const VectorView<MemorySpace>       &matrices,
    const VectorView<MemorySpace>       &normalizedWeights,
    const MatrixOffsetView<MemorySpace> &evalOffsets,
    const VectorView<MemorySpace>       &evalMat,
    const VectorOffsetView<MemorySpace> &outOffsets,
    const GlobalIDView<MemorySpace>     &globalOutIDs,
    VectorView<MemorySpace>              out,
    // For the polynomial required in addition
    const MeshView<MemorySpace>   &inMesh,
    const MeshView<MemorySpace>   &outMesh,
    const VectorView<MemorySpace> &qrMatrix,
    const VectorView<MemorySpace> &qrTau,
    const PivotView<MemorySpace>  &qrP)
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;

  using ScratchSpace = typename MemorySpace::scratch_memory_space;
  // Layout is important for how we use these matrices: we need to ensure that cols are contiguous in memory
  // We use the scratch memory for the input data and potentially for workspace required for the
  // polynomial or for the input mesh in case we have to compute the evaluation on the fly
  using ScratchView1d = Kokkos::View<double *[1], Kokkos::LayoutLeft, ScratchSpace, UnmanagedMemory>;
  using ScratchView4d = Kokkos::View<double *[4], Kokkos::LayoutLeft, ScratchSpace, UnmanagedMemory>;
  using ScratchVector = Kokkos::View<double *, ScratchSpace, UnmanagedMemory>;
  using ScratchMatrix = std::conditional_t<!evaluation_op_available || polynomial, ScratchView4d, ScratchView1d>;
  using ScratchMesh   = Kokkos::View<double **, Kokkos::LayoutRight, ScratchSpace, UnmanagedMemory>;

  // TODO: Add checks for the matrices (extent() > 0 with the template params)
  // PRECICE_ASSERT()

  const auto rbf_params = f.getFunctionParameters();
  // We define the lambda here such that we can query the recommended team size from Kokkos
  // Launch policy is then handled below
  auto kernel = KOKKOS_LAMBDA(const MemberType &team)
  {
    // Required for correct capturing (mostly by device compilers), as these variables are only conditionally used further down
    impl::capture_conditional_variables(dim, qrMatrix, qrTau, qrP, inMesh, outMesh, evalOffsets, evalMat, globalOutIDs, f, rbf_params, normalizedWeights, out);

    // Step 1: Define some pointers
    // TODO: We could potentially remove the rhsOffsets here and use a sqrt instead
    const int  batch    = team.league_rank();
    const auto inBegin  = rhsOffsets(batch);
    const int  inSize   = rhsOffsets(batch + 1) - inBegin;
    const auto outBegin = outOffsets(batch);
    const int  outSize  = outOffsets(batch + 1) - outBegin;

    // Step 2: Allocate shared memory for the team and fill it with the inData of this cluster
    // The scratch memory (shared memory for the device)
    ScratchMatrix work(team.team_scratch(0), Kokkos::max(4, inSize));
    auto          in = Kokkos::subview(work, std::pair<int, int>(0, inSize), 0);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, inSize), [&](int i) {
          auto globalID = globalRhsIDs(i + inBegin);
          in(i)         = rhs(globalID);
        });
    team.team_barrier();

    Kokkos::Array<double, 4> qrCoeffs = {0., 0., 0., 0.};

    // Step 3: Solve the polynomial QR system, if we have one
    if constexpr (polynomial) {

      // Step 3a: Backup the current in data, since we solve the QR in place
      // In principle, we need a vector here (just as in), but the ApplyQ routine expects a Rank2 matrix,
      // so we have to stick to this particular syntax (keeping it rank 2 with one column)
      auto in_cp = Kokkos::subview(work, std::pair<int, int>(0, inSize), std::pair<int, int>(1, 2));
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, inSize), [&](int i) { in_cp(i, 0) = in(i); });
      team.team_barrier();

      // Step 3b: Define pointers and matrices
      const int            matrixCols = dim + 1;
      const offset_1d_type qrBegin    = inBegin * matrixCols;
      const offset_1d_type tauBegin   = batch * matrixCols;
      const offset_1d_type PBegin     = batch * (matrixCols + 1);
      const int            rank       = qrP(PBegin + matrixCols);

      BatchMatrix<MemorySpace>           qr(&qrMatrix(qrBegin), inSize, matrixCols);
      BatchVector<double *, MemorySpace> tau(&qrTau(tauBegin), matrixCols);
      BatchVector<int *, MemorySpace>    P(&qrP(PBegin), matrixCols);

      // Step 3c: Apply Q on the left of in, i.e., y = Q^T * in
      if (team.team_rank() == 0) {

        // tmp size might be insufficient: there was no size requirement specified for the workspace
        // however, it needs to be contiguous
        auto tmp = Kokkos::subview(work, Kokkos::ALL, 2);
        KokkosBatched::ApplyQ<MemberType,
                              KokkosBatched::Side::Left,
                              KokkosBatched::Trans::Transpose,
                              KokkosBatched::Mode::Serial,
                              KokkosBatched::Algo::ApplyQ::Unblocked>::invoke(team, qr, tau, in_cp, tmp);

        // Step 3d: Solve triangular solve R z = y
        auto in_r = Kokkos::subview(in_cp, std::pair<int, int>(0, rank), 0);
        auto R    = Kokkos::subview(qr, std::pair<int, int>(0, rank), std::pair<int, int>(0, rank));

        KokkosBatched::Trsv<
            MemberType,
            KokkosBatched::Uplo::Upper,
            KokkosBatched::Trans::NoTranspose,
            KokkosBatched::Diag::NonUnit,
            KokkosBatched::Mode::Serial,
            KokkosBatched::Algo::Trsv::Unblocked>::invoke(team, 1.0, R, in_r);
      }

      team.team_barrier();

      // Step 3e: Copy the result over into memory for each thread
      for (int r = 0; r < rank; ++r) {
        qrCoeffs[r] = in_cp(r, 0);
      }

      // Step 3f: Apply pivoting x = P z
      // There is also convenience routines for the pivoting, but we let every thread just
      // apply the pivoting on its own, as it is more compact and the routine doesn't allow
      // using an Array The below is the equivalent of the following (but qrCoeffs would need to be a view):
      // KokkosBatched::TeamVectorApplyPivot<MemberType, Side::Left, Direct::Backward>::invoke(team, P, qrCoeffs);
      for (int i = (matrixCols - 1); i >= 0; --i) {
        Kokkos::kokkos_swap(qrCoeffs[i], qrCoeffs[i + P(i)]);
      }
    }

    // Step 3g: Subtract polynomial portion from the input data: in -= Q * p
    // In case we have to compute the evaluation operator further down, we
    // also pull the in mesh into local shared memory
    // threading over inSize
    // This vector uses memory we have in principle managed by "work" (for the QR solve as tmp),
    // but its layout is different (LayoutRight to have the coordinates aligned in memory).
    // We have to declare it here outsie the "if" to be able to use it further down then.
    // The memory here might point to null (or rather the end of "work"), in case polynomial = false
    // and the evaluation_op is available but then it also remains unused
    ScratchMesh localInMesh(&work(0, 1), inSize, dim);

    if constexpr (!evaluation_op_available || polynomial) {
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, inSize),
          [&](int i) {
            auto globalID = globalRhsIDs(i + inBegin);
            // The "1"/constant term is the last value in the result
            // dim is here matrixCols - 1
            double sum = qrCoeffs[dim];
            // ... and the linear polynomial
            for (int d = 0; d < dim; ++d) {
              sum += inMesh(globalID, d) * qrCoeffs[d];
              // Put it in shared memory as we later need it in the output evaluation
              localInMesh(i, d) = inMesh(globalID, d);
            }
            in(i) -= sum;
          });
      team.team_barrier();
    }

    // Step 4: Solve the LU decomposition
    // The lu inplace lu decomposition computed with KokkosBatched
    // There is also a convenience routine for the LU solve, but it
    // uses Trsm under the hood, which is quite a bit slower
    auto                     matStart = matrixOffsets(batch);
    BatchMatrix<MemorySpace> A(&matrices(matStart), inSize, inSize);

    // Forward substitution: solve L * y = b and
    KokkosBatched::Trsv<
        MemberType,
        KokkosBatched::Uplo::Lower,
        KokkosBatched::Trans::NoTranspose,
        KokkosBatched::Diag::Unit,
        KokkosBatched::Mode::Team,
        KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, A, in);

    team.team_barrier();

    // Backward substitution: solve U * x = y
    KokkosBatched::Trsv<
        MemberType,
        KokkosBatched::Uplo::Upper,
        KokkosBatched::Trans::NoTranspose,
        KokkosBatched::Diag::NonUnit,
        KokkosBatched::Mode::Team,
        KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, A, in);
    team.team_barrier();

    // Step 5: Apply the output operator
    // If we have the evaluation operator, we use it (might be slower though)
    if constexpr (evaluation_op_available) {
      // Step 5a: Allocate and zero out a local result vector (more of a safety feature)
      ScratchVector res(team.team_scratch(1), outSize);
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, outSize),
          [&](int i) { res(i) = 0; });
      team.team_barrier();

      // Step 5b: Multiply by the evaluation operator
      // the evaluation matrix
      auto                     startEval = evalOffsets(batch);
      BatchMatrix<MemorySpace> eval(&evalMat(startEval), outSize, inSize);
      // res := 1.0 * eval * b + 0.0 * res
      KokkosBlas::Experimental::Gemv<
          KokkosBlas::Mode::Team,
          KokkosBlas::Algo::Gemv::Blocked>::invoke(team, 'N', 1.0, eval, in, 0.0, res);

      team.team_barrier();

      // Step 5c: write the weightes result back to the global vector
      // potentially applying the polynomial term alongside
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, outSize),
          [&](int i) {
            auto   globalID = globalOutIDs(i + outBegin);
            double sum      = res(i);
            // Add polynomial portion to the output data: out += V * p
            if constexpr (polynomial) {
              // The "1"/constant term is the last value in the result
              // dim is here matrixCols - 1
              sum += qrCoeffs[dim];
              // ... and the linear polynomial
              for (int d = 0; d < dim; ++d) {
                sum += outMesh(globalID, d) * qrCoeffs[d];
              }
            }
            auto w = normalizedWeights(i + outBegin);
            Kokkos::atomic_add(&out(globalID), sum * w);
          }); // TeamThreadRange
    } else {
      // Alternative approach: do all in one go:
      // Step 5a: each thread takes care of one output vertex
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, outSize), [&](int r) {
            auto globalID = globalOutIDs(r + outBegin);

            // we first extract the vertex coordinates
            Kokkos::Array<double, 3> outVertex = {0., 0., 0.};
            for (int d = 0; d < dim; ++d) {
              outVertex[d] = outMesh(globalID, d);
            }

            // Step 5b: The matrix vector multiplication res = A * in
            // Accumulate partial dot product in a thread-parallel manner
            double sum = 0.0;
            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(team, inSize),
                [&](int c, double &localSum) {
                  // compute the local output coefficients
                  double dist = 0;
                  for (int d = 0; d < dim; ++d) {
                    double diff = outVertex[d] - localInMesh(c, d);
                    dist += diff * diff;
                  }
                  dist = Kokkos::sqrt(dist);
                  // Evaluate your RBF or similar function
                  double val = f(dist, rbf_params);
                  localSum += val * in(c);
                },
                sum); // ThreadVectorRange

            // Step 5c: if we have the polynomial as well, we have to apply it here
            if constexpr (polynomial) {
              // The "1"/constant term is the last value in the result
              // dim is here matrixCols - 1
              sum += qrCoeffs[dim];
              // ... and the linear polynomial
              for (int d = 0; d < dim; ++d) {
                sum += outVertex[d] * qrCoeffs[d];
              }
            }
            // Step 5d: Store final (weighted result) in the global out vector
            auto w = normalizedWeights(r + outBegin);
            Kokkos::atomic_add(&out(globalID), sum * w);
          }); // End TeamThreadRange
    } // end if-merged-evaluation
    // End Team parallel loop
  };

  // We allocate one vector for indata and once for outdata
  // We need at least four entries for the polynomial, seems unlikely for maxInClusterSize to be lower
  // but we should be certain
  auto inBytes  = ScratchVector::shmem_size(std::max(4, maxInClusterSize));
  auto outBytes = ScratchVector::shmem_size(maxOutClusterSize);
  if (!evaluation_op_available || polynomial) {
    // and additional storage if we have the polynomial
    inBytes = 4 * inBytes;
  }
  // We use the outBytes only for the per-cluster output vector, which
  // we don't need if we evaluate everything on the fly
  if (!evaluation_op_available) {
    outBytes = 0;
  }

  // We put the solution and the in data values into shared memory
  // TODO: Avoid the duplicate memory definitions here, currently needed as we
  // cannot change the Kokkos::AUTO to the actual recommendataion later
  auto tmpPol = TeamPolicy(nCluster, Kokkos::AUTO)
                    .set_scratch_size(
                        /* level = */ 0, Kokkos::PerTeam(inBytes))
                    .set_scratch_size(
                        /* level = */ 1, Kokkos::PerTeam(outBytes));

  auto teamSize = impl::findTeamSize<ExecSpace>(avgInClusterSize, kernel, tmpPol);

  auto policy = TeamPolicy(nCluster, teamSize)
                    .set_scratch_size(
                        /* level = */ 0, Kokkos::PerTeam(inBytes))
                    .set_scratch_size(
                        /* level = */ 1, Kokkos::PerTeam(outBytes));

  Kokkos::parallel_for("do_batched_solve", policy, kernel);
}

// Currently not used at all, solves the QR decomposition used for the polynomial
// but cannot be applied standalone, as we subtract the polynomial part on a per-cluster
// basis from the inData. Still useful for development purposes
template <typename MemorySpace>
void do_qr_solve(int                           nCluster,
                 int                           dim,
                 int                           maxInClusterSize,
                 VectorOffsetView<MemorySpace> inOffsets,
                 GlobalIDView<MemorySpace>     globalInIDs,
                 VectorView<MemorySpace>       inData,
                 MeshView<MemorySpace>         inMesh,
                 VectorView<MemorySpace>       qrMatrix,
                 VectorView<MemorySpace>       qrTau,
                 PivotView<MemorySpace>        qrP,
                 const VectorView<MemorySpace> weights,
                 VectorOffsetView<MemorySpace> outOffsets,
                 GlobalIDView<MemorySpace>     globalOutIDs,
                 VectorView<MemorySpace>       outData,
                 MeshView<MemorySpace>         outMesh)
{
  using TeamPolicy  = Kokkos::TeamPolicy<MemorySpace>;
  using MemberType  = typename TeamPolicy::member_type;
  using ScratchView = Kokkos::View<double *[2], typename MemorySpace::scratch_memory_space, UnmanagedMemory>;

  // Used for the inData solution vector and the QR solving
  auto scratchSize = ScratchView::shmem_size(2 * maxInClusterSize);
  Kokkos::parallel_for("do_qr_solve", TeamPolicy(nCluster, Kokkos::AUTO).set_scratch_size(
                                          /* level = */ 0, Kokkos::PerTeam(scratchSize)),
                       KOKKOS_LAMBDA(const MemberType &team) {
                         // Step 1: Some pointers we need
                         const int batch = team.league_rank();
                         // Ranges
                         const auto inBegin  = inOffsets(batch);
                         const auto inEnd    = inOffsets(batch + 1);
                         const int  inSize   = inEnd - inBegin;
                         const auto outBegin = outOffsets(batch);
                         const auto outEnd   = outOffsets(batch + 1);
                         const int  outSize  = outEnd - outBegin;

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
                         const int            matrixCols  = dim + 1;
                         const offset_1d_type matrixBegin = inBegin * matrixCols;
                         const offset_1d_type tauBegin    = batch * matrixCols;
                         const offset_1d_type PBegin      = batch * (matrixCols + 1);
                         const int            rank        = qrP(PBegin + matrixCols);

                         BatchMatrix<MemorySpace>           qr(&qrMatrix(matrixBegin), inSize, matrixCols);
                         BatchVector<double *, MemorySpace> tau(&qrTau(tauBegin), matrixCols);
                         BatchVector<int *, MemorySpace>    P(&qrP(PBegin), matrixCols);

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
                             [&](offset_1d_type i) {
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
} // namespace precice::mapping::kernel
