#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>
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

  // Required as workspace for the pivoted QR, see code comment
  // workspace (norm and householder application, 2 * max(m,n) is needed)
  auto scratchSize = ScratchView::shmem_size(2 * maxClusterSize);
  Kokkos::parallel_for("do_batched_qr", TeamPolicy(nCluster, Kokkos::AUTO).set_scratch_size(
                                            /* level = */ 0, Kokkos::PerTeam(scratchSize)),
                       KOKKOS_LAMBDA(const MemberType &team) {
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
                       });
}

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

// TODO: Using Kokkos::LayoutRight for the Coords performs a bit better for the assembly,
// but it might deteriorate performance related to the polynomials. Especially for the gemv
// for the polynomial contributions etc
// For the full GPU porting, the Layout can only be LayoutRight, as we don't access the coordinates
// coalesced, but rather first pick what we need
template <typename EvalFunctionType, typename MemorySpace>
void do_batched_assembly(
    int                                       nCluster, // Number of local systems
    int                                       dim,      // Dimension of points
    EvalFunctionType                          f,
    ::precice::mapping::RadialBasisParameters rbf_params,
    const VectorOffsetView<MemorySpace>      &inOffsets, // vertex offsets (length N+1)
    const GlobalIDView<MemorySpace>          &globalInIDs,
    const MeshView<MemorySpace>              &inCoords, // meshes
    const VectorOffsetView<MemorySpace>      &targetOffsets,
    const GlobalIDView<MemorySpace>          &globalTargetIDs,
    const MeshView<MemorySpace>              &targetCoords,
    const MatrixOffsetView<MemorySpace>      &matrixOffsets,
    VectorView<MemorySpace>                   matrices) // 1D view of batched matrices
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using TeamMember = typename TeamPolicy::member_type;

  // We launch one team per local system
  Kokkos::parallel_for("do_batched_assembly", TeamPolicy(nCluster, Kokkos::AUTO), KOKKOS_LAMBDA(const TeamMember &team) {
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
        Kokkos::TeamThreadRange(team, nrows),
        [=](int r) {
          Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(team, ncols),
              [=](int c) {
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
        });       // TeamThreadRange
  });
}

template <typename MemorySpace>
void do_batched_lu(
    int                                  nCluster,
    const MatrixOffsetView<MemorySpace> &matrixOffsets,
    VectorView<MemorySpace>              matrices)
{
  using ExecSpace  = typename MemorySpace::execution_space;
  using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType = typename TeamPolicy::member_type;

  Kokkos::parallel_for("do_batched_lu", TeamPolicy(nCluster, Kokkos::AUTO), KOKKOS_LAMBDA(const MemberType &team) {
    const int i     = team.league_rank();
    auto      start = matrixOffsets(i);
    auto      end   = matrixOffsets(i + 1);
    int       n     = static_cast<int>(Kokkos::sqrt(end - start));

    BatchMatrix<MemorySpace> A(&matrices(start), n, n);

    KokkosBatched::TeamLU<MemberType, KokkosBatched::Algo::LU::Blocked>::invoke(team, A);
    // Parallel end
  });
}

template <bool polynomial, typename MemorySpace>
void do_batched_solve(
    int                                  nCluster,
    int                                  dim,
    int                                  maxInClusterSize,
    int                                  maxOutClusterSize,
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

  using ScratchSpace  = typename MemorySpace::scratch_memory_space;
  using ScratchView1d = Kokkos::View<double *[1], ScratchSpace, UnmanagedMemory>;
  using ScratchView4d = Kokkos::View<double *[4], ScratchSpace, UnmanagedMemory>;
  using ScratchVector = Kokkos::View<double *, ScratchSpace, UnmanagedMemory>;
  using ScratchMatrix = std::conditional_t<polynomial, ScratchView4d, ScratchView1d>;

  // Once for indata and once for outdata
  // We need at least four entries for the polynomial, seems unlikely for maxInClusterSize to be lower
  // but we should be certain
  auto       inBytes  = ScratchVector::shmem_size(std::max(4, maxInClusterSize));
  auto       outBytes = ScratchVector::shmem_size(maxOutClusterSize);
  TeamPolicy policy(nCluster, Kokkos::AUTO);

  // We put the solution and the in data values into shared memory
  if (polynomial) {
    policy.set_scratch_size(
        /* level = */ 0, Kokkos::PerTeam(4 * inBytes));
  } else {
    policy.set_scratch_size(
        /* level = */ 0, Kokkos::PerTeam(inBytes));
  }
  // Put the out vector into level 1
  policy.set_scratch_size(
      /* level = */ 1, Kokkos::PerTeam(outBytes));
  Kokkos::parallel_for("do_batched_solve", policy, KOKKOS_LAMBDA(const MemberType &team) {
    // Required for correct capturing, as these variables are only conditionally used further down
    (void) dim;
    (void) qrMatrix;
    (void) qrTau;
    (void) qrP;
    (void) inMesh;
    (void) outMesh;

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
      auto in_cp = Kokkos::subview(work, std::pair<int, int>(0, inSize), 1);
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, inSize), [&](int i) { in_cp(i) = in(i); });
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
      // tmp size might be insufficient: there was no size requirement specified for the workspace
      auto tmp = Kokkos::subview(work, Kokkos::ALL, 2);
      KokkosBatched::TeamVectorApplyQ<MemberType,
                                      KokkosBatched::Side::Left,
                                      KokkosBatched::Trans::Transpose,
                                      KokkosBatched::Algo::ApplyQ::Unblocked>::invoke(team, qr, tau, in_cp, tmp);
      team.team_barrier();

      // Step 3d: Solve triangular solve R z = y
      auto in_r = Kokkos::subview(in_cp, std::pair<int, int>(0, rank));
      auto R    = Kokkos::subview(qr, std::pair<int, int>(0, rank), std::pair<int, int>(0, rank));
      KokkosBatched::Trsv<
          MemberType,
          KokkosBatched::Uplo::Upper,
          KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Diag::NonUnit,
          KokkosBatched::Mode::Team,
          KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, R, in_r);
      team.team_barrier();

      // Steo 3e: zero out the entries which are not within the rank region
      auto poly = Kokkos::subview(work, std::pair<int, int>(0, matrixCols), 1);
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, rank, matrixCols), [&](int i) { poly(i) = 0; });

      team.team_barrier();

      // Step 3f: Apply pivoting x = P z
      KokkosBatched::TeamVectorApplyPivot<MemberType,
                                          KokkosBatched::Side::Left,
                                          KokkosBatched::Direct::Backward>::invoke(team, P, poly);
      team.team_barrier();

      // Step 3g: store the result in the array
      for (int c = 0; c < matrixCols; ++c) {
        qrCoeffs[c] = poly(c);
      }
      team.team_barrier();

      // Step 3h: Subtract polynomial portion from the input data: in -= Q * p
      // threading over inSize
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, inSize),
          [&](int i) {
            auto globalID = globalRhsIDs(i + inBegin);
            // The "1"/constant term is the last value in the result
            // dim is here matrixCols - 1
            double sum = qrCoeffs[dim];
            // ... and the linear polynomial
            for (int d = 0; d < dim; ++d)
              sum += inMesh(globalID, d) * qrCoeffs[d];
            in(i) -= sum;
          });
      team.team_barrier();
    }

    // Step 4: Solve the LU decomposition
    // The lu inplace lu decomposition computed with KokkosBatched
    auto                     matStart = matrixOffsets(batch);
    BatchMatrix<MemorySpace> A(&matrices(matStart), inSize, inSize);

    // Convenience routine for LU, which performs the
    // Forward substitution: solve L * y = b and
    // Backward substitution: solve U * x = y inplace
    KokkosBatched::SolveLU<MemberType,
                           KokkosBatched::Trans::NoTranspose,
                           KokkosBatched::Mode::Team,
                           KokkosBatched::Algo::SolveLU::Blocked>::invoke(team, A, in);

    team.team_barrier();

    // Step 5: Allocate and zero out a local result vector (more of a safety feature)
    ScratchVector res(team.team_scratch(1), outSize);
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, outSize),
        [&](int i) { res(i) = 0; });
    team.team_barrier();

    // Step 6: Multiply by the evaluation operator
    // the evaluation matrix
    auto                     startEval = evalOffsets(batch);
    BatchMatrix<MemorySpace> eval(&evalMat(startEval), outSize, inSize);

    // res := 1.0 * eval * b + 0.0 * res
    KokkosBlas::Experimental::Gemv<
        KokkosBlas::Mode::Team,
        KokkosBlas::Algo::Gemv::Blocked>::invoke(team, 'N', 1.0, eval, in, 0.0, res);

    team.team_barrier();

    // Step 7: write the result back to the global vector
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
    // End Team parallel loop
  });
}
} // namespace precice::mapping::kernel
