#pragma once

#include <Kokkos_Core.hpp>
#include "mapping/impl/BasisFunctions.hpp"

namespace precice::mapping::kernel {

void compute_offsets(const Kokkos::View<int *> src1, const Kokkos::View<int *> src2,
                     Kokkos::View<std::size_t *> dst, int N);

// returns true, if successful
bool compute_weights(const std::size_t                                                           nCenters,
                     const std::size_t                                                           nWeights,
                     const std::size_t                                                           nMeshVertices,
                     const int                                                                   dim,
                     Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          offsets,
                     Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> centers,
                     Kokkos::View<int *, Kokkos::DefaultExecutionSpace>                          globalIDs,
                     Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> mesh,
                     const CompactPolynomialC2                                                  &w,
                     Kokkos::View<double *, Kokkos::DefaultExecutionSpace>                       normalizedWeights);

template <typename EvalFunctionType, typename MemorySpace>
void do_batched_assembly(int                                                              N,   // Number of local systems
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
                         Kokkos::View<double *, MemorySpace>                              matrices);

template <typename MemorySpace>
void do_batched_qr(std::size_t                                               nCluster,
                   int                                                       dim,
                   int                                                       maxClusterSize,
                   Kokkos::View<int *, MemorySpace>                          inOffsets,
                   Kokkos::View<int *, MemorySpace>                          globalInIDs,
                   Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> inMesh,
                   Kokkos::View<double *, MemorySpace>                       qrMatrix,
                   Kokkos::View<double *, MemorySpace>                       qrTau,
                   Kokkos::View<int *, MemorySpace>                          qrP);

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
                 Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> outMesh);

template <typename MemorySpace>
void do_batched_lu(
    int                                        N,
    const Kokkos::View<size_t *, MemorySpace> &matrixOffsets,
    Kokkos::View<double *, MemorySpace>        matrices);

template <bool polynomial, typename MemorySpace>
void do_batched_solve(
    std::size_t                                                      N,
    int                                                              dim,
    std::size_t                                                      maxInClusterSize,
    std::size_t                                                      maxOutClusterSize,
    const Kokkos::View<int *, MemorySpace>                          &rhsOffsets,
    const Kokkos::View<int *, MemorySpace>                          &globalRhsIDs,
    Kokkos::View<double *, MemorySpace>                              rhs,
    const Kokkos::View<size_t *, MemorySpace>                       &matrixOffsets,
    const Kokkos::View<double *, MemorySpace>                       &matrices,
    const Kokkos::View<double *, MemorySpace>                       &normalizedWeights,
    const Kokkos::View<size_t *, MemorySpace>                       &evalOffsets,
    const Kokkos::View<double *, MemorySpace>                       &evalMat,
    const Kokkos::View<int *, MemorySpace>                          &outOffsets,
    const Kokkos::View<int *, MemorySpace>                          &globalOutIDs,
    Kokkos::View<double *, MemorySpace>                              out,
    const Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> &inMesh,
    const Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> &outMesh,
    const Kokkos::View<double *, MemorySpace>                       &qrMatrix,
    const Kokkos::View<double *, MemorySpace>                       &qrTau,
    const Kokkos::View<int *, MemorySpace>                          &qrP);

} // namespace precice::mapping::kernel
