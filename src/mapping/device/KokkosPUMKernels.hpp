#pragma once

#include <Kokkos_Core.hpp>
#include "mapping/device/KokkosTypes.hpp"
#include "mapping/impl/BasisFunctions.hpp"
namespace precice::mapping::kernel {

template <typename MemorySpace>
void compute_offsets(const VectorOffsetView<MemorySpace> src1, const VectorOffsetView<MemorySpace> src2,
                     MatrixOffsetView<MemorySpace> dst, int nCluster);

// returns true, if successful, currently not tuned as it is not performance critical
template <typename MemorySpace>
bool compute_weights(const int                     nCluster,
                     const offset_1d_type          nWeights,
                     const int                     nMeshVertices,
                     const int                     dim,
                     VectorOffsetView<MemorySpace> offsets,
                     MeshView<MemorySpace>         centers,
                     GlobalIDView<MemorySpace>     globalIDs,
                     MeshView<MemorySpace>         mesh,
                     const CompactPolynomialC2    &w,
                     VectorView<MemorySpace>       normalizedWeights);

// the input assembly
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
    VectorView<MemorySpace>              matrices);

template <typename EvalFunctionType, typename MemorySpace>
void do_batched_assembly(int                                  nCluster, // Number of local systems
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
                         VectorView<MemorySpace>              matrices);

template <typename MemorySpace>
void do_batched_qr(int                           nCluster,
                   int                           dim,
                   int                           avgClusterSize,
                   int                           maxClusterSize,
                   VectorOffsetView<MemorySpace> inOffsets,
                   GlobalIDView<MemorySpace>     globalInIDs,
                   MeshView<MemorySpace>         inMesh,
                   VectorView<MemorySpace>       qrMatrix,
                   VectorView<MemorySpace>       qrTau,
                   PivotView<MemorySpace>        qrP);

// currently unused
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
                 MeshView<MemorySpace>         outMesh);

template <typename MemorySpace>
void do_batched_lu(
    int                                  nCluster,
    int                                  avgClusterSize,
    const MatrixOffsetView<MemorySpace> &matrixOffsets,
    VectorView<MemorySpace>              matrices);

template <bool polynomial, typename MemorySpace>
void do_batched_solve(
    int                                  nCluster,
    int                                  dim,
    int                                  avgInClusterSize,
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
    const MeshView<MemorySpace>         &inMesh,
    const MeshView<MemorySpace>         &outMesh,
    const VectorView<MemorySpace>       &qrMatrix,
    const VectorView<MemorySpace>       &qrTau,
    const PivotView<MemorySpace>        &qrP);

} // namespace precice::mapping::kernel

#include "mapping/device/KokkosPUMKernels_Impl.hpp"
