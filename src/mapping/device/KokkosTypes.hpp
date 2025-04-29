#pragma once

#include <Kokkos_Core.hpp>
#include <cstdint>

namespace precice::mapping {

// Alias definitions

using ExecutionSpace  = Kokkos::DefaultExecutionSpace;
using UnmanagedMemory = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

using offset_1d_type = ExecutionSpace::size_type;
using offset_2d_type = ExecutionSpace::size_type;


// For the meshes, we keep the last index contiguous in memory, PUM has anyway random access
// in the global vectors
template <typename MemorySpace = ExecutionSpace>
using MeshView = Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace>;

// For all global data structures we store
template <typename MemorySpace = ExecutionSpace>
using VectorView = Kokkos::View<double *, MemorySpace>;

template <typename MemorySpace = ExecutionSpace>
using PivotView = Kokkos::View<int *, MemorySpace>;

// The VertexIDs mapping the cluster vertices to global IDs
template <typename MemorySpace = ExecutionSpace>
using GlobalIDView = Kokkos::View<VertexID *, MemorySpace>;

// Offsets for the flat views for 1D data (essentially for double **, where the last * is for the batches)
template <typename MemorySpace = ExecutionSpace>
using VectorOffsetView = Kokkos::View<offset_1d_type *, MemorySpace>;

// Offsets for the flat views for 2D data (essentially for double ***, where the last * is for the batches)
template <typename MemorySpace = ExecutionSpace>
using MatrixOffsetView = Kokkos::View<offset_2d_type *, MemorySpace>;
} // namespace precice::mapping
