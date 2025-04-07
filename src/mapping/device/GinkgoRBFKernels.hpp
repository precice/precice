#pragma once

#include <Kokkos_Core.hpp>
#include "mapping/GinkgoDefinitions.hpp"
#include "mapping/impl/BasisFunctions.hpp"

namespace precice {
namespace mapping {

std::shared_ptr<gko::Executor> create_device_executor(const std::string &execName, bool enableUnifiedMemory);

namespace kernel {

template <typename EvalFunctionType>
void create_rbf_system_matrix(std::shared_ptr<const gko::Executor> exec,
                              bool                                 unifiedMemory,
                              gko::ptr_param<GinkgoMatrix> mtx, const std::array<bool, 3> activeAxis,
                              gko::ptr_param<GinkgoMatrix> supportPoints, gko::ptr_param<GinkgoMatrix> targetPoints,
                              EvalFunctionType f, ::precice::mapping::RadialBasisParameters rbf_params, bool addPolynomial, unsigned int extraDims = 0);

void fill_polynomial_matrix(std::shared_ptr<const gko::Executor> exec,
                            bool                                 unifiedMemory,
                            gko::ptr_param<GinkgoMatrix> mtx, gko::ptr_param<const GinkgoMatrix> x, const unsigned int dims);

void compute_offsets(const Kokkos::View<int *> src1, const Kokkos::View<int *> src2,
                     Kokkos::View<std::size_t *> dst, int N);

template <typename EvalFunctionType, typename MemorySpace>
void do_batched_assembly(int                                                              N,   // Number of local systems
                         int                                                              dim, // Dimension of points
                         EvalFunctionType                                                 f,
                         ::precice::mapping::RadialBasisParameters                        rbf_params,
                         const Kokkos::View<int *, MemorySpace>                          &inOffsets, // vertex offsets (length N+1)
                         const Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> &inCoords,  // meshes
                         const Kokkos::View<int *, MemorySpace>                          &targetOffsets,
                         const Kokkos::View<double **, Kokkos::LayoutRight, MemorySpace> &targetCoords,
                         const Kokkos::View<size_t *, MemorySpace>                       &matrixOffsets,
                         Kokkos::View<double *, MemorySpace>                              matrices);

template <typename MemorySpace>
void do_batched_lu(
    int                                        N,
    const Kokkos::View<size_t *, MemorySpace> &matrixOffsets,
    Kokkos::View<double *, MemorySpace>        matrices);

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
    Kokkos::View<double *, MemorySpace>        out);

} // namespace kernel
} // namespace mapping
} // namespace precice
