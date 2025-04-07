#pragma once
#ifndef PRECICE_NO_GINKGO

#include <array>
#include <cmath>
#include <functional>
// #include <ginkgo/extensions/kokkos.hpp>
#include <Kokkos_Core.hpp>
#include <numeric>
#include "mapping/GinkgoDefinitions.hpp"
#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/device/Ginkgo.hpp"
#include "mapping/device/GinkgoRBFKernels.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/impl/SphericalVertexCluster.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#ifdef PRECICE_WITH_HIP
#include "mapping/device/HipQRSolver.hip.hpp"
#endif
#ifdef PRECICE_WITH_CUDA
#include "mapping/device/CudaQRSolver.cuh"
#endif
#ifdef PRECICE_WITH_OPENMP
#include <omp.h>
#endif

using precice::mapping::RadialBasisParameters;

namespace precice {
namespace mapping {

// Runtime lookups as suggested by Ginkgo

/**
 * This class assembles and solves an RBF system, given an input mesh and an output mesh with relevant vertex IDs.
 * It uses iterative solvers (CG, GMRES) and preconditioners ((Block-)Jacobi, Cholesky, Ilu) to solve the interpolation
 * systems. Furthermore, it optionally does that on Nvidia or AMD GPUs which provides significant speedup over (single-threaded)
 * CPU implementations.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class BatchedRBFSolver {
public:
  using RBF_T = RADIAL_BASIS_FUNCTION_T;

  /// Essentially an initialize of the solver: allocates device memory
  /// and computes the Cholesky decompositions
  BatchedRBFSolver(RBF_T                                             basisFunction,
                   mesh::PtrMesh                                     inMesh,
                   mesh::PtrMesh                                     outMesh,
                   const std::vector<SphericalVertexCluster<RBF_T>> &clusters,
                   Polynomial                                        polynomial,
                   MappingConfiguration::GinkgoParameter             ginkgoParameter);

  void clear();

private:
  mutable precice::logging::Logger _log{"mapping::BatchedRBFSolver"};

  // std::shared_ptr<gko::Executor>                     _deviceExecutor;
  // std::shared_ptr<gko::Executor>                     _hostExecutor = gko::ReferenceExecutor::create();
  Kokkos::View<int *, Kokkos::DefaultExecutionSpace> _inOffsets;
  Kokkos::View<int *, Kokkos::DefaultExecutionSpace> _outOffsets;

  Kokkos::View<std::size_t *, Kokkos::DefaultExecutionSpace> _kernelOffsets;
  Kokkos::View<std::size_t *, Kokkos::DefaultExecutionSpace> _evaluationOffsets;

  // Essentially a MatrixXd, where the last index is contiguous in memory (important for the assembly kernel)
  Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> _inMesh;
  Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> _outMesh;

  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _kernelMatrices;

  Polynomial _polynomial;
  void       _solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const;
  // MappingConfiguration::GinkgoParameter _ginkgoParameter;
};

template <typename RADIAL_BASIS_FUNCTION_T>
BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::BatchedRBFSolver(RBF_T                                             basisFunction,
                                                            mesh::PtrMesh                                     inMesh,
                                                            mesh::PtrMesh                                     outMesh,
                                                            const std::vector<SphericalVertexCluster<RBF_T>> &clusters,
                                                            Polynomial                                        polynomial,
                                                            MappingConfiguration::GinkgoParameter             ginkgoParameter)
    : _polynomial(polynomial)
{
  PRECICE_TRACE();
  PRECICE_CHECK(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && _polynomial == Polynomial::ON), "Not supported.");
  PRECICE_CHECK(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite(), "Using a Cholesky decomposition.");
  PRECICE_CHECK(!(inMesh->vertices().empty() || outMesh->vertices().empty()), "One of the meshes in the batched solvers is empty, which is invalid.");
  PRECICE_CHECK(inMesh->getDimensions() == outMesh->getDimensions(), "Incompatible dimensions passed to the batched solver.");

  // We have to initialize Kokkos and Ginkgo here, as the initialization call allocates memory
  // in the current setup, this will only initialize the device (and allocate memory) on the primary rank
  device::Ginkgo::initialize(ginkgoParameter.nThreads, ginkgoParameter.deviceId);
  PRECICE_INFO("Using Batched solver on executor {}", ginkgoParameter.executor);
  // _deviceExecutor = create_device_executor(ginkgoParameter.executor, ginkgoParameter.enableUnifiedMemory);
#ifdef PRECICE_WITH_OPENMP
  if (ginkgoParameter.nThreads > 0 && ginkgoParameter.executor == "omp-executor")
    omp_set_num_threads(ginkgoParameter.nThreads);
#endif

  // Step 1: compute offsets, which hold the range for each cluster
  PRECICE_DEBUG("Computing mesh offsets");
  const std::size_t nCluster = clusters.size();

  // Initialize the view
  _inOffsets  = Kokkos::View<int *>("inOffsets", nCluster + 1);
  _outOffsets = Kokkos::View<int *>("outOffsets", nCluster + 1);

  // we fill this view on the host side
  auto hostIn  = Kokkos::create_mirror_view(_inOffsets);
  auto hostOut = Kokkos::create_mirror_view(_outOffsets);

  // has to be a separate loop, as we first need to gather knowledge about
  // the shape for the meshes
  hostIn(0)  = 0;
  hostOut(0) = 0;
  for (std::size_t i = 0; i < nCluster; ++i) {
    hostIn(i + 1)  = hostIn(i) + clusters[i].getNumberOfInputVertices();
    hostOut(i + 1) = hostOut(i) + clusters[i].getNumberOfOutputVertices();
  }

  // Copy offsets onto the device
  Kokkos::deep_copy(_inOffsets, hostIn);
  Kokkos::deep_copy(_outOffsets, hostOut);

  // Step 2: Compute the matrix offsets on the device
  PRECICE_DEBUG("Computing matrix offsets");
  _kernelOffsets     = Kokkos::View<std::size_t *>("kernelOffsets", nCluster + 1);
  _evaluationOffsets = Kokkos::View<std::size_t *>("evaluationOffsets", nCluster + 1);

  Kokkos::deep_copy(_kernelOffsets, 0);
  Kokkos::deep_copy(_evaluationOffsets, 0);

  // We use a parallel scan for that
  kernel::compute_offsets(_inOffsets, _inOffsets, _kernelOffsets, nCluster);
  kernel::compute_offsets(_inOffsets, _outOffsets, _evaluationOffsets, nCluster);

  // Step 3: Handle the mesh data structure and copy over to the device
  PRECICE_DEBUG("Computing mesh data on the device");
  const auto dim = inMesh->getDimensions();
  _inMesh        = Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>("inMesh", hostIn(nCluster), dim);
  _outMesh       = Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>("outMesh", hostOut(nCluster), dim);

  auto hostInMesh  = Kokkos::create_mirror_view(_inMesh);
  auto hostOutMesh = Kokkos::create_mirror_view(_outMesh);

  // fill to account for local dead axis
  Kokkos::deep_copy(hostInMesh, 0.0);
  Kokkos::deep_copy(hostOutMesh, 0.0);

  Eigen::Index inIndex  = 0;
  Eigen::Index outIndex = 0;
  for (std::size_t i = 0; i < nCluster; ++i) {
    const Eigen::MatrixXd Q = clusters[i].getLocalPolynomialInputMatrix(inMesh);

    for (int i = 0; i < Q.rows(); ++i, ++inIndex) {
      for (int d = 0; d < dim; ++d) {
        hostInMesh(inIndex, d) = Q(i, d + 1);
      }
    }

    const Eigen::MatrixXd V = clusters[i].getLocalPolynomialOutputMatrix(outMesh);
    PRECICE_ASSERT(Q.cols() == V.cols());

    for (int i = 0; i < V.rows(); ++i, ++outIndex) {
      for (int d = 0; d < dim; ++d) {
        hostOutMesh(outIndex, d) = V(i, d + 1);
      }
    }
  }

  // Copy to device
  Kokkos::deep_copy(_inMesh, hostInMesh);
  Kokkos::deep_copy(_outMesh, hostOutMesh);

  // Step 4: Launch the parallel kernel to assemble the kernel matrices
  PRECICE_DEBUG("Assemble batched matrices");
  std::size_t unrolledSize   = 0;
  auto        last_elem_view = Kokkos::subview(_kernelOffsets, nCluster);
  Kokkos::deep_copy(unrolledSize, last_elem_view);
  _kernelMatrices = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>("batched_kernels", unrolledSize);

  kernel::do_batched_assembly(nCluster, dim, basisFunction, basisFunction.getFunctionParameters(),
                              _inOffsets, _inMesh, _inOffsets, _inMesh, _kernelOffsets, _kernelMatrices);

  kernel::do_batched_lu(nCluster, _kernelOffsets, _kernelMatrices);
  // Maybe simply add a function to the cluster to give as the polynomial matrices
  // maybe compute those matrices during the initialization (as we have the mesh then)

  //   const std::size_t deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  //   const std::size_t dimensions     = 3;
  //   const std::size_t polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  //   // Add linear polynom degrees if polynomial requires this
  //   const auto inputSize  = inputIDs.size();
  //   const auto outputSize = outputIDs.size();
  //   const auto n          = inputSize + polyparams;

  //   PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  //   PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  //   const std::size_t inputMeshSize  = inputMesh.nVertices();
  //   const std::size_t outputMeshSize = outputMesh.nVertices();
  //   const std::size_t meshDim        = inputMesh.vertex(0).getDimensions();

  //   _scalarOne         = gko::share(gko::initialize<GinkgoScalar>({1.0}, _deviceExecutor));
  //   _scalarNegativeOne = gko::share(gko::initialize<GinkgoScalar>({-1.0}, _deviceExecutor));

  //   // Now we fill the RBF system matrix on the GPU (or any other selected device)
  //   precice::profiling::Event _allocCopyEvent{"map.rbf.ginkgo.memoryAllocAndCopy"};
  //   _rbfCoefficients = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{n, 1}));
  //   _allocCopyEvent.stop();
  //   // Initial guess is required since uninitialized memory could lead to a never converging system
  //   _rbfCoefficients->fill(0.0);

  //   // We need to copy the input data into a CPU stored vector first and copy it to the GPU afterwards
  //   // To allow for coalesced memory accesses on the GPU, we need to store them in transposed order IFF the backend is the GPU
  //   // However, the CPU does not need that; in fact, it would make it slower
  //   std::size_t inputVerticesM, inputVerticesN, outputVerticesM, outputVerticesN;

  //   if ("cuda-executor" == ginkgoParameter.executor || "hip-executor" == ginkgoParameter.executor) {
  //     inputVerticesM  = meshDim;
  //     inputVerticesN  = inputMeshSize;
  //     outputVerticesM = meshDim;
  //     outputVerticesN = outputMeshSize;
  //   } else {
  //     inputVerticesM  = inputMeshSize;
  //     inputVerticesN  = meshDim;
  //     outputVerticesM = outputMeshSize;
  //     outputVerticesN = meshDim;
  //   }

  //   auto inputVertices  = gko::share(GinkgoMatrix::create(_hostExecutor, gko::dim<2>{inputVerticesM, inputVerticesN}));
  //   auto outputVertices = gko::share(GinkgoMatrix::create(_hostExecutor, gko::dim<2>{outputVerticesM, outputVerticesN}));
  //   for (std::size_t i = 0; i < inputMeshSize; ++i) {
  //     for (std::size_t j = 0; j < meshDim; ++j) {
  //       if ("cuda-executor" == ginkgoParameter.executor || "hip-executor" == ginkgoParameter.executor) {
  //         inputVertices->at(j, i) = inputMesh.vertex(i).coord(j);
  //       } else {
  //         inputVertices->at(i, j) = inputMesh.vertex(i).coord(j);
  //       }
  //     }
  //   }
  //   for (std::size_t i = 0; i < outputMeshSize; ++i) {
  //     for (std::size_t j = 0; j < meshDim; ++j) {
  //       if ("cuda-executor" == ginkgoParameter.executor || "hip-executor" == ginkgoParameter.executor) {
  //         outputVertices->at(j, i) = outputMesh.vertex(i).coord(j);
  //       } else {
  //         outputVertices->at(i, j) = outputMesh.vertex(i).coord(j);
  //       }
  //     }
  //   }

  //   _allocCopyEvent.start();

  //   auto dInputVertices  = gko::clone(_deviceExecutor, inputVertices);
  //   auto dOutputVertices = gko::clone(_deviceExecutor, outputVertices);
  //   inputVertices->clear();
  //   outputVertices->clear();

  //   _deviceExecutor->synchronize();

  //   _rbfSystemMatrix = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{n, n}));
  //   _matrixA         = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{outputSize, n}));

  //   _allocCopyEvent.stop();

  //   // Launch RBF fill kernel on device
  //   precice::profiling::Event _assemblyEvent{"map.rbf.ginkgo.assembleMatrices"};
  //   precice::profiling::Event systemMatrixAssemblyEvent{"map.rbf.ginkgo.assembleSystemMatrix"};
  //   kernel::create_rbf_system_matrix(_deviceExecutor, ginkgoParameter.enableUnifiedMemory, _rbfSystemMatrix, activeAxis, dInputVertices, dInputVertices, basisFunction,
  //                                    basisFunction.getFunctionParameters(), Polynomial::ON == polynomial,
  //                                    polyparams); // polynomial evaluates to true only if ON is set
  //   _deviceExecutor->synchronize();
  //   systemMatrixAssemblyEvent.stop();

  //   precice::profiling::Event outputMatrixAssemblyEvent{"map.rbf.ginkgo.assembleOutputMatrix"};
  //   kernel::create_rbf_system_matrix(_deviceExecutor, ginkgoParameter.enableUnifiedMemory, _matrixA, activeAxis, dInputVertices, dOutputVertices, basisFunction,
  //                                    basisFunction.getFunctionParameters(), Polynomial::ON == polynomial, polyparams);

  //   // Wait for the kernels to finish
  //   _deviceExecutor->synchronize();
  //   outputMatrixAssemblyEvent.stop();
  //   _assemblyEvent.stop();

  //   dInputVertices->clear();
  //   dOutputVertices->clear();

  //   const std::size_t M = _rbfSystemMatrix->get_size()[0];
  //   const std::size_t N = _rbfSystemMatrix->get_size()[1];
  //   _decompMatrixQ_T    = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>(N, M)));
  //   _decompMatrixR      = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>(N, N)));

  //   if ("cuda-executor" == ginkgoParameter.executor) {
  // #ifdef PRECICE_WITH_CUDA
  //     // _rbfSystemMatrix will be overridden into Q
  //     computeQRDecompositionCuda(_deviceExecutor, _rbfSystemMatrix.get(), _decompMatrixR.get());
  // #endif
  //   } else if ("hip-executor" == ginkgoParameter.executor) {
  // #ifdef PRECICE_WITH_HIP
  //     // _rbfSystemMatrix will be overridden into Q
  //     computeQRDecompositionHip(_deviceExecutor, _rbfSystemMatrix.get(), _decompMatrixR.get());
  // #endif
  //   } else {
  //     PRECICE_UNREACHABLE("Not implemented");
  //   }
  //   _rbfSystemMatrix->transpose(_decompMatrixQ_T);
  //   _dQ_T_Rhs = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_decompMatrixQ_T->get_size()[0], 1}));
}

template <typename RADIAL_BASIS_FUNCTION_T>
void BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::_solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const
{
  PRECICE_TRACE();
}

template <typename RADIAL_BASIS_FUNCTION_T>
void BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::clear()
{
}

} // namespace mapping
} // namespace precice

#endif // PRECICE_NO_GINKGO
