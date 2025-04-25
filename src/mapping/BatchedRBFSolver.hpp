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
  BatchedRBFSolver(RBF_T                                 basisFunction,
                   mesh::PtrMesh                         inMesh,
                   mesh::PtrMesh                         outMesh,
                   const std::vector<mesh::Vertex>      &centers,
                   double                                clusterRadius,
                   Polynomial                            polynomial,
                   MappingConfiguration::GinkgoParameter ginkgoParameter);

  void clear();

  void solveConsistent(const time::Sample &globalIn, Eigen::VectorXd &globalOut);

private:
  mutable precice::logging::Logger _log{"mapping::BatchedRBFSolver"};

  // std::shared_ptr<gko::Executor>                     _deviceExecutor;
  // std::shared_ptr<gko::Executor>                     _hostExecutor = gko::ReferenceExecutor::create();
  Kokkos::View<int *, Kokkos::DefaultExecutionSpace> _inOffsets;
  Kokkos::View<int *, Kokkos::DefaultExecutionSpace> _outOffsets;

  Kokkos::View<int *, Kokkos::DefaultExecutionSpace> _globalInIDs;
  Kokkos::View<int *, Kokkos::DefaultExecutionSpace> _globalOutIDs;

  Kokkos::View<std::size_t *, Kokkos::DefaultExecutionSpace> _kernelOffsets;
  Kokkos::View<std::size_t *, Kokkos::DefaultExecutionSpace> _evaluationOffsets;

  // Essentially a MatrixXd, where the last index is contiguous in memory (important for the assembly kernel)
  Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> _inMesh;
  Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> _outMesh;

  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _qrMatrix; // flat view of (nCluster x verticesPerCluster_i x (dim + 1) = nCluster x verticesPerCluster_i x polyParams)
  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _qrTau;    // flat view of Householder tau (nCluster x (dim + 1) = nCluster x polyParams)
  Kokkos::View<int *, Kokkos::DefaultExecutionSpace>    _qrP;      // flat view of Permutation and rank (nCluster x (dim + 2) = nCluster x (polyParams + rank))

  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _kernelMatrices;

  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _evalMatrices;
  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _normalizedWeights;

  // Currently only scalar data
  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _inData;
  // Kokkos::View<double *>::HostMirror                    _inDataMirror;
  Kokkos::View<double *, Kokkos::DefaultExecutionSpace> _outData;
  // Kokkos::View<double *>::HostMirror                    _outDataMirror;

  std::size_t       _maxInClusterSize;
  std::size_t       _maxOutClusterSize;
  Polynomial        _polynomial;
  const std::size_t _nCluster;
  const int         _dim; // Mesh dimension
  // MappingConfiguration::GinkgoParameter _ginkgoParameter;
};

template <typename RADIAL_BASIS_FUNCTION_T>
BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::BatchedRBFSolver(RBF_T                                 basisFunction,
                                                            mesh::PtrMesh                         inMesh,
                                                            mesh::PtrMesh                         outMesh,
                                                            const std::vector<mesh::Vertex>      &centers,
                                                            double                                clusterRadius,
                                                            Polynomial                            polynomial,
                                                            MappingConfiguration::GinkgoParameter ginkgoParameter)
    : _polynomial(polynomial), _nCluster(centers.size()), _dim(inMesh->getDimensions())
{
  PRECICE_TRACE();
  PRECICE_CHECK(_polynomial != Polynomial::ON, "Not supported.");
  // TODO: check the positive definiteness again later
  // PRECICE_CHECK(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite(), "Using a LU decomposition.");
  PRECICE_CHECK(!(inMesh->vertices().empty() || outMesh->vertices().empty()), "One of the meshes in the batched solvers is empty, which is invalid.");
  PRECICE_CHECK(inMesh->getDimensions() == outMesh->getDimensions(), "Incompatible dimensions passed to the batched solver.");

  precice::profiling::Event eInit("map.pou.gpu.initializeKokko");
  // We have to initialize Kokkos and Ginkgo here, as the initialization call allocates memory
  // in the current setup, this will only initialize the device (and allocate memory) on the primary rank
  device::Ginkgo::initialize(ginkgoParameter.nThreads, ginkgoParameter.deviceId);
  PRECICE_INFO("Using Batched solver on executor {} for {} PU-RBF clusters", ginkgoParameter.executor, centers.size());
  // _deviceExecutor = create_device_executor(ginkgoParameter.executor, ginkgoParameter.enableUnifiedMemory);
#ifdef PRECICE_WITH_OPENMP
  if (ginkgoParameter.nThreads > 0 && ginkgoParameter.executor == "omp-executor")
    omp_set_num_threads(ginkgoParameter.nThreads);
#endif
  eInit.stop();

// Assumption of the algorithm
#ifndef NDEBUG
  for (std::size_t i = 0; i < inMesh->nVertices(); ++i) {
    PRECICE_ASSERT(inMesh->vertices()[i].getID() == i);
  }
  for (std::size_t i = 0; i < outMesh->nVertices(); ++i) {
    PRECICE_ASSERT(outMesh->vertices()[i].getID() == i);
  }
#endif

  precice::profiling::Event eNearestNeighbors("map.pou.gpu.queryVertices");

  // Step 1:  Query n-nearest neighbors and compute offsets, which hold the range for each cluster

  PRECICE_DEBUG("Computing cluster association on the GPU");

  // Initialize the view for the GPU offsets
  _inOffsets  = Kokkos::View<int *>("inOffsets", _nCluster + 1);
  _outOffsets = Kokkos::View<int *>("outOffsets", _nCluster + 1);

  // we fill this view on the host side
  auto hostIn  = Kokkos::create_mirror_view(_inOffsets);
  auto hostOut = Kokkos::create_mirror_view(_outOffsets);

  // for the global IDs, we use now a std::vector which we emplace
  // we need at least contiguous memory here
  // TODO: Check the performance of reallocations
  std::vector<int> globalInIDs;
  std::vector<int> globalOutIDs;

  // has to be a separate loop, as we first need to gather knowledge about
  // the shape for the meshes
  hostIn(0)         = 0;
  hostOut(0)        = 0;
  _maxInClusterSize = _maxOutClusterSize = 0;
  for (std::size_t i = 0; i < _nCluster; ++i) {
    const auto &center = centers[i];

    // First we handle the input side
    auto inIDs        = inMesh->index().getVerticesInsideBox(center, clusterRadius);
    _maxInClusterSize = std::max(_maxInClusterSize, inIDs.size());
    hostIn(i + 1)     = hostIn(i) + inIDs.size();
    std::copy(inIDs.begin(), inIDs.end(), std::back_inserter(globalInIDs));

    // ... and the same for the output side
    auto outIDs        = outMesh->index().getVerticesInsideBox(center, clusterRadius - math::NUMERICAL_ZERO_DIFFERENCE);
    _maxOutClusterSize = std::max(_maxOutClusterSize, outIDs.size());
    hostOut(i + 1)     = hostOut(i) + outIDs.size();
    std::copy(outIDs.begin(), outIDs.end(), std::back_inserter(globalOutIDs));
  }

  // Copy offsets onto the device
  Kokkos::deep_copy(_inOffsets, hostIn);
  Kokkos::deep_copy(_outOffsets, hostOut);

  // ... now that we have the sizes, we transfer the map onto the device
  _globalInIDs  = Kokkos::View<int *>("globalInIDs", globalInIDs.size());
  _globalOutIDs = Kokkos::View<int *>("globalOutIDs", globalOutIDs.size());

  // Wrap in a view to perform deep copies further down
  Kokkos::View<int *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      tmpIn(globalInIDs.data(), globalInIDs.size());
  Kokkos::View<int *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      tmpOut(globalOutIDs.data(), globalOutIDs.size());

  Kokkos::deep_copy(_globalInIDs, tmpIn);
  Kokkos::deep_copy(_globalOutIDs, tmpOut);

  Kokkos::fence();
  eNearestNeighbors.stop();

  precice::profiling::Event eOff2d("map.pou.gpu.compute2DOffsets");

  // Step 2: Compute the matrix offsets on the device
  PRECICE_DEBUG("Computing matrix offsets");
  _kernelOffsets     = Kokkos::View<std::size_t *>("kernelOffsets", _nCluster + 1);
  _evaluationOffsets = Kokkos::View<std::size_t *>("evaluationOffsets", _nCluster + 1);

  Kokkos::deep_copy(_kernelOffsets, 0);
  Kokkos::deep_copy(_evaluationOffsets, 0);

  // We use a parallel scan for that
  kernel::compute_offsets(_inOffsets, _inOffsets, _kernelOffsets, _nCluster);
  kernel::compute_offsets(_inOffsets, _outOffsets, _evaluationOffsets, _nCluster);

  Kokkos::fence();
  eOff2d.stop();
  precice::profiling::Event eMesh("map.pou.gpu.copyMeshes");
  // Step 3: Handle the mesh data structure and copy over to the device
  PRECICE_DEBUG("Computing mesh data on the device");

  _inMesh  = Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>("inMesh", inMesh->nVertices(), _dim);
  _outMesh = Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>("outMesh", outMesh->nVertices(), _dim);

  auto hostInMesh  = Kokkos::create_mirror_view(_inMesh);
  auto hostOutMesh = Kokkos::create_mirror_view(_outMesh);

  for (std::size_t i = 0; i < inMesh->nVertices(); ++i) {
    const auto &v = inMesh->vertex(i);
    for (int d = 0; d < _dim; ++d) {
      hostInMesh(i, d) = v.rawCoords()[d];
    }
  }
  for (std::size_t i = 0; i < outMesh->nVertices(); ++i) {
    const auto &v = outMesh->vertex(i);
    for (int d = 0; d < _dim; ++d) {
      hostOutMesh(i, d) = v.rawCoords()[d];
    }
  }
  // Copy to device
  Kokkos::deep_copy(_inMesh, hostInMesh);
  Kokkos::deep_copy(_outMesh, hostOutMesh);
  Kokkos::fence();
  eMesh.stop();
  {
    precice::profiling::Event eWeights("map.pou.gpu.computeWeights");

    // Step 4: Compute the weights for each vertex
    // we first need to transfer the center coordinates and the meshes onto the device
    Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> centerMesh("centerMesh", _nCluster, _dim);
    auto                                                                        hostCenterMesh = Kokkos::create_mirror_view(centerMesh);
    for (std::size_t i = 0; i < _nCluster; ++i) {
      const auto &v = centers[i];
      for (int d = 0; d < _dim; ++d) {
        hostCenterMesh(i, d) = v.rawCoords()[d];
      }
    }
    Kokkos::deep_copy(centerMesh, hostCenterMesh);

    _normalizedWeights = Kokkos::View<double *>("normalizedWeights", globalOutIDs.size());
    CompactPolynomialC2 weightingFunction(clusterRadius);
    bool                success = kernel::compute_weights(_nCluster, globalOutIDs.size(), outMesh->nVertices(), _dim, _outOffsets, centerMesh, _globalOutIDs, _outMesh, weightingFunction, _normalizedWeights);
    PRECICE_CHECK(success, "Clustering resulted in unassigned vertices for the output mesh \"{}\".", outMesh->getName());
  }

  if (_polynomial == Polynomial::SEPARATE) {
    precice::profiling::Event ePoly("map.pou.gpu.computePolynomials");
    _qrMatrix = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>("qrMatrix", globalInIDs.size() * (_dim + 1)); // = nCluster x verticesPerCluster_i x polyParams
    _qrTau    = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>("qrTau", _nCluster * (_dim + 1));             // = nCluster x polyParams
    _qrP      = Kokkos::View<int *, Kokkos::DefaultExecutionSpace>("qrP", _nCluster * (_dim + 2));                  //  = nCluster x (polyParams + rank)
    kernel::do_batched_qr(_nCluster, _dim, _maxInClusterSize, _inOffsets, _globalInIDs, _inMesh, _qrMatrix, _qrTau, _qrP);
  }
  precice::profiling::Event eMatr("map.pou.gpu.assembleMatrices");
  // Step 6: Launch the parallel kernel to assemble the kernel matrices
  PRECICE_DEBUG("Assemble batched matrices");
  std::size_t unrolledSize   = 0;
  auto        last_elem_view = Kokkos::subview(_kernelOffsets, _nCluster);
  Kokkos::deep_copy(unrolledSize, last_elem_view);
  _kernelMatrices = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>("kernelMatrices", unrolledSize);

  kernel::do_batched_assembly(_nCluster, _dim, basisFunction, basisFunction.getFunctionParameters(),
                              _inOffsets, _globalInIDs, _inMesh, _inOffsets, _globalInIDs, _inMesh, _kernelOffsets, _kernelMatrices);

  // The eval matrices ///////////////
  std::size_t evalSize        = 0;
  auto        last_elem_view2 = Kokkos::subview(_evaluationOffsets, _nCluster);
  Kokkos::deep_copy(evalSize, last_elem_view2);
  _evalMatrices = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>("evalMatrices", evalSize);

  kernel::do_batched_assembly(_nCluster, _dim, basisFunction, basisFunction.getFunctionParameters(),
                              _inOffsets, _globalInIDs, _inMesh, _outOffsets, _globalOutIDs, _outMesh, _evaluationOffsets, _evalMatrices);

  Kokkos::fence();
  eMatr.stop();
  precice::profiling::Event eLU("map.pou.gpu.compute.lu");
  // Step 7: Compute batched lu
  PRECICE_DEBUG("Compute batched lu");
  kernel::do_batched_lu(_nCluster, _kernelOffsets, _kernelMatrices);
  Kokkos::fence();
  eLU.stop();
  precice::profiling::Event eAllo("map.pou.gpu.allocateData");
  // Step 8: Allocate memory for data transfer
  PRECICE_DEBUG("Allocate data containers for data transfer");

  _inData  = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>("inData", inMesh->nVertices());
  _outData = Kokkos::View<double *, Kokkos::DefaultExecutionSpace>("outData", outMesh->nVertices());
  // _inDataMirror  = Kokkos::create_mirror_view(_inData);
  // _outDataMirror = Kokkos::create_mirror_view(_outData);
  Kokkos::fence();
}

template <typename RADIAL_BASIS_FUNCTION_T>
void BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(const time::Sample &globalIn, Eigen::VectorXd &globalOut)
{
  PRECICE_ASSERT(globalIn.dataDims == 1, "Not implemented.");
  Kokkos::View<const double *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      inView(globalIn.values.data(), globalIn.values.size());

  // Step 1: Copy over
  precice::profiling::Event e1("map.pou.gpu.copyHostToDevice");
  Kokkos::deep_copy(_inData, inView);
  // Reset output data
  Kokkos::deep_copy(_outData, 0.0);

  Kokkos::fence();
  e1.stop();

  // Step 2: Launch kernel
  precice::profiling::Event e3("map.pou.gpu.BatchedSolve");
  if (_polynomial == Polynomial::SEPARATE) {
    kernel::do_batched_solve<true>(_nCluster, _dim, _maxInClusterSize, _maxOutClusterSize,
                                   _inOffsets, _globalInIDs, _inData, _kernelOffsets, _kernelMatrices, _normalizedWeights,
                                   _evaluationOffsets, _evalMatrices, _outOffsets, _globalOutIDs, _outData,
                                   _inMesh, _outMesh, _qrMatrix, _qrTau, _qrP);
  } else {
    kernel::do_batched_solve<false>(_nCluster, _dim, _maxInClusterSize, _maxOutClusterSize,
                                    _inOffsets, _globalInIDs, _inData, _kernelOffsets, _kernelMatrices, _normalizedWeights,
                                    _evaluationOffsets, _evalMatrices, _outOffsets, _globalOutIDs, _outData,
                                    _inMesh, _outMesh, _qrMatrix, _qrTau, _qrP);
  }

  Kokkos::fence();
  e3.stop();
  precice::profiling::Event e4("map.pou.gpu.copyDeviceToHost");
  Kokkos::View<double *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      outView(globalOut.data(), globalOut.size());

  Kokkos::deep_copy(outView, _outData);
  Kokkos::fence();
  e4.stop();
}

template <typename RADIAL_BASIS_FUNCTION_T>
void BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::clear()
{
}

} // namespace mapping
} // namespace precice

#endif // PRECICE_NO_GINKGO
