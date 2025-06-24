#pragma once
#ifndef PRECICE_NO_KOKKOS_KERNELS

#include <Kokkos_Core.hpp>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/device/Device.hpp"
#include "mapping/device/KokkosPUMKernels.hpp"
#include "mapping/device/KokkosTypes.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/impl/SphericalVertexCluster.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"

using precice::mapping::RadialBasisParameters;

namespace precice::mapping {

/**
 * This class solves the PU-RBF interpolation in a batched manner, i.e.,
 * clusters are computed in parallel. The parallelization can be handled
 * through GPUs or OpenMP. It uses Kokkos-kernels dispatch the kernels on
 * the selected execution backend. Only the mesh indexing is handled on the
 * CPU through boost geometry, all other components are handled in parallel.
 * The solver is compatible with the MPI-parallel layout of the PURBF class
 * itself, which means in this case that the problem is solved purely local,
 * but each rank can instantiate and launch its own solver.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class BatchedRBFSolver {
public:
  using RBF_T = RADIAL_BASIS_FUNCTION_T;

  /// Essentially an initialize of the solver: allocates device memory
  /// and computes the LU decompositions etc
  BatchedRBFSolver(RBF_T                                 basisFunction,
                   mesh::PtrMesh                         inMesh,
                   mesh::PtrMesh                         outMesh,
                   const std::vector<mesh::Vertex>      &centers,
                   double                                clusterRadius,
                   Polynomial                            polynomial,
                   bool                                  computeEvaluationOffline,
                   MappingConfiguration::GinkgoParameter ginkgoParameter);

  void solveConsistent(const time::Sample &globalIn, Eigen::VectorXd &globalOut);

private:
  mutable precice::logging::Logger _log{"mapping::BatchedRBFSolver"};

  // Helper to dispatch the actual kernel. Needed to help with the template parameters
  template <typename... Args>
  void _dispatch_solve_kernel(bool polynomial, bool evaluation_op_available, Args &&...args);

  // Linear offsets for each cluster, i.e., all cluster sizes
  VectorOffsetView<> _inOffsets;
  VectorOffsetView<> _outOffsets;

  // Stores for each cluster the VertexIDs
  GlobalIDView<> _globalInIDs;
  GlobalIDView<> _globalOutIDs;

  MatrixOffsetView<> _kernelOffsets;
  MatrixOffsetView<> _evaluationOffsets;

  MeshView<> _inMesh;
  MeshView<> _outMesh;

  VectorView<> _qrMatrix; // flat view of (nCluster x verticesPerCluster_i x (dim + 1) = nCluster x verticesPerCluster_i x polyParams)
  VectorView<> _qrTau;    // flat view of Householder tau (nCluster x (dim + 1) = nCluster x polyParams)
  PivotView<>  _qrP;      // flat view of Permutation and rank (nCluster x (dim + 2) = nCluster x (polyParams + rank))

  VectorView<> _kernelMatrices;

  VectorView<> _evalMatrices;
  VectorView<> _normalizedWeights;

  // Currently only scalar data
  VectorView<> _inData;
  // Kokkos::View<double *>::HostMirror                    _inDataMirror;
  VectorView<> _outData;
  // Kokkos::View<double *>::HostMirror                    _outDataMirror;

  int _maxInClusterSize;
  int _maxOutClusterSize;

  RBF_T      _basisFunction;
  Polynomial _polynomial;
  const int  _nCluster;
  const int  _dim; // Mesh dimension
  int        _avgClusterSize{};
  const bool _computeEvaluationOffline;
  // MappingConfiguration::GinkgoParameter _ginkgoParameter;
};

template <typename RADIAL_BASIS_FUNCTION_T>
BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::BatchedRBFSolver(RBF_T                                 basisFunction,
                                                            mesh::PtrMesh                         inMesh,
                                                            mesh::PtrMesh                         outMesh,
                                                            const std::vector<mesh::Vertex>      &centers,
                                                            double                                clusterRadius,
                                                            Polynomial                            polynomial,
                                                            bool                                  computeEvaluationOffline,
                                                            MappingConfiguration::GinkgoParameter ginkgoParameter)
    : _basisFunction(basisFunction), _polynomial(polynomial), _nCluster(static_cast<int>(centers.size())), _dim(inMesh->getDimensions()), _computeEvaluationOffline(computeEvaluationOffline)
{
  PRECICE_TRACE();
  PRECICE_CHECK(_polynomial != Polynomial::ON, "Setting polynomial to \"on\" for the mapping between \"{}\" and \"{}\" is not supported", inMesh->getName(), outMesh->getName());
  // The LU decomposition uses no pivoting, which leads to divisions by zero if the diagonal contains zero entries, which is the case for our basis functions
  PRECICE_CHECK(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite(), "batched solver is only available for positive definite basis functions, i.e., compact-polynomial functions and Gaussian.");

  PRECICE_CHECK(!(inMesh->vertices().empty() || outMesh->vertices().empty()), "One of the meshes in the batched solvers is empty, which is invalid.");
  PRECICE_CHECK(inMesh->getDimensions() == outMesh->getDimensions(), "Incompatible dimensions passed to the batched solver.");

  precice::profiling::Event eInit("solver.initializeKokkos");
  // We have to initialize Kokkos and Ginkgo here, as the initialization call allocates memory
  // in the current setup, this will only initialize the device (and allocate memory) on the primary rank
  // TODO: Document restriction: all mappings must use the same executor configuration within one participant
  device::Device::initialize(ginkgoParameter.nThreads, ginkgoParameter.deviceId);
  PRECICE_DEBUG("Using batched PU-RBF solver on executor \"{}\" for \"{}\" PU-RBF clusters in execution mode {}.", ginkgoParameter.executor, centers.size(), _computeEvaluationOffline ? "\"minimal-compute\" (evaluation offline)" : "\"minimal-memory\" (evaluation online)");
  Kokkos::fence();
  eInit.stop();

// General assumption of the algorithm
#ifndef NDEBUG
  for (int i = 0; i < inMesh->nVertices(); ++i) {
    PRECICE_ASSERT(inMesh->vertices()[i].getID() == i);
  }
  for (int i = 0; i < outMesh->nVertices(); ++i) {
    PRECICE_ASSERT(outMesh->vertices()[i].getID() == i);
  }
#endif

  precice::profiling::Event eNearestNeighbors("solver.queryVertices");

  // Step 1:  Query n-nearest neighbors and compute offsets, which hold the range for each cluster

  PRECICE_DEBUG("Computing cluster association on the GPU");

  // Initialize the view for the GPU offsets
  _inOffsets  = VectorOffsetView<>("inOffsets", _nCluster + 1);
  _outOffsets = VectorOffsetView<>("outOffsets", _nCluster + 1);

  // we fill this view on the host side
  auto hostIn  = Kokkos::create_mirror_view(_inOffsets);
  auto hostOut = Kokkos::create_mirror_view(_outOffsets);

  // for the global IDs, we use now a std::vector which we emplace
  // we need at least contiguous memory here
  // TODO: Check the performance of reallocations
  std::vector<VertexID> globalInIDs;
  std::vector<VertexID> globalOutIDs;

  // has to be a separate loop, as we first need to gather knowledge about
  // the shape for the meshes
  hostIn(0)  = 0;
  hostOut(0) = 0;
  // To detect overflows
  std::uint64_t inCheck  = 0;
  std::uint64_t outCheck = 0;
  _maxInClusterSize = _maxOutClusterSize = 0;
  for (int i = 0; i < _nCluster; ++i) {
    const auto &center = centers[i];

    // First we handle the input side
    auto inIDs          = inMesh->index().getVerticesInsideBox(center, clusterRadius);
    _maxInClusterSize   = std::max(_maxInClusterSize, static_cast<int>(inIDs.size()));
    std::uint64_t tmpIn = hostIn(i) + inIDs.size();

    // Check overflows
    if constexpr (std::numeric_limits<offset_1d_type>::digits < std::numeric_limits<std::uint64_t>::digits) {
      PRECICE_CHECK(tmpIn < std::numeric_limits<offset_1d_type>::max(),
                    "The selected integer precision for the (input) vector offsets (\"offset_1d_type\") overflow. You might want to change the precision specified in \"device/KokkosTypes.hpp\"");
    }
    if constexpr (std::numeric_limits<offset_2d_type>::digits < std::numeric_limits<std::uint64_t>::digits) {
      inCheck += static_cast<std::uint64_t>(inIDs.size() * inIDs.size());
      PRECICE_CHECK(inCheck < std::numeric_limits<offset_2d_type>::max(),
                    "The selected integer precision for the (input) matrix offsets (\"offset_2d_type\") overflow. You might want to change the precision specified in \"device/KokkosTypes.hpp\"");
    }
    hostIn(i + 1) = static_cast<offset_1d_type>(tmpIn);
    std::copy(inIDs.begin(), inIDs.end(), std::back_inserter(globalInIDs));

    // ... and the same for the output side
    auto outIDs          = outMesh->index().getVerticesInsideBox(center, clusterRadius - math::NUMERICAL_ZERO_DIFFERENCE);
    _maxOutClusterSize   = std::max(_maxOutClusterSize, static_cast<int>(outIDs.size()));
    std::uint64_t tmpOut = hostOut(i) + outIDs.size();

    // Check overflows
    if constexpr (std::numeric_limits<offset_1d_type>::digits < std::numeric_limits<std::uint64_t>::digits) {
      PRECICE_CHECK(tmpOut < std::numeric_limits<offset_1d_type>::max(),
                    "The selected integer precision for the (output) vector offsets (\"offset_1d_type\") overflow. You might want to change the precision specified in \"device/KokkosTypes.hpp\"");
    }
    if (_computeEvaluationOffline) {
      if constexpr (std::numeric_limits<offset_2d_type>::digits < std::numeric_limits<std::uint64_t>::digits) {
        outCheck += static_cast<std::uint64_t>(outIDs.size() * inIDs.size()); // is in x out
        PRECICE_CHECK(outCheck < std::numeric_limits<offset_2d_type>::max(),
                      "The selected integer precision for the (output) matrix offsets (\"offset_2d_type\") overflow. You might want to change the precision specified in \"device/KokkosTypes.hpp\"");
      }
    }
    hostOut(i + 1) = static_cast<offset_1d_type>(tmpOut);
    std::copy(outIDs.begin(), outIDs.end(), std::back_inserter(globalOutIDs));
  }

  _avgClusterSize = hostIn(_nCluster /* = hostIn.extent(0) - 1 */) / _nCluster;
  PRECICE_DEBUG("Average cluster size used to find a good team size of the kernel execution: {}", _avgClusterSize);

  // Copy offsets onto the device
  Kokkos::deep_copy(_inOffsets, hostIn);
  Kokkos::deep_copy(_outOffsets, hostOut);

  // ... now that we have the sizes, we transfer the map onto the device
  _globalInIDs  = GlobalIDView<>("globalInIDs", globalInIDs.size());
  _globalOutIDs = GlobalIDView<>("globalOutIDs", globalOutIDs.size());

  // Wrap in a view to perform deep copies further down
  Kokkos::View<VertexID *, Kokkos::HostSpace, UnmanagedMemory>
      tmpIn(globalInIDs.data(), globalInIDs.size());
  Kokkos::View<VertexID *, Kokkos::HostSpace, UnmanagedMemory>
      tmpOut(globalOutIDs.data(), globalOutIDs.size());

  Kokkos::deep_copy(_globalInIDs, tmpIn);
  Kokkos::deep_copy(_globalOutIDs, tmpOut);

  Kokkos::fence();
  eNearestNeighbors.stop();

  precice::profiling::Event eOff2d("solver.kernel.compute2DOffsets");

  // Step 2: Compute the matrix offsets on the device
  PRECICE_DEBUG("Computing matrix offsets");
  // We use a parallel scan for that
  _kernelOffsets = MatrixOffsetView<>("kernelOffsets", _nCluster + 1);
  Kokkos::deep_copy(_kernelOffsets, 0);
  kernel::compute_offsets(_inOffsets, _inOffsets, _kernelOffsets, _nCluster);

  if (_computeEvaluationOffline) {
    _evaluationOffsets = MatrixOffsetView<>("evaluationOffsets", _nCluster + 1);
    Kokkos::deep_copy(_evaluationOffsets, 0);
    kernel::compute_offsets(_inOffsets, _outOffsets, _evaluationOffsets, _nCluster);
  }
  Kokkos::fence();
  eOff2d.stop();
  precice::profiling::Event eMesh("solver.copyMeshes");
  // Step 3: Handle the mesh data structure and copy over to the device
  PRECICE_DEBUG("Computing mesh data on the device");

  _inMesh  = MeshView<>("inMesh", inMesh->nVertices(), _dim);
  _outMesh = MeshView<>("outMesh", outMesh->nVertices(), _dim);

  auto hostInMesh  = Kokkos::create_mirror_view(_inMesh);
  auto hostOutMesh = Kokkos::create_mirror_view(_outMesh);

  for (int i = 0; i < inMesh->nVertices(); ++i) {
    const auto &v = inMesh->vertex(i);
    for (int d = 0; d < _dim; ++d) {
      hostInMesh(i, d) = v.rawCoords()[d];
    }
  }
  for (int i = 0; i < outMesh->nVertices(); ++i) {
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
    PRECICE_DEBUG("Computing PU-RBF weights");
    precice::profiling::Event eWeights("solver.kernel.computeWeights");

    // Step 4: Compute the weights for each vertex
    // we first need to transfer the center coordinates and the meshes onto the device
    MeshView<> centerMesh("centerMesh", _nCluster, _dim);
    auto       hostCenterMesh = Kokkos::create_mirror_view(centerMesh);
    for (int i = 0; i < _nCluster; ++i) {
      const auto &v = centers[i];
      for (int d = 0; d < _dim; ++d) {
        hostCenterMesh(i, d) = v.rawCoords()[d];
      }
    }
    Kokkos::deep_copy(centerMesh, hostCenterMesh);

    _normalizedWeights = VectorView<>("normalizedWeights", globalOutIDs.size());
    CompactPolynomialC2 weightingFunction(clusterRadius);
    // Computing the weights parallelizes over the number of output mesh vertices
    int  avgOutClusterSize = hostOut(_nCluster /* = hostOut.extent(0) - 1 */) / _nCluster;
    bool success           = kernel::compute_weights(_nCluster, avgOutClusterSize, globalOutIDs.size(), outMesh->nVertices(), _dim, _outOffsets,
                                                     centerMesh, _globalOutIDs, _outMesh, weightingFunction, _normalizedWeights);
    PRECICE_CHECK(success, "Clustering resulted in unassigned vertices for the output mesh \"{}\".", outMesh->getName());
    Kokkos::fence();
  }

  PRECICE_ASSERT(_avgClusterSize > 0);
  if (_polynomial == Polynomial::SEPARATE) {
    PRECICE_DEBUG("Computing polynomial QR");
    precice::profiling::Event ePoly("solver.kernel.computePolynomialQR");
    _qrMatrix = VectorView<>("qrMatrix", globalInIDs.size() * (_dim + 1)); // = nCluster x verticesPerCluster_i x polyParams
    _qrTau    = VectorView<>("qrTau", _nCluster * (_dim + 1));             // = nCluster x polyParams
    _qrP      = PivotView<>("qrP", _nCluster * (_dim + 2));                //  = nCluster x (polyParams + rank)
    kernel::do_batched_qr(_nCluster, _dim, _avgClusterSize, _maxInClusterSize, _inOffsets, _globalInIDs, _inMesh, _qrMatrix, _qrTau, _qrP);
    Kokkos::fence();
  }
  precice::profiling::Event eMatr("solver.kernel.assembleInputMatrices");
  // Step 6: Launch the parallel kernel to assemble the kernel matrices
  PRECICE_DEBUG("Assemble batched matrices");
  // The kernel matrices /////////////
  offset_2d_type unrolledSize   = 0;
  auto           last_elem_view = Kokkos::subview(_kernelOffsets, _nCluster);
  Kokkos::deep_copy(unrolledSize, last_elem_view);
  _kernelMatrices = VectorView<>("kernelMatrices", unrolledSize);

  kernel::do_input_assembly(_nCluster, _dim, _avgClusterSize, _maxInClusterSize, basisFunction,
                            _inOffsets, _globalInIDs, _inMesh, _kernelOffsets, _kernelMatrices);

  Kokkos::fence();
  eMatr.stop();

  if (_computeEvaluationOffline) {
    // The eval matrices ///////////////
    precice::profiling::Event eMatrOut("solver.kernel.assembleOutputMatrices");
    offset_2d_type            evalSize        = 0;
    auto                      last_elem_view2 = Kokkos::subview(_evaluationOffsets, _nCluster);
    Kokkos::deep_copy(evalSize, last_elem_view2);
    _evalMatrices = VectorView<>("evalMatrices", evalSize);

    kernel::do_batched_assembly(_nCluster, _dim, _avgClusterSize, basisFunction,
                                _inOffsets, _globalInIDs, _inMesh, _outOffsets, _globalOutIDs, _outMesh, _evaluationOffsets, _evalMatrices);
    Kokkos::fence();
  }

  precice::profiling::Event eLU("solver.kernel.lu");
  // Step 7: Compute batched lu
  PRECICE_DEBUG("Compute batched lu");
  kernel::do_batched_lu(_nCluster, _avgClusterSize, _kernelOffsets, _kernelMatrices);
  Kokkos::fence();
  eLU.stop();
  precice::profiling::Event eAllo("solver.allocateData");
  // Step 8: Allocate memory for data transfer
  PRECICE_DEBUG("Allocate data containers for data transfer");

  _inData  = VectorView<>("inData", inMesh->nVertices());
  _outData = VectorView<>("outData", outMesh->nVertices());
  Kokkos::fence();
}

template <typename RADIAL_BASIS_FUNCTION_T>
void BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(const time::Sample &globalIn, Eigen::VectorXd &globalOut)
{
  auto solve_component =
      [&](const double *inPtr, Eigen::Index inSize, double *outPtr, Eigen::Index outSize) {
        // Step 1: Wrap memory into an unmanaged view
        Kokkos::View<const double *, Kokkos::HostSpace, UnmanagedMemory>
            inView(inPtr, inSize);

        // Step 2: Copy over
        precice::profiling::Event e1("solver.copyHostToDevice");
        Kokkos::deep_copy(_inData, inView);
        Kokkos::deep_copy(_outData, 0.0); // Reset output data

        Kokkos::fence();
        e1.stop();

        // Step 3: Launch the kernel
        precice::profiling::Event e2("solver.kernel.batchedSolve");
        _dispatch_solve_kernel(_polynomial == Polynomial::SEPARATE, _computeEvaluationOffline,
                               _nCluster, _dim, _avgClusterSize, _maxInClusterSize, _maxOutClusterSize, _basisFunction,
                               _inOffsets, _globalInIDs, _inData, _kernelOffsets, _kernelMatrices, _normalizedWeights,
                               _evaluationOffsets, _evalMatrices, _outOffsets, _globalOutIDs, _outData,
                               _inMesh, _outMesh, _qrMatrix, _qrTau, _qrP);

        Kokkos::fence();
        e2.stop();

        // Step 4: Copy back
        precice::profiling::Event e3("solver.copyDeviceToHost");
        Kokkos::View<double *, Kokkos::HostSpace, UnmanagedMemory>
            outView(outPtr, outSize);
        Kokkos::deep_copy(outView, _outData);
        Kokkos::fence();
        e3.stop();
      };

  const int nComponents = globalIn.dataDims;

  // If we have just one component, we can directly copy the data over and solve
  if (nComponents == 1) {
    solve_component(globalIn.values.data(), globalIn.values.size(), globalOut.data(), globalOut.size());
  } else {
    // Otherwise, we map the data to a component-wise matrix
    Eigen::Map<const Eigen::MatrixXd> inMatrix(globalIn.values.data(), nComponents, globalIn.values.size() / nComponents);
    Eigen::Map<Eigen::MatrixXd>       outMatrix(globalOut.data(), nComponents, globalOut.size() / nComponents);

    // ... and solve component-wise. This requires each component to be contiguous in memory
    Eigen::VectorXd tmpIn(inMatrix.cols());
    Eigen::VectorXd tmpOut(outMatrix.cols());
    for (int c = 0; c < nComponents; ++c) {
      tmpIn = inMatrix.row(c);
      solve_component(tmpIn.data(), tmpIn.size(), tmpOut.data(), tmpOut.size());
      outMatrix.row(c) = tmpOut;
    }
  }
}

// Forwarding dispatcher for the actual implementation to help with the template parameters:
template <typename RADIAL_BASIS_FUNCTION_T>
template <typename... Args>
void BatchedRBFSolver<RADIAL_BASIS_FUNCTION_T>::_dispatch_solve_kernel(bool polynomial, bool evaluation_op_available,
                                                                       Args &&...args)
{
  if (polynomial) {
    if (evaluation_op_available)
      kernel::do_batched_solve<true, true>(std::forward<Args>(args)...);
    else
      kernel::do_batched_solve<true, false>(std::forward<Args>(args)...);
  } else {
    if (evaluation_op_available)
      kernel::do_batched_solve<false, true>(std::forward<Args>(args)...);
    else
      kernel::do_batched_solve<false, false>(std::forward<Args>(args)...);
  }
}
} // namespace precice::mapping

#else

/// Stub class in case we compile without the feature
#include "mapping/config/MappingConfiguration.hpp"
namespace precice::mapping {

template <typename RADIAL_BASIS_FUNCTION_T>
class BatchedRBFSolver {
public:
  BatchedRBFSolver(RADIAL_BASIS_FUNCTION_T,
                   mesh::PtrMesh,
                   mesh::PtrMesh,
                   const std::vector<mesh::Vertex> &,
                   double,
                   Polynomial,
                   bool,
                   MappingConfiguration::GinkgoParameter) {}

  void solveConsistent(const time::Sample &, Eigen::VectorXd &) {}
};
} // namespace precice::mapping
#endif // PRECICE_NO_KOKKOS_KERNELS
