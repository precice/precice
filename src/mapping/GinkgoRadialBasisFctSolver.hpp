#pragma once
#ifndef PRECICE_NO_GINKGO

#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include "mapping/GinkgoDefinitions.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/device/GinkgoRBFKernels.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "utils/Ginkgo.hpp"
#ifdef PRECICE_WITH_HIP
#include "mapping/device/HipQRSolver.hip.hpp"
#endif
#ifdef PRECICE_WITH_CUDA
#include "mapping/device/CudaQRSolver.cuh"
#endif
#ifdef PRECICE_WITH_OMP
#include <omp.h>
#endif

using precice::mapping::RadialBasisParameters;

namespace precice {
namespace mapping {

enum class GinkgoSolverType {
  CG,
  GMRES,
  QR
};

enum class GinkgoPreconditionerType {
  Jacobi,
  Cholesky,
  None
};

// Runtime lookups as suggested by Ginkgo

const std::map<std::string, GinkgoSolverType> solverTypeLookup{
    {"cg-solver", GinkgoSolverType::CG},
    {"gmres-solver", GinkgoSolverType::GMRES},
    {"qr-solver", GinkgoSolverType::QR}};

const std::map<std::string, GinkgoPreconditionerType> preconditionerTypeLookup{
    {"jacobi-preconditioner", GinkgoPreconditionerType::Jacobi},
    {"cholesky-preconditioner", GinkgoPreconditionerType::Cholesky},
    {"no-preconditioner", GinkgoPreconditionerType::None}};

/**
 * This class assembles and solves an RBF system, given an input mesh and an output mesh with relevant vertex IDs.
 * It uses iterative solvers (CG, GMRES) and preconditioners ((Block-)Jacobi, Cholesky, Ilu) to solve the interpolation
 * systems. Furthermore, it optionally does that on Nvidia or AMD GPUs which provides significant speedup over (single-threaded)
 * CPU implementations.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class GinkgoRadialBasisFctSolver {
public:
  using BASIS_FUNCTION_T       = RADIAL_BASIS_FUNCTION_T;
  GinkgoRadialBasisFctSolver() = default;

  /// Assembles the system matrices and computes the decomposition of the interpolation matrix
  template <typename IndexContainer>
  GinkgoRadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                             const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial,
                             MappingConfiguration::GinkgoParameter ginkgoParameter);

  /// Maps the given input data
  Eigen::VectorXd solveConsistent(const Eigen::VectorXd &inputData, Polynomial polynomial);

  /// Maps the given input data
  Eigen::VectorXd solveConservative(const Eigen::VectorXd &inputData, Polynomial polynomial);

  void clear();

  Eigen::Index getInputSize() const;

  Eigen::Index getOutputSize() const;

  std::shared_ptr<gko::Executor> getReferenceExecutor() const;

private:
  mutable precice::logging::Logger _log{"mapping::GinkgoRadialBasisFctSolver"};

  std::shared_ptr<gko::Executor> _deviceExecutor;
  std::shared_ptr<gko::Executor> _hostExecutor = gko::ReferenceExecutor::create();

  // Stores the RBF interpolation matrix
  std::shared_ptr<GinkgoMatrix> _rbfSystemMatrix;

  /// Evaluation matrix (output x input)
  std::shared_ptr<GinkgoMatrix> _matrixA;

  /// Polynomial matrix of the input mesh (for separate polynomial)
  std::shared_ptr<GinkgoMatrix> _matrixQ;

  /// Transposed Polynomial matrix of the input mesh (for separate polynomial) (to solve Q^T*Q*x=Q^T*b)
  std::shared_ptr<gko::LinOp> _matrixQ_T;

  /// Product Q^T*Q (to solve Q^TQx=Q^Tb)
  std::shared_ptr<gko::LinOp> _matrixQ_TQ;

  /// Product Q*Q^T
  std::shared_ptr<gko::LinOp> _matrixQQ_T;

  /// Right-hand side of the polynomial system
  std::shared_ptr<GinkgoVector> _polynomialRhs;

  /// Subtraction of the polynomial contribution
  std::shared_ptr<GinkgoVector> _subPolynomialContribution;

  /// Addition of the polynomial contribution
  std::shared_ptr<GinkgoVector> _addPolynomialContribution;

  /// Polynomial matrix of the output mesh (for separate polynomial)
  std::shared_ptr<GinkgoMatrix> _matrixV;

  /// Stores the calculated coefficients of the RBF interpolation
  std::shared_ptr<GinkgoVector> _rbfCoefficients;

  std::shared_ptr<GinkgoVector> _polynomialContribution;

  /// Matrix Q^T of QR decomposition
  std::shared_ptr<GinkgoMatrix> _decompMatrixQ_T;

  /// Q^T * b of QR decomposition
  std::shared_ptr<GinkgoMatrix> _dQ_T_Rhs;

  /// Matrix R of QR decomposition
  std::shared_ptr<GinkgoMatrix> _decompMatrixR;

  /// Backwards Solver
  std::shared_ptr<triangular> _triangularSolver;

  /// QR Solver
  // std::unique_ptr<QRSolver> _qrSolver;

  // Solver used for iteratively solving linear systems of equations
  std::shared_ptr<cg>    _cgSolver    = nullptr;
  std::shared_ptr<gmres> _gmresSolver = nullptr;

  std::shared_ptr<cg> _polynomialSolver = nullptr;

  GinkgoSolverType _solverType = GinkgoSolverType::CG;

  GinkgoPreconditionerType _preconditionerType;

  // 1x1 identity matrix used for AXPY operations
  std::shared_ptr<GinkgoScalar> _scalarOne;
  std::shared_ptr<GinkgoScalar> _scalarNegativeOne;

  void _solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const;

  std::shared_ptr<gko::stop::Iteration::Factory> _iterationCriterion;

  std::shared_ptr<gko::stop::ResidualNorm<>::Factory> _residualCriterion;

  std::shared_ptr<gko::stop::ResidualNorm<>::Factory> _absoluteResidualCriterion;

  MappingConfiguration::GinkgoParameter _ginkgoParameter;
};

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::GinkgoRadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                                                                                const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial,
                                                                                MappingConfiguration::GinkgoParameter ginkgoParameter)
    : _ginkgoParameter(ginkgoParameter)
{
  PRECICE_TRACE();
  PRECICE_INFO("Using Ginkgo solver {} on executor {} with max. iterations {} and residual reduction {}",
               ginkgoParameter.solver,
               ginkgoParameter.executor,
               ginkgoParameter.maxIterations,
               ginkgoParameter.residualNorm);
  _deviceExecutor = create_device_executor(ginkgoParameter.enableUnifiedMemory);
#ifdef PRECICE_WITH_OMP
  if (_ginkgoParameter.nThreads > 0 && _ginkgoParameter.executor == "omp-executor")
    omp_set_num_threads(_ginkgoParameter.nThreads);
#endif
  _solverType         = solverTypeLookup.at(ginkgoParameter.solver);
  _preconditionerType = preconditionerTypeLookup.at(ginkgoParameter.preconditioner);

  PRECICE_CHECK(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && polynomial == Polynomial::ON), "The integrated polynomial (polynomial=\"on\") is not supported for the selected radial-basis function. Please select another radial-basis function or change the polynomial configuration.");
  // Convert dead axis vector into an active axis array so that we can handle the reduction more easily
  std::array<bool, 3> activeAxis({{false, false, false}});
  std::transform(deadAxis.begin(), deadAxis.end(), activeAxis.begin(), [](const auto ax) { return !ax; });

  const std::size_t deadDimensions = std::count(activeAxis.begin(), activeAxis.end(), false);
  const std::size_t dimensions     = 3;
  const std::size_t polyparams     = polynomial == Polynomial::ON ? 1 + dimensions - deadDimensions : 0;

  // Add linear polynom degrees if polynomial requires this
  const auto inputSize  = inputIDs.size();
  const auto outputSize = outputIDs.size();
  const auto n          = inputSize + polyparams;

  PRECICE_ASSERT((inputMesh.getDimensions() == 3) || activeAxis[2] == false);
  PRECICE_ASSERT((inputSize >= 1 + polyparams) || polynomial != Polynomial::ON, inputSize);

  const std::size_t inputMeshSize  = inputMesh.nVertices();
  const std::size_t outputMeshSize = outputMesh.nVertices();
  const std::size_t meshDim        = inputMesh.vertex(0).getDimensions();

  _scalarOne         = gko::share(gko::initialize<GinkgoScalar>({1.0}, _deviceExecutor));
  _scalarNegativeOne = gko::share(gko::initialize<GinkgoScalar>({-1.0}, _deviceExecutor));

  // Now we fill the RBF system matrix on the GPU (or any other selected device)
  precice::profiling::Event _allocCopyEvent{"map.rbf.ginkgo.memoryAllocAndCopy"};
  _rbfCoefficients = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{n, 1}));
  _allocCopyEvent.stop();
  // Initial guess is required since uninitialized memory could lead to a never converging system
  _rbfCoefficients->fill(0.0);

  // We need to copy the input data into a CPU stored vector first and copy it to the GPU afterwards
  // To allow for coalesced memory accesses on the GPU, we need to store them in transposed order IFF the backend is the GPU
  // However, the CPU does not need that; in fact, it would make it slower
  std::size_t inputVerticesM, inputVerticesN, outputVerticesM, outputVerticesN;

  if ("cuda-executor" == ginkgoParameter.executor || "hip-executor" == ginkgoParameter.executor) {
    inputVerticesM  = meshDim;
    inputVerticesN  = inputMeshSize;
    outputVerticesM = meshDim;
    outputVerticesN = outputMeshSize;
  } else {
    inputVerticesM  = inputMeshSize;
    inputVerticesN  = meshDim;
    outputVerticesM = outputMeshSize;
    outputVerticesN = meshDim;
  }

  auto inputVertices  = gko::share(GinkgoMatrix::create(_hostExecutor, gko::dim<2>{inputVerticesM, inputVerticesN}));
  auto outputVertices = gko::share(GinkgoMatrix::create(_hostExecutor, gko::dim<2>{outputVerticesM, outputVerticesN}));
  for (std::size_t i = 0; i < inputMeshSize; ++i) {
    for (std::size_t j = 0; j < meshDim; ++j) {
      if ("cuda-executor" == ginkgoParameter.executor || "hip-executor" == ginkgoParameter.executor) {
        inputVertices->at(j, i) = inputMesh.vertex(i).coord(j);
      } else {
        inputVertices->at(i, j) = inputMesh.vertex(i).coord(j);
      }
    }
  }
  for (std::size_t i = 0; i < outputMeshSize; ++i) {
    for (std::size_t j = 0; j < meshDim; ++j) {
      if ("cuda-executor" == ginkgoParameter.executor || "hip-executor" == ginkgoParameter.executor) {
        outputVertices->at(j, i) = outputMesh.vertex(i).coord(j);
      } else {
        outputVertices->at(i, j) = outputMesh.vertex(i).coord(j);
      }
    }
  }

  _allocCopyEvent.start();

  auto dInputVertices  = gko::clone(_deviceExecutor, inputVertices);
  auto dOutputVertices = gko::clone(_deviceExecutor, outputVertices);
  inputVertices->clear();
  outputVertices->clear();

  _deviceExecutor->synchronize();

  _rbfSystemMatrix = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{n, n}));
  _matrixA         = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{outputSize, n}));

  _allocCopyEvent.stop();

  if (polynomial == Polynomial::SEPARATE) {
    const unsigned int separatePolyParams = 4 - std::count(activeAxis.begin(), activeAxis.end(), false);
    _allocCopyEvent.start();
    _matrixQ = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{n, separatePolyParams}));
    _matrixV = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{outputSize, separatePolyParams}));
    _allocCopyEvent.stop();

    _matrixQ->fill(0.0);
    _matrixV->fill(0.0);

    precice::profiling::Event _assemblyEvent{"map.rbf.ginkgo.assembleMatrices"};
    kernel::fill_polynomial_matrix(_deviceExecutor, _matrixQ, dInputVertices, separatePolyParams);
    kernel::fill_polynomial_matrix(_deviceExecutor, _matrixV, dOutputVertices, separatePolyParams);
    _assemblyEvent.stop();

    _deviceExecutor->synchronize();

    _matrixQ_T = gko::share(_matrixQ->transpose());

    _allocCopyEvent.start();
    _matrixQ_TQ                = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{_matrixQ_T->get_size()[0], _matrixQ->get_size()[1]}));
    _polynomialRhs             = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ_T->get_size()[0], 1}));
    _subPolynomialContribution = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ->get_size()[0], 1}));
    _addPolynomialContribution = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixV->get_size()[0], 1}));
    _allocCopyEvent.stop();

    _matrixQ_T->apply(_matrixQ, _matrixQ_TQ);

    auto polynomialSolverFactory = cg::build()
                                       .with_criteria(gko::stop::Iteration::build()
                                                          .with_max_iters(static_cast<std::size_t>(40))
                                                          .on(_deviceExecutor),
                                                      gko::stop::ResidualNorm<>::build()
                                                          .with_reduction_factor(1e-6)
                                                          .with_baseline(gko::stop::mode::initial_resnorm)
                                                          .on(_deviceExecutor))
                                       .on(_deviceExecutor);

    _polynomialSolver = polynomialSolverFactory->generate(_matrixQ_TQ);
  }

  // Launch RBF fill kernel on device
  precice::profiling::Event _assemblyEvent{"map.rbf.ginkgo.assembleMatrices"};
  precice::profiling::Event systemMatrixAssemblyEvent{"map.rbf.ginkgo.assembleSystemMatrix"};
  kernel::create_rbf_system_matrix(_deviceExecutor, _rbfSystemMatrix, activeAxis, dInputVertices, dInputVertices, basisFunction,
                                   basisFunction.getFunctionParameters(), Polynomial::ON == polynomial,
                                   polyparams); // polynomial evaluates to true only if ON is set
  _deviceExecutor->synchronize();
  systemMatrixAssemblyEvent.stop();

  precice::profiling::Event outputMatrixAssemblyEvent{"map.rbf.ginkgo.assembleOutputMatrix"};
  kernel::create_rbf_system_matrix(_deviceExecutor, _matrixA, activeAxis, dInputVertices, dOutputVertices, basisFunction,
                                   basisFunction.getFunctionParameters(), Polynomial::ON == polynomial, polyparams);

  // Wait for the kernels to finish
  _deviceExecutor->synchronize();
  outputMatrixAssemblyEvent.stop();
  _assemblyEvent.stop();

  dInputVertices->clear();
  dOutputVertices->clear();

  _iterationCriterion = gko::share(gko::stop::Iteration::build()
                                       .with_max_iters(ginkgoParameter.maxIterations)
                                       .on(_deviceExecutor));

  _residualCriterion = gko::share(gko::stop::ResidualNorm<>::build()
                                      .with_reduction_factor(ginkgoParameter.residualNorm)
                                      .with_baseline(gko::stop::mode::initial_resnorm)
                                      .on(_deviceExecutor));

  // For cases where we reach a stationary solution such that the coupling data doesn't change (or map zero data)
  _absoluteResidualCriterion = gko::share(gko::stop::ResidualNorm<>::build()
                                              .with_reduction_factor(1e-30)
                                              .with_baseline(gko::stop::mode::absolute)
                                              .on(_deviceExecutor));

  if (_solverType == GinkgoSolverType::CG) {

    if (GinkgoPreconditionerType::None != _preconditionerType && ginkgoParameter.usePreconditioner) {
      auto solverFactoryWithPreconditioner = [preconditionerType = _preconditionerType, executor = _deviceExecutor, &ginkgoParameter]() {
        if (preconditionerType == GinkgoPreconditionerType::Jacobi) {
          return cg::build().with_preconditioner(jacobi::build().with_max_block_size(ginkgoParameter.jacobiBlockSize).on(executor));
        } else {
          return cg::build().with_preconditioner(cholesky::build().on(executor));
        }
      }();

      auto solverFactory = solverFactoryWithPreconditioner
                               .with_criteria(_iterationCriterion, _residualCriterion, _absoluteResidualCriterion)
                               .on(_deviceExecutor);

      _cgSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
    } else {
      auto solverFactory = cg::build()
                               .with_criteria(_iterationCriterion, _residualCriterion, _absoluteResidualCriterion)
                               .on(_deviceExecutor);

      _cgSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
    }

  } else if (_solverType == GinkgoSolverType::GMRES) {

    if (GinkgoPreconditionerType::None != _preconditionerType && ginkgoParameter.usePreconditioner) {
      auto solverFactoryWithPreconditioner = [preconditionerType = _preconditionerType, executor = _deviceExecutor, &ginkgoParameter]() {
        if (preconditionerType == GinkgoPreconditionerType::Jacobi) {
          return gmres::build().with_preconditioner(jacobi::build().with_max_block_size(ginkgoParameter.jacobiBlockSize).on(executor));
        } else {
          return gmres::build().with_preconditioner(cholesky::build().on(executor));
        }
      }();

      auto solverFactory = solverFactoryWithPreconditioner
                               .with_criteria(_iterationCriterion, _residualCriterion, _absoluteResidualCriterion)
                               .on(_deviceExecutor);

      _gmresSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
    } else {
      auto solverFactory = gmres::build()
                               .with_criteria(_iterationCriterion, _residualCriterion, _absoluteResidualCriterion)
                               .on(_deviceExecutor);

      _gmresSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
    }
  } else if (_solverType == GinkgoSolverType::QR) {
    const std::size_t M = _rbfSystemMatrix->get_size()[0];
    const std::size_t N = _rbfSystemMatrix->get_size()[1];
    _decompMatrixQ_T    = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>(N, M)));
    _decompMatrixR      = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>(N, N)));

    if ("cuda-executor" == ginkgoParameter.executor) {
#ifdef PRECICE_WITH_CUDA
      // _rbfSystemMatrix will be overridden into Q
      computeQRDecompositionCuda(ginkgoParameter.deviceId, _deviceExecutor, _rbfSystemMatrix, _decompMatrixR);
#endif
    } else if ("hip-executor" == ginkgoParameter.executor) {
#ifdef PRECICE_WITH_HIP
      // _rbfSystemMatrix will be overridden into Q
      computeQRDecompositionHip(ginkgoParameter.deviceId, _deviceExecutor, _rbfSystemMatrix.get(), _decompMatrixR.get());
#endif
    } else {
      PRECICE_UNREACHABLE("Not implemented");
    }
    _rbfSystemMatrix->transpose(_decompMatrixQ_T);

    _dQ_T_Rhs = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_decompMatrixQ_T->get_size()[0], 1}));

    auto triangularSolverFactory = triangular::build().on(_deviceExecutor);
    _triangularSolver            = gko::share(triangularSolverFactory->generate(_decompMatrixR));
  } else {
    PRECICE_UNREACHABLE("Unknown solver type");
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::_solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const
{
  PRECICE_TRACE();
  auto logger = gko::share(gko::log::Convergence<>::create(gko::log::Logger::all_events_mask));

  _iterationCriterion->add_logger(logger);
  _residualCriterion->add_logger(logger);
  _absoluteResidualCriterion->add_logger(logger);

  precice::profiling::Event solverEvent("map.rbf.ginkgo.solveSystemMatrix");
  if (_solverType == GinkgoSolverType::CG) {
    _cgSolver->apply(rhs, _rbfCoefficients);
  } else if (_solverType == GinkgoSolverType::GMRES) {
    _gmresSolver->apply(rhs, _rbfCoefficients);
  }
  solverEvent.stop();
  PRECICE_INFO("The iterative solver stopped after {} iterations.", logger->get_num_iterations());

// Only compute time-consuming statistics in debug mode
#ifndef NDEBUG
  auto dResidual = gko::initialize<GinkgoScalar>({0.0}, _deviceExecutor);
  _rbfSystemMatrix->apply(_scalarOne, _rbfCoefficients, _scalarNegativeOne, rhs);
  rhs->compute_norm2(dResidual);
  auto residual = gko::clone(_hostExecutor, dResidual);
  PRECICE_INFO("Ginkgo Solver Final Residual: {}", residual->at(0, 0));
#endif

  _iterationCriterion->clear_loggers();
  _residualCriterion->clear_loggers();
  _absoluteResidualCriterion->clear_loggers();
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(const Eigen::VectorXd &rhsValues, Polynomial polynomial)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(rhsValues.cols() == 1);
  // Copy rhs vector onto GPU by creating a Ginkgo Vector
  auto rhs = gko::share(GinkgoVector::create(_hostExecutor, gko::dim<2>{static_cast<unsigned long>(rhsValues.rows()), 1}));

  for (Eigen::Index i = 0; i < rhsValues.rows(); ++i) {
    rhs->at(i, 0) = rhsValues(i, 0);
  }

  precice::profiling::Event _allocCopyEvent{"map.rbf.ginkgo.memoryAllocAndCopy"};
  auto                      dRhs = gko::share(gko::clone(_deviceExecutor, rhs));
  rhs->clear();
  _allocCopyEvent.stop();

  if (polynomial == Polynomial::SEPARATE) {
    _allocCopyEvent.start();
    _polynomialContribution = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ_TQ->get_size()[1], 1}));
    _allocCopyEvent.stop();
    _polynomialContribution->fill(0.0);

    _matrixQ_T->apply(dRhs, _polynomialRhs);
    _polynomialSolver->apply(_polynomialRhs, _polynomialContribution);

    _matrixQ->apply(_polynomialContribution, _subPolynomialContribution);
    dRhs->sub_scaled(_scalarOne, _subPolynomialContribution);
  }

  if (GinkgoSolverType::QR == _solverType) {
    _decompMatrixQ_T->apply(dRhs, _dQ_T_Rhs);
    _triangularSolver->apply(_dQ_T_Rhs, _rbfCoefficients);
  } else {
    _solveRBFSystem(dRhs);
  }

  dRhs->clear();

  _allocCopyEvent.start();
  auto dOutput = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixA->get_size()[0], _rbfCoefficients->get_size()[1]}));
  _allocCopyEvent.stop();

  _matrixA->apply(_rbfCoefficients, dOutput);

  if (polynomial == Polynomial::SEPARATE) {
    _matrixV->apply(_polynomialContribution, _addPolynomialContribution);
    dOutput->add_scaled(_scalarOne, _addPolynomialContribution);
  }

  _allocCopyEvent.start();
  auto output = gko::clone(_hostExecutor, dOutput);
  _allocCopyEvent.stop();

  Eigen::VectorXd result(output->get_size()[0], 1);

  for (Eigen::Index i = 0; i < result.rows(); ++i) {
    result(i, 0) = output->at(i, 0);
  }

  return result;
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConservative(const Eigen::VectorXd &rhsValues, Polynomial polynomial)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(rhsValues.cols() == 1);
  // Copy rhs vector onto GPU by creating a Ginkgo Vector
  auto rhs = gko::share(GinkgoVector::create(_hostExecutor, gko::dim<2>{static_cast<unsigned long>(rhsValues.rows()), 1}));

  for (Eigen::Index i = 0; i < rhsValues.rows(); ++i) {
    rhs->at(i, 0) = rhsValues(i, 0);
  }

  precice::profiling::Event _allocCopyEvent{"map.rbf.ginkgo.memoryAllocAndCopy"};
  auto                      dRhs = gko::share(gko::clone(_deviceExecutor, rhs));
  rhs->clear();
  _allocCopyEvent.stop();

  auto dAu = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixA->get_size()[1], dRhs->get_size()[1]}));

  _matrixA->transpose()->apply(dRhs, dAu);

  if (GinkgoSolverType::QR == _solverType) {
    _decompMatrixQ_T->apply(dAu, _dQ_T_Rhs);
    _triangularSolver->apply(_dQ_T_Rhs, _rbfCoefficients);
  } else {
    _solveRBFSystem(dAu);
  }

  auto dOutput = gko::clone(_deviceExecutor, _rbfCoefficients);

  if (polynomial == Polynomial::SEPARATE) {
    auto dEpsilon = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixV->get_size()[1], dRhs->get_size()[1]}));
    _matrixV->transpose()->apply(dRhs, dEpsilon);

    auto dTmp = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ->get_size()[1], _rbfCoefficients->get_size()[1]}));
    _matrixQ->transpose()->apply(dOutput, dTmp);

    // epsilon -= tmp
    dEpsilon->sub_scaled(_scalarOne, dTmp);

    // Since this class is constructed for consistent mapping per default, we have to delete unused memory and initialize conservative variables
    if (nullptr == _matrixQQ_T) {
      _matrixQ_TQ->clear();
      _deviceExecutor->synchronize();
      _matrixQQ_T = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{_matrixQ->get_size()[0], _matrixQ_T->get_size()[1]}));

      _matrixQ->apply(_matrixQ_T, _matrixQQ_T);

      auto polynomialSolverFactory = cg::build()
                                         .with_criteria(gko::stop::Iteration::build()
                                                            .with_max_iters(static_cast<std::size_t>(40))
                                                            .on(_deviceExecutor),
                                                        gko::stop::ResidualNorm<>::build()
                                                            .with_reduction_factor(1e-6)
                                                            .with_baseline(gko::stop::mode::initial_resnorm)
                                                            .on(_deviceExecutor))
                                         .on(_deviceExecutor);

      _polynomialSolver = polynomialSolverFactory->generate(_matrixQQ_T);

      _polynomialRhs->clear();
      _deviceExecutor->synchronize();
    }

    _polynomialContribution = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQQ_T->get_size()[1], 1}));
    _polynomialContribution->fill(0.0);

    dEpsilon->scale(_scalarNegativeOne);

    _polynomialRhs = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ->get_size()[0], dEpsilon->get_size()[1]}));

    _matrixQ->apply(dEpsilon, _polynomialRhs);

    _polynomialSolver->apply(_polynomialRhs, _polynomialContribution);

    // out -= poly
    dOutput->sub_scaled(_scalarOne, _polynomialContribution);
  }

  _allocCopyEvent.start();
  auto output = gko::clone(_hostExecutor, dOutput);
  _allocCopyEvent.stop();

  Eigen::VectorXd result(output->get_size()[0], 1);

  for (Eigen::Index i = 0; i < result.rows(); ++i) {
    result(i, 0) = output->at(i, 0);
  }

  return result;
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::shared_ptr<gko::Executor> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getReferenceExecutor() const
{
  return _hostExecutor;
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::Index GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getInputSize() const
{
  return _matrixA->get_size()[1];
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::Index GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getOutputSize() const
{
  return _matrixA->get_size()[0];
}

template <typename RADIAL_BASIS_FUNCTION_T>
void GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::clear()
{
  if (nullptr != _rbfSystemMatrix) {
    _rbfSystemMatrix->clear();
  }
  if (nullptr != _matrixA) {
    _matrixA->clear();
  }
  if (nullptr != _matrixV) {
    _matrixV->clear();
  }
  if (nullptr != _matrixQ) {
    _matrixQ->clear();
  }
  if (nullptr != _matrixQ_T) {
    _matrixQ_T->clear();
  }
  if (nullptr != _matrixQ_TQ) {
    _matrixQ_TQ->clear();
  }
  if (nullptr != _rbfCoefficients) {
    _rbfCoefficients->clear();
  }
  if (nullptr != _polynomialRhs) {
    _polynomialRhs->clear();
  }
  if (nullptr != _subPolynomialContribution) {
    _subPolynomialContribution->clear();
  }
  if (nullptr != _addPolynomialContribution) {
    _addPolynomialContribution->clear();
  }
  if (nullptr != _polynomialContribution) {
    _polynomialContribution->clear();
  }
}

} // namespace mapping
} // namespace precice

#endif // PRECICE_NO_GINKGO
