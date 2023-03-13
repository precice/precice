#pragma once
#ifndef PRECICE_NO_GINKGO

#include <array>
#include <cmath>
#include <functional>
#include <ginkgo/ginkgo.hpp>
#include <ginkgo/kernels/kernel_declaration.hpp>
#include <numeric>
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"

// Every class uses Ginkgo's default_precision = double
// Ginkgo Data Structures
using GinkgoVector = gko::matrix::Dense<>;
using GinkgoMatrix = gko::matrix::Dense<>;
using GinkgoScalar = gko::matrix::Dense<>;
// Ginkgo Solver
using cg         = gko::solver::Cg<>;
using gmres      = gko::solver::Gmres<>;
using mg         = gko::solver::Multigrid;
using ir         = gko::solver::Ir<>;
using triangular = gko::solver::UpperTrs<>;
// Ginkgo Preconditioner
using jacobi   = gko::preconditioner::Jacobi<>;
using cholesky = gko::preconditioner::Ic<>;
using ilu      = gko::preconditioner::Ilu<>;

// Ginkgo Helpers
using amgx_pgm = gko::multigrid::AmgxPgm<>; // TODO: It was later renamed to Pgm so this needs to be fixed as soon as switching to newer Ginkgo version is done

using precice::mapping::RadialBasisParameters;

// Declare Ginkgo Kernels as required by Ginkgo's unified kernel interface
GKO_DECLARE_UNIFIED(template <typename ValueType, typename EvalFunctionType> void create_rbf_system_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, const std::size_t dataDimensionality, const std::array<bool, 3> activeAxis, ValueType *mtx, ValueType *supportPoints,
    ValueType *targetPoints, EvalFunctionType f, const RadialBasisParameters rbf_params, const std::size_t inputRowLength, const std::size_t outputRowLength,
    const bool addPolynomial, const unsigned int extraDims = 0));

GKO_DECLARE_UNIFIED(template <typename ValueType> void fill_polynomial_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *x, const std::size_t supportPointsRowLength, const unsigned int dims = 4));

GKO_REGISTER_UNIFIED_OPERATION(rbf_fill_operation, create_rbf_system_matrix);
GKO_REGISTER_UNIFIED_OPERATION(polynomial_fill_operation, fill_polynomial_matrix);

extern void initCuda(const int deviceId = 0);
extern void deInitCuda();
extern void computeQR(const std::shared_ptr<gko::Executor> &exec, GinkgoMatrix *A_Q, GinkgoMatrix *R);

namespace precice {
namespace mapping {

enum class GinkgoSolverType {
  CG,
  GMRES,
  MG,
  QR
};

enum class GinkgoPreconditionerType {
  Jacobi,
  Cholesky,
  Ilu,
  None
};

// Runtime lookups as suggested by Ginkgo

const std::map<std::string, GinkgoSolverType> solverTypeLookup{
    {"cg-solver", GinkgoSolverType::CG},
    {"gmres-solver", GinkgoSolverType::GMRES},
    {"mg-solver", GinkgoSolverType::MG},
    {"qr-solver", GinkgoSolverType::QR}};

const std::map<std::string, GinkgoPreconditionerType> preconditionerTypeLookup{
    {"jacobi-preconditioner", GinkgoPreconditionerType::Jacobi},
    {"cholesky-preconditioner", GinkgoPreconditionerType::Cholesky},
    {"ilu-preconditioner", GinkgoPreconditionerType::Ilu},
    {"no-preconditioner", GinkgoPreconditionerType::None}};

const std::map<std::string, std::function<std::shared_ptr<gko::Executor>(const unsigned int, const bool)>> ginkgoExecutorLookup{{"reference-executor", [](auto unused, auto unused2) { return gko::ReferenceExecutor::create(); }},
                                                                                                                                {"omp-executor", [](auto unused, auto unused2) { return gko::OmpExecutor::create(); }},
                                                                                                                                {"cuda-executor", [](auto deviceId, auto enableUnifiedMemory) { if(enableUnifiedMemory) return gko::CudaExecutor::create(deviceId, gko::OmpExecutor::create(), true, gko::allocation_mode::unified_global); else return gko::CudaExecutor::create(deviceId, gko::OmpExecutor::create(), true, gko::allocation_mode::device); }},
                                                                                                                                {"hip-executor", [](auto deviceId, auto unused) { return gko::HipExecutor::create(0, gko::OmpExecutor::create(), true); }}};

/**
 * This class assembles and solves an RBF system, given an input mesh and an output mesh with relevant vertex IDs.
 * It uses iterative solvers (CG, GMRES) and preconditioners ((Block-)Jacobi, Cholesky, Ilu) to solve the interpolation
 * systems. Furthermore, it optionally does that on Nvidia or AMD GPUs which provides significant speedup over (single-threaded)
 * CPU implementations.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class GinkgoRadialBasisFctSolver {
public:
  GinkgoRadialBasisFctSolver() = default;

  GinkgoRadialBasisFctSolver(const MappingConfiguration::GinkgoParameter &ginkgoParameter);

  /// Assembles the system matrices and computes the decomposition of the interpolation matrix
  template <typename IndexContainer>
  GinkgoRadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                             const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial,
                             const MappingConfiguration::GinkgoParameter &ginkgoParameter);

  ~GinkgoRadialBasisFctSolver();

  GinkgoRadialBasisFctSolver(const GinkgoRadialBasisFctSolver &solver) = delete;
  GinkgoRadialBasisFctSolver &operator=(const GinkgoRadialBasisFctSolver &solver) = delete;

  /// Maps the given input data
  Eigen::VectorXd solveConsistent(const Eigen::VectorXd &rhsValues, Polynomial polynomial);

  /// Maps the given input data
  Eigen::VectorXd solveConservative(const Eigen::VectorXd &inputData, Polynomial polynomial);

  void clear();

  // Access to the evaluation matrix (output x input)
  const std::shared_ptr<GinkgoMatrix> getEvaluationMatrix() const;

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

  /// Right-hand side of the polynomial system
  std::shared_ptr<GinkgoVector> _polynomialRhs;

  /// Subtraction of the polynomial contribution
  std::shared_ptr<GinkgoVector> _subPolynomialContribution;

  /// Addition of the polynomial contribution
  std::shared_ptr<GinkgoVector> _addPolynomialContribution;

  /// Polynomial matrix of the output mesh (for separate polynomial)
  std::shared_ptr<GinkgoMatrix> _matrixV;

  /// Stores the calculated cofficients of the RBF interpolation
  std::shared_ptr<GinkgoVector> _rbfCoefficients;

  std::shared_ptr<GinkgoVector> _polynomialContribution;

  /// Matrix Q^T of QR decomposition
  std::shared_ptr<GinkgoMatrix> _decompMatrixQ_T;

  /// Matrix R of QR decomposition
  std::shared_ptr<GinkgoMatrix> _decompMatrixR;

  /// Backwards Solver
  std::shared_ptr<triangular> _triangularSolver;

  // Solver used for iteratively solving linear systems of equations
  std::shared_ptr<cg>    _cgSolver    = nullptr;
  std::shared_ptr<gmres> _gmresSolver = nullptr;
  std::shared_ptr<mg>    _mgSolver    = nullptr;

  std::shared_ptr<cg> _polynomialSolver = nullptr;

  GinkgoSolverType _solverType;

  GinkgoPreconditionerType _preconditionerType;

  // 1x1 identity matrix used for AXPY operations
  std::shared_ptr<GinkgoScalar> _scalarOne;
  std::shared_ptr<GinkgoScalar> _scalarNegativeOne;

  void _solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const;

  precice::utils::Event _allocCopyEvent{"map.rbf.ginkgo.memoryAllocAndCopy", false, false};

  precice::utils::Event _assemblyEvent{"map.rbf.ginkgo.assembleMatrices", false, false};

  std::shared_ptr<gko::log::Convergence<>> _logger;
};

template <typename RADIAL_BASIS_FUNCTION_T>
GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::GinkgoRadialBasisFctSolver(const MappingConfiguration::GinkgoParameter &ginkgoParameter)
{
  _deviceExecutor = ginkgoExecutorLookup.at(ginkgoParameter.executor)(ginkgoParameter.deviceId, ginkgoParameter.enableUnifiedMemory);
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::GinkgoRadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                                                                                const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial,
                                                                                const MappingConfiguration::GinkgoParameter &ginkgoParameter)
{
  PRECICE_INFO("Using Ginkgo solver {} on executor {} with max. iterations {} and residual reduction {}", ginkgoParameter.solver, ginkgoParameter.executor, ginkgoParameter.maxIterations, ginkgoParameter.residualNorm);
  _deviceExecutor = ginkgoExecutorLookup.at(ginkgoParameter.executor)(ginkgoParameter.deviceId, ginkgoParameter.enableUnifiedMemory);

  _solverType         = solverTypeLookup.at(ginkgoParameter.solver);
  _preconditionerType = preconditionerTypeLookup.at(ginkgoParameter.preconditioner);

  _logger = gko::share(gko::log::Convergence<>::create(_deviceExecutor, gko::log::Logger::all_events_mask));

  if (GinkgoSolverType::QR == _solverType) {
    PRECICE_ASSERT("cuda-executor" == ginkgoParameter.executor, "The QR decomposition is only available on CUDA yet.");
    initCuda(ginkgoParameter.deviceId);
  }

  PRECICE_ASSERT(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && polynomial == Polynomial::ON), "The integrated polynomial (polynomial=\"on\") is not supported for the selected radial-basis function. Please select another radial-basis function or change the polynomial configuration.");
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

  const std::size_t inputMeshSize  = inputMesh.vertices().size();
  const std::size_t outputMeshSize = outputMesh.vertices().size();
  const std::size_t meshDim        = inputMesh.vertices().at(0).getDimensions();

  _scalarOne         = gko::share(gko::initialize<GinkgoScalar>({1.0}, _deviceExecutor));
  _scalarNegativeOne = gko::share(gko::initialize<GinkgoScalar>({1.0}, _deviceExecutor));

  // Now we fill the RBF system matrix on the GPU (or any other selected device)
  _allocCopyEvent.start();
  _rbfCoefficients = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{n, 1}));
  _allocCopyEvent.pause();
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
        inputVertices->at(j, i) = inputMesh.vertices().at(i).rawCoords()[j];
      } else {
        inputVertices->at(i, j) = inputMesh.vertices().at(i).rawCoords()[j];
      }
    }
  }
  for (std::size_t i = 0; i < outputMeshSize; ++i) {
    for (std::size_t j = 0; j < meshDim; ++j) {
      if ("cuda-executor" == ginkgoParameter.executor || "hip-executor" == ginkgoParameter.executor) {
        outputVertices->at(j, i) = outputMesh.vertices().at(i).rawCoords()[j];
      } else {
        outputVertices->at(i, j) = outputMesh.vertices().at(i).rawCoords()[j];
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

  _allocCopyEvent.pause();

  if (polynomial == Polynomial::SEPARATE) {
    const unsigned int separatePolyParams = 4 - std::count(activeAxis.begin(), activeAxis.end(), false);
    _allocCopyEvent.start();
    _matrixQ = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{n, separatePolyParams}));
    _matrixV = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{outputSize, separatePolyParams}));
    _allocCopyEvent.pause();

    _assemblyEvent.start();
    _deviceExecutor->run(make_polynomial_fill_operation(_matrixQ->get_size()[0], _matrixQ->get_size()[1], _matrixQ->get_values(), dInputVertices->get_values(), dInputVertices->get_size()[1], separatePolyParams));
    _deviceExecutor->run(make_polynomial_fill_operation(_matrixV->get_size()[0], _matrixV->get_size()[1], _matrixV->get_values(), dOutputVertices->get_values(), dOutputVertices->get_size()[1], separatePolyParams));
    _assemblyEvent.pause();

    _deviceExecutor->synchronize();

    _matrixQ_T = gko::share(_matrixQ->transpose());

    _allocCopyEvent.start();
    _matrixQ_TQ                = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>{_matrixQ_T->get_size()[0], _matrixQ->get_size()[1]}));
    _polynomialRhs             = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ_T->get_size()[0], 1}));
    _subPolynomialContribution = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ->get_size()[0], 1}));
    _addPolynomialContribution = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixV->get_size()[0], 1}));
    _allocCopyEvent.pause();

    _matrixQ_T->apply(gko::lend(_matrixQ), gko::lend(_matrixQ_TQ));

    auto polynomialSolverFactory = cg::build()
                                       .with_criteria(gko::stop::Iteration::build()
                                                          .with_max_iters(static_cast<std::size_t>(1e6))
                                                          .on(_deviceExecutor),
                                                      gko::stop::ResidualNormReduction<>::build()
                                                          .with_reduction_factor(1e-4)
                                                          .on(_deviceExecutor))
                                       .on(_deviceExecutor);

    _polynomialSolver = polynomialSolverFactory->generate(_matrixQ_TQ);
  }

  // Launch RBF fill kernel on device
  _assemblyEvent.start();
  precice::utils::Event systemMatrixAssemblyEvent{"map.rbf.ginkgo.assembleSystemMatrix", false};
  _deviceExecutor->run(make_rbf_fill_operation(_rbfSystemMatrix->get_size()[0], _rbfSystemMatrix->get_size()[1], meshDim, activeAxis, _rbfSystemMatrix->get_values(), dInputVertices->get_values(), dInputVertices->get_values(), basisFunction, basisFunction.getFunctionParameters(), dInputVertices->get_size()[1], dInputVertices->get_size()[1], Polynomial::ON == polynomial, polyparams)); // polynomial evaluates to true only if ON is set
  _deviceExecutor->synchronize();
  systemMatrixAssemblyEvent.stop();

  precice::utils::Event outputMatrixAssemblyEvent{"map.rbf.ginkgo.assembleOutputMatrix", false};
  _deviceExecutor->run(make_rbf_fill_operation(_matrixA->get_size()[0], _matrixA->get_size()[1], meshDim, activeAxis, _matrixA->get_values(), dInputVertices->get_values(), dOutputVertices->get_values(), basisFunction, basisFunction.getFunctionParameters(), dInputVertices->get_size()[1], dOutputVertices->get_size()[1], Polynomial::ON == polynomial, polyparams));

  // Wait for the kernels to finish
  _deviceExecutor->synchronize();
  outputMatrixAssemblyEvent.stop();
  _assemblyEvent.stop();

  dInputVertices->clear();
  dOutputVertices->clear();

  auto iterationCriterion = gko::share(gko::stop::Iteration::build()
                                           .with_max_iters(ginkgoParameter.maxIterations)
                                           .on(_deviceExecutor));

  auto residualCriterion = gko::share(gko::stop::ResidualNormReduction<>::build()
                                          .with_reduction_factor(ginkgoParameter.residualNorm)
                                          .on(_deviceExecutor));

  iterationCriterion->add_logger(_logger);
  residualCriterion->add_logger(_logger);

  if (_solverType == GinkgoSolverType::MG) {

    auto smootherFactory = gko::share(
        ir::build()
            .with_solver(jacobi::build().with_max_block_size(1u).on(_deviceExecutor))
            .with_relaxation_factor(0.9)
            .with_criteria(
                gko::stop::Iteration::build().with_max_iters(2u).on(_deviceExecutor))
            .on(_deviceExecutor));

    auto mgLevelFactory = amgx_pgm::build().with_deterministic(false).on(_deviceExecutor);

    auto coarsestFactory =
        cg::build()
            .with_criteria(
                gko::stop::Iteration::build().with_max_iters(4u).on(_deviceExecutor))
            .on(_deviceExecutor);

    auto multigridFactory =
        mg::build()
            .with_max_levels(2u)
            .with_min_coarse_rows(2u)
            .with_pre_smoother(gko::share(smootherFactory))
            .with_post_uses_pre(true)
            .with_mg_level(gko::share(mgLevelFactory))
            .with_coarsest_solver(
                gko::share(jacobi::build().with_max_block_size(1u).on(_deviceExecutor)))
            .with_criteria(residualCriterion)
            .on(_deviceExecutor);

    _mgSolver = gko::share(multigridFactory->generate(_rbfSystemMatrix));
    _mgSolver->add_logger(_logger);

  } else if (_solverType == GinkgoSolverType::CG) {

    if (GinkgoPreconditionerType::None != _preconditionerType && ginkgoParameter.usePreconditioner) {
      auto solverFactoryWithPreconditioner = [preconditionerType = _preconditionerType, executor = _deviceExecutor, &ginkgoParameter]() {
        if (preconditionerType == GinkgoPreconditionerType::Jacobi) {
          return cg::build().with_preconditioner(jacobi::build().with_max_block_size(ginkgoParameter.jacobiBlockSize).on(executor));
        } else if (preconditionerType == GinkgoPreconditionerType::Cholesky) {
          return cg::build().with_preconditioner(cholesky::build().on(executor));
        } else {
          return cg::build().with_preconditioner(ilu::build().on(executor));
        }
      }();

      auto solverFactory = solverFactoryWithPreconditioner
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(_deviceExecutor);

      _cgSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
      _cgSolver->add_logger(_logger);
    }

    else {
      auto solverFactory = cg::build()
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(_deviceExecutor);

      _cgSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
      _cgSolver->add_logger(_logger);
    }

  } else if (_solverType == GinkgoSolverType::GMRES) {

    if (GinkgoPreconditionerType::None != _preconditionerType && ginkgoParameter.usePreconditioner) {
      auto solverFactoryWithPreconditioner = [preconditionerType = _preconditionerType, executor = _deviceExecutor, &ginkgoParameter]() {
        if (preconditionerType == GinkgoPreconditionerType::Jacobi) {
          return gmres::build().with_preconditioner(jacobi::build().with_max_block_size(ginkgoParameter.jacobiBlockSize).on(executor));
        } else if (preconditionerType == GinkgoPreconditionerType::Cholesky) {
          return gmres::build().with_preconditioner(cholesky::build().on(executor));
        } else {
          return gmres::build().with_preconditioner(ilu::build().on(executor));
        }
      }();

      auto solverFactory = solverFactoryWithPreconditioner
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(_deviceExecutor);

      _gmresSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
      _gmresSolver->add_logger(_logger);
    } else {
      auto solverFactory = gmres::build()
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(_deviceExecutor);

      _gmresSolver = gko::share(solverFactory->generate(_rbfSystemMatrix));
      _gmresSolver->add_logger(_logger);
    }
  } else if (_solverType == GinkgoSolverType::QR) {
    const std::size_t M = _rbfSystemMatrix->get_size()[0];
    const std::size_t N = _rbfSystemMatrix->get_size()[1];
    _decompMatrixQ_T    = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>(N, M)));
    _decompMatrixR      = gko::share(GinkgoMatrix::create(_deviceExecutor, gko::dim<2>(N, N)));

    // _rbfSystemMatrix will be overridden into Q
    computeQR(_deviceExecutor, gko::lend(_rbfSystemMatrix), gko::lend(_decompMatrixR));

    _rbfSystemMatrix->transpose(gko::lend(_decompMatrixQ_T));

    auto triangularSolverFactory = triangular::build().on(_deviceExecutor);
    _triangularSolver            = gko::share(triangularSolverFactory->generate(_decompMatrixR));
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::_solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const
{
  precice::utils::Event solverEvent("map.rbf.ginkgo.solveSystemMatrix");
  if (_solverType == GinkgoSolverType::CG) {
    _cgSolver->apply(gko::lend(rhs), gko::lend(_rbfCoefficients));
  } else if (_solverType == GinkgoSolverType::GMRES) {
    _gmresSolver->apply(gko::lend(rhs), gko::lend(_rbfCoefficients));
  } else if (_solverType == GinkgoSolverType::MG) {
    _mgSolver->apply(gko::lend(rhs), gko::lend(_rbfCoefficients));
  }
  solverEvent.stop();
  PRECICE_INFO("The iterative solver stopped after {} iterations.", _logger->get_num_iterations());

// Only compute time-consuming statistics in debug mode
#ifndef NDEBUG

  auto dResidual = gko::initialize<GinkgoScalar>({0.0}, _deviceExecutor);
  _rbfSystemMatrix->apply(gko::lend(_scalarOne), gko::lend(_rbfCoefficients), gko::lend(_scalarNegativeOne), gko::lend(rhs));
  rhs->compute_norm2(gko::lend(dResidual));
  auto residual = gko::clone(_hostExecutor, dResidual);
  PRECICE_INFO("Ginkgo Solver Final Residual: {}", residual->at(0, 0));

#endif
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(const Eigen::VectorXd &rhsValues, Polynomial polynomial)
{
  PRECICE_ASSERT(rhsValues.cols() == 1);
  // Copy rhs vector onto GPU by creating a Ginkgo Vector
  auto rhs = gko::share(GinkgoVector::create(_hostExecutor, gko::dim<2>{static_cast<unsigned long>(rhsValues.rows()), 1}));

  for (Eigen::Index i = 0; i < rhsValues.rows(); ++i) {
    rhs->at(i, 0) = rhsValues(i, 0);
  }

  _allocCopyEvent.start();
  auto dRhs = gko::share(gko::clone(_deviceExecutor, rhs));
  rhs->clear();
  _allocCopyEvent.pause();

  if (polynomial == Polynomial::SEPARATE) {
    _allocCopyEvent.start();
    _polynomialContribution = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixQ_TQ->get_size()[1], 1}));
    _allocCopyEvent.pause();
    _polynomialContribution->fill(0.0);

    _matrixQ_T->apply(gko::lend(dRhs), gko::lend(_polynomialRhs));
    _polynomialSolver->apply(gko::lend(_polynomialRhs), gko::lend(_polynomialContribution));

    _matrixQ->apply(gko::lend(_polynomialContribution), gko::lend(_subPolynomialContribution));
    dRhs->sub_scaled(gko::lend(_scalarOne), gko::lend(_subPolynomialContribution));
  }

  if (GinkgoSolverType::QR == _solverType) {
    // New rhs Q^T * b
    auto dQ_T_Rhs = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_decompMatrixQ_T->get_size()[0], 1}));
    _decompMatrixQ_T->apply(gko::lend(dRhs), gko::lend(dQ_T_Rhs));
    _triangularSolver->apply(gko::lend(dQ_T_Rhs), gko::lend(_rbfCoefficients));
  } else {
    _solveRBFSystem(dRhs);
  }

  dRhs->clear();

  _allocCopyEvent.start();
  auto dOutput = gko::share(GinkgoVector::create(_deviceExecutor, gko::dim<2>{_matrixA->get_size()[0], _rbfCoefficients->get_size()[1]}));
  _allocCopyEvent.pause();

  _matrixA->apply(gko::lend(_rbfCoefficients), gko::lend(dOutput));

  if (polynomial == Polynomial::SEPARATE) {
    _matrixV->apply(gko::lend(_polynomialContribution), gko::lend(_addPolynomialContribution));
    dOutput->add_scaled(gko::lend(_scalarOne), gko::lend(_addPolynomialContribution));
  }

  _allocCopyEvent.start();
  auto output = gko::clone(_hostExecutor, dOutput);
  _allocCopyEvent.pause();

  Eigen::VectorXd result(output->get_size()[0], 1);

  for (Eigen::Index i = 0; i < result.rows(); ++i) {
    result(i, 0) = output->at(i, 0);
  }

  return result;
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConservative(const Eigen::VectorXd &rhsValues, Polynomial polynomial)
{
  return Eigen::VectorXd(1, 1);
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::shared_ptr<gko::Executor> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getReferenceExecutor() const
{
  return _hostExecutor;
}

template <typename RADIAL_BASIS_FUNCTION_T>
const std::shared_ptr<GinkgoMatrix> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getEvaluationMatrix() const
{
  return _matrixA;
}

template <typename RADIAL_BASIS_FUNCTION_T>
GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::~GinkgoRadialBasisFctSolver()
{
  _allocCopyEvent.stop();

  if (GinkgoSolverType::QR == _solverType) {
    deInitCuda();
  }

  clear();
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
