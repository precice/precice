#pragma once
#ifndef PRECICE_NO_GINKGO

#include <array>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/irange.hpp>
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

// Declare Ginkgo Kernels
GKO_DECLARE_UNIFIED(template <typename ValueType, typename EvalFunctionType> void create_rbf_system_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, const std::size_t dataDimensionality, const std::array<bool, 3> activeAxis, ValueType *mtx, ValueType *supportPoints,
    ValueType *targetPoints, EvalFunctionType f, const std::array<ValueType, 3> rbf_params, const std::size_t inputRowLength, const std::size_t outputRowLength,
    const bool addPolynomial, const unsigned int extraDims = 0));

GKO_DECLARE_UNIFIED(template <typename ValueType> void fill_polynomial_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *x, const std::size_t supportPointsRowLength, const unsigned int dims = 4));

GKO_REGISTER_UNIFIED_OPERATION(rbf_fill_operation, create_rbf_system_matrix);
GKO_REGISTER_UNIFIED_OPERATION(polynomial_fill_operation, fill_polynomial_matrix);

namespace precice {
namespace mapping {

// Every class uses Ginkgo's default_precision = double
// Ginkgo Data Structures
using GinkgoVector = gko::matrix::Dense<>;
using GinkgoMatrix = gko::matrix::Dense<>;
using GinkgoScalar = gko::matrix::Dense<>;
// Ginkgo Solver
using cg    = gko::solver::Cg<>;
using gmres = gko::solver::Gmres<>;
using mg    = gko::solver::Multigrid;
using ir    = gko::solver::Ir<>;
// Ginkgo Preconditioner
using jacobi   = gko::preconditioner::Jacobi<>;
using cholesky = gko::preconditioner::Ic<>;
using ilu      = gko::preconditioner::Ilu<>;

// Ginkgo Helpers
using amgx_pgm = gko::multigrid::AmgxPgm<>; // TODO: It was later renamed to Pgm so this needs to be fixed as soon as switching to newer Ginkgo version is done

enum SolverType {
  CG,
  GMRES,
  MG
};

enum PreconditionerType {
  Jacobi,
  Cholesky,
  Ilu,
  None
};

const std::map<std::string, SolverType> solverTypeLookup{
    {"cg-solver", SolverType::CG},
    {"gmres-solver", SolverType::GMRES},
    {"mg-solver", SolverType::MG}};

const std::map<std::string, PreconditionerType> preconditionerTypeLookup{
    {"jacobi-preconditioner", PreconditionerType::Jacobi},
    {"cholesky-preconditioner", PreconditionerType::Cholesky},
    {"ilu-preconditioner", PreconditionerType::Ilu},
    {"no-preconditioner", PreconditionerType::None}};

const std::map<std::string, std::function<std::shared_ptr<gko::Executor>()>> ginkgoExecutorLookup{{"reference-executor", [] { return gko::ReferenceExecutor::create(); }},
                                                                                                  {"omp-executor", [] { return gko::OmpExecutor::create(); }},
                                                                                                  {"cuda-executor", [] { return gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true, gko::allocation_mode::device); }},
                                                                                                  {"hip-executor", [] { return gko::HipExecutor::create(0, gko::OmpExecutor::create(), true); }}};

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

  ~GinkgoRadialBasisFctSolver() = default;

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

  std::shared_ptr<gko::Executor>        _deviceExecutor;
  static std::shared_ptr<gko::Executor> _hostExecutor;
  static std::shared_ptr<gko::Executor> _ompExecutor;

  // Stores the RBF interpolation matrix
  std::shared_ptr<GinkgoMatrix> _rbfSystemMatrix;

  /// Decomposition of the interpolation matrix
  std::shared_ptr<GinkgoMatrix> _decMatrixC;

  /// Evaluation matrix (output x input)
  std::shared_ptr<GinkgoMatrix> _matrixA;

  /// Decomposition of the polynomial (for separate polynomial)
  // Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qrMatrixQ;

  /// Polynomial matrix of the input mesh (for separate polynomial)
  std::shared_ptr<GinkgoMatrix> _matrixQ;

  /// Polynomial matrix of the output mesh (for separate polynomial)
  std::shared_ptr<GinkgoMatrix> _matrixV;

  /// @brief Stores the calculated cofficients of the RBF interpolation
  std::shared_ptr<GinkgoVector> _rbfCoefficients;

  std::shared_ptr<GinkgoVector> _polynomialContribution;

  // Solver used for iteratively solving linear systems of equations TODO: Find out how to make dynamic for different solver families
  std::shared_ptr<precice::mapping::cg>    _cgSolver;
  std::shared_ptr<precice::mapping::gmres> _gmresSolver;
  std::shared_ptr<precice::mapping::mg>    _mgSolver;

  SolverType _solverType;

  PreconditionerType _preconditionerType;

  // 1x1 identity matrix used for AXPY operations
  std::shared_ptr<GinkgoScalar> _scalarOne;
  std::shared_ptr<GinkgoScalar> _scalarNegativeOne;

  void _solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const;

  precice::utils::Event _copyEvent{"map.rbf.ginkgo.memCopy", false, false};

  precice::utils::Event _assemblyEvent{"map.rbf.ginkgo.assembleMatrices", false, false};

  std::shared_ptr<gko::log::Convergence<>> _logger;
};

template <typename RADIAL_BASIS_FUNCTION_T>
GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::GinkgoRadialBasisFctSolver(const MappingConfiguration::GinkgoParameter &ginkgoParameter)
{
  this->_deviceExecutor = ginkgoExecutorLookup.at(ginkgoParameter.executor)();
}

template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::GinkgoRadialBasisFctSolver(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const IndexContainer &inputIDs,
                                                                                const mesh::Mesh &outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial,
                                                                                const MappingConfiguration::GinkgoParameter &ginkgoParameter)
{
  this->_deviceExecutor = ginkgoExecutorLookup.at(ginkgoParameter.executor)();

  this->_solverType         = solverTypeLookup.at(ginkgoParameter.solver);
  this->_preconditionerType = preconditionerTypeLookup.at(ginkgoParameter.preconditioner);

  this->_logger = gko::share(gko::log::Convergence<>::create(this->_deviceExecutor, gko::log::Logger::all_events_mask));

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

  this->_scalarOne         = gko::share(gko::initialize<GinkgoScalar>({1.0}, this->_deviceExecutor));
  this->_scalarNegativeOne = gko::share(gko::initialize<GinkgoScalar>({1.0}, this->_deviceExecutor));

  // Now we fill the RBF system matrix on the GPU (or any other selected device)
  this->_rbfCoefficients = gko::share(GinkgoVector::create(this->_hostExecutor, gko::dim<2>{n, 1}));

  // We need to copy the input data into a CPU stored vector first and copy it to the GPU afterwards
  // To allow for coalesced memory accesses on the GPU, we need to store them in transposed order IFF the backend is the GPU
  // However, the CPU does not need that; in fact, it would make it slower
  std::size_t inputVerticesM, inputVerticesN, outputVerticesM, outputVerticesN;

  if ("cuda-executor" == ginkgoParameter.executor) {
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

  auto inputVertices  = gko::share(GinkgoMatrix::create(this->_hostExecutor, gko::dim<2>{inputVerticesM, inputVerticesN}));
  auto outputVertices = gko::share(GinkgoMatrix::create(this->_hostExecutor, gko::dim<2>{outputVerticesM, outputVerticesN}));
  for (std::size_t i = 0; i < inputMeshSize; ++i) {
    for (std::size_t j = 0; j < meshDim; ++j) {
      if ("cuda-executor" == ginkgoParameter.executor) {
        inputVertices->at(j, i) = inputMesh.vertices().at(i).rawCoords()[j];
      } else {
        inputVertices->at(i, j) = inputMesh.vertices().at(i).rawCoords()[j];
      }
    }
    // Initial guess is required since memory chunk could lead to never converging system
    if (i < n) {
      this->_rbfCoefficients->at(i, 0) = 0.0;
    }
  }
  for (std::size_t i = 0; i < outputMeshSize; ++i) {
    for (std::size_t j = 0; j < meshDim; ++j) {
      if ("cuda-executor" == ginkgoParameter.executor) {
        outputVertices->at(j, i) = outputMesh.vertices().at(i).rawCoords()[j];
      } else {
        outputVertices->at(i, j) = outputMesh.vertices().at(i).rawCoords()[j];
      }
    }
  }

  this->_copyEvent.start();
  //  TODO: Check how to circumvent if executor is OMP -> CPU side
  inputVertices          = gko::clone(this->_deviceExecutor, inputVertices);
  outputVertices         = gko::clone(this->_deviceExecutor, outputVertices);
  this->_rbfCoefficients = gko::clone(this->_deviceExecutor, this->_rbfCoefficients); // TODO: Check if Ginkgo supports zero vector creation
  this->_deviceExecutor->synchronize();
  this->_copyEvent.pause();

  this->_rbfSystemMatrix = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{n, n}));
  this->_matrixA         = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{outputSize, n}));

  if (polynomial == Polynomial::SEPARATE) {
    const unsigned int polyParams = 4 - std::count(activeAxis.begin(), activeAxis.end(), false);
    this->_matrixQ                = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{n, polyParams}));
    this->_matrixV                = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{outputSize, polyParams}));

    this->_assemblyEvent.start();
    this->_deviceExecutor->run(make_polynomial_fill_operation(this->_matrixQ->get_size()[0], this->_matrixQ->get_size()[1], this->_matrixQ->get_values(), inputVertices->get_values(), inputVertices->get_size()[1], polyParams));
    this->_deviceExecutor->run(make_polynomial_fill_operation(this->_matrixV->get_size()[0], this->_matrixV->get_size()[1], this->_matrixV->get_values(), outputVertices->get_values(), outputVertices->get_size()[1], polyParams));
    this->_assemblyEvent.pause();

    this->_deviceExecutor->synchronize();
  }

  // Launch RBF fill kernel on device
  this->_assemblyEvent.start();
  this->_deviceExecutor->run(make_rbf_fill_operation(this->_rbfSystemMatrix->get_size()[0], this->_rbfSystemMatrix->get_size()[1], meshDim, activeAxis, this->_rbfSystemMatrix->get_values(), inputVertices->get_values(), inputVertices->get_values(), basisFunction, basisFunction.getFunctionParameters(), inputVertices->get_size()[1], inputVertices->get_size()[1], Polynomial::ON == polynomial, polyparams)); // polynomial evaluates to true only if ON is set
  this->_deviceExecutor->run(make_rbf_fill_operation(this->_matrixA->get_size()[0], this->_matrixA->get_size()[1], meshDim, activeAxis, this->_matrixA->get_values(), inputVertices->get_values(), outputVertices->get_values(), basisFunction, basisFunction.getFunctionParameters(), inputVertices->get_size()[1], outputVertices->get_size()[1], Polynomial::ON == polynomial, polyparams));

  // Wait for the kernels to finish
  this->_deviceExecutor->synchronize();
  this->_assemblyEvent.stop();

  // TODO: Add Polynomial == SEPARATE case

  auto iterationCriterion = gko::share(gko::stop::Iteration::build()
                                           .with_max_iters(static_cast<std::size_t>(1e6))
                                           .on(this->_deviceExecutor));

  auto residualCriterion = gko::share(gko::stop::ResidualNormReduction<>::build()
                                          .with_reduction_factor(ginkgoParameter.residualNorm)
                                          .on(this->_deviceExecutor));

  iterationCriterion->add_logger(this->_logger);
  residualCriterion->add_logger(this->_logger);

  if (this->_solverType == SolverType::MG) {

    // TODO: Add loggers to each step here

    auto smootherFactory = gko::share(
        ir::build()
            .with_solver(jacobi::build().with_max_block_size(1u).on(this->_deviceExecutor))
            .with_relaxation_factor(0.9)
            .with_criteria(
                gko::stop::Iteration::build().with_max_iters(2u).on(this->_deviceExecutor))
            .on(this->_deviceExecutor));

    auto mgLevelFactory = amgx_pgm::build().with_deterministic(false).on(this->_deviceExecutor);

    auto coarsestFactory =
        cg::build()
            .with_criteria(
                gko::stop::Iteration::build().with_max_iters(4u).on(this->_deviceExecutor))
            .on(this->_deviceExecutor);

    auto multigridFactory =
        mg::build()
            .with_max_levels(2u)
            .with_min_coarse_rows(2u) // TODO: Check how to configure best
            .with_pre_smoother(gko::share(smootherFactory))
            .with_post_uses_pre(true)
            .with_mg_level(gko::share(mgLevelFactory))
            .with_coarsest_solver(
                gko::share(jacobi::build().with_max_block_size(1u).on(this->_deviceExecutor)))
            .with_criteria(residualCriterion)
            .on(this->_deviceExecutor);

    this->_mgSolver = gko::share(multigridFactory->generate(this->_rbfSystemMatrix));
    this->_mgSolver->add_logger(this->_logger);

  } else if (this->_solverType == SolverType::CG) {

    if (PreconditionerType::None != this->_preconditionerType) {
      auto solverFactoryWithPreconditioner = [preconditionerType = this->_preconditionerType, executor = this->_deviceExecutor, &ginkgoParameter]() {
        if (preconditionerType == PreconditionerType::Jacobi) {
          return cg::build().with_preconditioner(jacobi::build().with_max_block_size(ginkgoParameter.jacobiBlockSize).on(executor));
        } else if (preconditionerType == PreconditionerType::Cholesky) {
          return cg::build().with_preconditioner(cholesky::build().on(executor));
        } else {
          return cg::build().with_preconditioner(ilu::build().on(executor));
        }
      }();

      auto solverFactory = solverFactoryWithPreconditioner
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(this->_deviceExecutor);

      this->_cgSolver = gko::share(solverFactory->generate(this->_rbfSystemMatrix));
      this->_cgSolver->add_logger(this->_logger);
    }

    else {
      auto solverFactory = cg::build()
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(this->_deviceExecutor);

      this->_cgSolver = gko::share(solverFactory->generate(this->_rbfSystemMatrix));
      this->_cgSolver->add_logger(this->_logger);
    }

  } else if (this->_solverType == SolverType::GMRES) {

    if (PreconditionerType::None != this->_preconditionerType) {
      auto solverFactoryWithPreconditioner = [preconditionerType = this->_preconditionerType, executor = _deviceExecutor, &ginkgoParameter]() {
        if (preconditionerType == PreconditionerType::Jacobi) {
          return gmres::build().with_preconditioner(jacobi::build().with_max_block_size(ginkgoParameter.jacobiBlockSize).on(executor));
        } else if (preconditionerType == PreconditionerType::Cholesky) {
          return gmres::build().with_preconditioner(cholesky::build().on(executor));
        } else {
          return gmres::build().with_preconditioner(ilu::build().on(executor));
        }
      }();

      auto solverFactory = solverFactoryWithPreconditioner
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(this->_deviceExecutor);

      this->_gmresSolver = gko::share(solverFactory->generate(this->_rbfSystemMatrix));
      this->_gmresSolver->add_logger(this->_logger);
    } else {
      auto solverFactory = gmres::build()
                               .with_criteria(iterationCriterion, residualCriterion)
                               .on(this->_deviceExecutor);

      this->_gmresSolver = gko::share(solverFactory->generate(this->_rbfSystemMatrix));
      this->_gmresSolver->add_logger(this->_logger);
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::_solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const
{
  precice::utils::Event e("map.rbf.ginkgo.solveSystemMatrix");
  if (this->_solverType == SolverType::CG) {
    this->_cgSolver->apply(gko::lend(rhs), gko::lend(this->_rbfCoefficients));
  } else if (this->_solverType == SolverType::GMRES) {
    this->_gmresSolver->apply(gko::lend(rhs), gko::lend(this->_rbfCoefficients));
  } else if (this->_solverType == SolverType::MG) {
    this->_mgSolver->apply(gko::lend(rhs), gko::lend(this->_rbfCoefficients));
  }

// Only compute time-consuming statistics in debug mode
#ifndef NDEBUG

  auto residual = gko::initialize<GinkgoScalar>({0.0}, this->_deviceExecutor);
  this->_rbfSystemMatrix->apply(gko::lend(this->_scalarOne), gko::lend(this->_rbfCoefficients), gko::lend(this->_scalarNegativeOne), gko::lend(rhs));
  rhs->compute_norm2(gko::lend(residual));
  residual = gko::clone(this->_hostExecutor, residual);
  PRECICE_INFO("Ginkgo Solver Iteration Count: {}", this->_logger->get_num_iterations());
  PRECICE_INFO("Ginkgo Solver Final Residual: {}", residual->at(0, 0));

#endif
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(const Eigen::VectorXd &rhsValues, Polynomial polynomial)
{
  PRECICE_ASSERT(rhsValues.cols() == 1);
  // Copy rhs vector onto GPU by creating a Ginkgo Vector
  auto rhs = gko::share(GinkgoVector::create(this->_hostExecutor, gko::dim<2>{static_cast<unsigned long>(rhsValues.rows()), 1}));

  for (Eigen::Index i = 0; i < rhsValues.rows(); ++i) {
    rhs->at(i, 0) = rhsValues(i, 0);
  }

  this->_copyEvent.start();
  rhs = gko::clone(this->_deviceExecutor, rhs);
  this->_copyEvent.pause();

  if (polynomial == Polynomial::SEPARATE) {
    // TODO: Check if there is least squares solver
    auto polynomialSolverFactory = cg::build()
                                       .with_criteria(gko::stop::Iteration::build()
                                                          .with_max_iters(static_cast<std::size_t>(1e6))
                                                          .on(this->_deviceExecutor),
                                                      gko::stop::ResidualNormReduction<>::build()
                                                          .with_reduction_factor(1e-4)
                                                          .on(this->_deviceExecutor))
                                       .on(this->_deviceExecutor);

    auto matrixQ_T     = gko::share(this->_matrixQ->transpose());
    auto matrixQTQ     = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{matrixQ_T->get_size()[0], this->_matrixQ->get_size()[1]}));
    auto polynomialRHS = gko::share(GinkgoVector::create(this->_deviceExecutor, gko::dim<2>{matrixQ_T->get_size()[0], 1}));
    matrixQ_T->apply(gko::lend(this->_matrixQ), gko::lend(matrixQTQ));
    matrixQ_T->apply(gko::lend(rhs), gko::lend(polynomialRHS));

    auto polynomialSolver         = polynomialSolverFactory->generate(matrixQTQ);
    this->_polynomialContribution = gko::share(GinkgoVector::create(this->_deviceExecutor, gko::dim<2>{matrixQTQ->get_size()[1], 1}));
    polynomialSolver->apply(gko::lend(polynomialRHS), gko::lend(this->_polynomialContribution));

    auto subPolynomialContribution = gko::share(GinkgoVector::create(this->_deviceExecutor, gko::dim<2>{this->_matrixQ->get_size()[0], 1}));
    this->_matrixQ->apply(gko::lend(this->_polynomialContribution), gko::lend(subPolynomialContribution));
    rhs->sub_scaled(gko::lend(this->_scalarOne), gko::lend(subPolynomialContribution));
  }

  this->_solveRBFSystem(rhs);

  auto output = gko::share(GinkgoVector::create(this->_deviceExecutor, gko::dim<2>{this->_matrixA->get_size()[0], this->_rbfCoefficients->get_size()[1]}));

  this->_matrixA->apply(gko::lend(this->_rbfCoefficients), gko::lend(output));

  if (polynomial == Polynomial::SEPARATE) {
    auto addPolynomialContribution = gko::share(GinkgoVector::create(this->_deviceExecutor, gko::dim<2>{this->_matrixV->get_size()[0], 1}));
    this->_matrixV->apply(gko::lend(this->_polynomialContribution), gko::lend(addPolynomialContribution));
    output->add_scaled(gko::lend(this->_scalarOne), gko::lend(addPolynomialContribution));
  }

  this->_copyEvent.start();
  output = gko::clone(this->_hostExecutor, output);

  // TODO: Check if rather use Ginkgo throughout process instead of Eigen
  Eigen::VectorXd result(output->get_size()[0], 1);

  // TODO: Check if rather put this procedure into function
  for (Eigen::Index i = 0; i < result.rows(); ++i) {
    result(i, 0) = output->at(i, 0);
  }
  this->_copyEvent.pause();

  return result;
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::solveConservative(const Eigen::VectorXd &rhsValues, Polynomial polynomial)
{
  return Eigen::VectorXd(1, 1);
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::shared_ptr<gko::Executor> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::_hostExecutor = gko::ReferenceExecutor::create();

template <typename RADIAL_BASIS_FUNCTION_T>
std::shared_ptr<gko::Executor> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::_ompExecutor = gko::OmpExecutor::create();

template <typename RADIAL_BASIS_FUNCTION_T>
std::shared_ptr<gko::Executor> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getReferenceExecutor() const
{
  return this->_hostExecutor;
}

template <typename RADIAL_BASIS_FUNCTION_T>
const std::shared_ptr<GinkgoMatrix> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::getEvaluationMatrix() const
{
  return this->_matrixA;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::clear()
{
  if (nullptr == this->_rbfSystemMatrix) {
    this->_rbfSystemMatrix = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{0, 0}));
    this->_matrixA         = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{0, 0}));
  } else {
    this->_rbfSystemMatrix = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{this->_rbfSystemMatrix->get_size()[0], this->_rbfSystemMatrix->get_size()[1]}));
    this->_matrixA         = gko::share(GinkgoMatrix::create(this->_deviceExecutor, gko::dim<2>{this->_matrixA->get_size()[0], this->_matrixA->get_size()[1]}));
  }
}

} // namespace mapping
} // namespace precice

#endif // PRECICE_NO_GINKGO
