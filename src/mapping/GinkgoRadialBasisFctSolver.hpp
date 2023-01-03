#pragma once
#ifndef PRECICE_NO_GINKGO

#include <array>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/irange.hpp>
#include <functional>
#include <ginkgo/ginkgo.hpp>
#include <ginkgo/kernels/kernel_declaration.hpp>
#include <numeric>
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/impl/DeviceBasisFunctions.cuh"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"

// Declare Ginkgo Kernels
GKO_DECLARE_UNIFIED(template <typename ValueType, typename EvalFunctionType> void create_rbf_system_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *supportPoints,
    ValueType *targetPoints, EvalFunctionType f, const std::array<ValueType, 3> rbf_params,
    const bool addPolynomial, const unsigned int extraDims = 0));

GKO_DECLARE_UNIFIED(template <typename ValueType> void fill_polynomial_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *x, const unsigned int dims = 4));

GKO_DECLARE_UNIFIED(template <typename ValueType> void extract_upper_triangular(
    std::shared_ptr<const DefaultExecutor> exec,
    ValueType *src, ValueType *dest,
    const std::size_t i, const std::size_t j, const std::size_t N));

GKO_REGISTER_UNIFIED_OPERATION(rbf_fill_operation, create_rbf_system_matrix);
GKO_REGISTER_UNIFIED_OPERATION(polynomial_fill_operation, fill_polynomial_matrix);
GKO_REGISTER_UNIFIED_OPERATION(tril_operation, extract_upper_triangular);

namespace precice {

extern bool syncMode;
namespace mapping {

// Ginkgo Data Structures
using GinkgoVector = gko::matrix::Dense<double>;
using GinkgoMatrix = gko::matrix::Dense<double>;
// Ginkgo Solver
using cg    = gko::solver::Cg<>;
using gmres = gko::solver::Gmres<>;
using mg    = gko::solver::Multigrid;
// Ginkgo Preconditioner
using jacobi   = gko::preconditioner::Jacobi<>;
using cholesky = gko::preconditioner::Ic<>;
using ilu      = gko::preconditioner::Ilu<>;

enum SolverType {
  CG,
  GMRES,
  MG
};

enum PreconditionerType {
  Jacobi,
  Cholesky,
  Ilu
};

const std::map<std::string, SolverType> solverTypeLookup{
    {"ginkgo-cg-solver", SolverType::CG},
    {"ginkgo-gmres-solver", SolverType::GMRES},
    {"ginkgo-mg-solver", SolverType::MG}};

const std::map<std::string, PreconditionerType> preconditionerTypeLookup{
    {"ginkgo-jacobi-preconditioner", PreconditionerType::Jacobi},
    {"ginkgo-cholesky-preconditioner", PreconditionerType::Cholesky},
    {"ginkgo-ilu-preconditioner", PreconditionerType::Ilu}};

const std::map<std::string, std::function<std::shared_ptr<gko::Executor>()>> ginkgoExecutorLookup{{"ginkgo-reference-executor", [] { return gko::ReferenceExecutor::create(); }},
                                                                                                  {"ginkgo-omp-executor", [] { return gko::OmpExecutor::create(); }},
                                                                                                  {"ginkgo-cuda-executor", [] { return gko::CudaExecutor::create(0, gko::OmpExecutor::create(), true); }},
                                                                                                  {"ginkgo-hip-executor", [] { return gko::HipExecutor::create(0, gko::OmpExecutor::create(), true); }}};

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
  precice::logging::Logger _log{"mapping::GinkgoRadialBasisFctSolver"};

  std::shared_ptr<gko::Executor>        _deviceExecutor;
  static std::shared_ptr<gko::Executor> _hostExecutor;

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
  std::shared_ptr<GinkgoMatrix> _scalarOne;

  void _solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const;

  precice::utils::Event _copyEvent{"map.rbf.ginkgo.memCopy", precice::syncMode, false};
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

  this->_scalarOne           = gko::share(GinkgoMatrix::create(this->_hostExecutor, gko::dim<2>{1, 1}));
  this->_scalarOne->at(0, 0) = 1.0;
  // Copying it to device
  this->_scalarOne = gko::clone(this->_deviceExecutor, this->_scalarOne);

  // Now we fill the RBF system matrix on the GPU (or any other selected device)
  this->_rbfCoefficients = gko::share(GinkgoVector::create(this->_hostExecutor, gko::dim<2>{n, 1}));

  // We need to copy the input data into a CPU stored vector first and copy it to the GPU afterwards
  auto inputVertices = gko::share(GinkgoMatrix::create(this->_hostExecutor, gko::dim<2>{inputMeshSize, meshDim}));
  for (std::size_t i = 0; i < inputMeshSize; ++i) {
    for (std::size_t j = 0; j < meshDim; ++j) {
      inputVertices->at(i, j) = inputMesh.vertices().at(i).rawCoords()[j];
    }
    // Initial guess is required since memory chunk could lead to never converging system
    if (i < n) {
      this->_rbfCoefficients->at(i, 0) = 0.0;
    }
  }

  auto outputVertices = gko::share(GinkgoMatrix::create(this->_hostExecutor, gko::dim<2>{outputMeshSize, meshDim}));
  for (std::size_t i = 0; i < outputMeshSize; ++i) {
    for (std::size_t j = 0; j < meshDim; ++j) {
      outputVertices->at(i, j) = outputMesh.vertices().at(i).rawCoords()[j];
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

    this->_deviceExecutor->run(make_polynomial_fill_operation(this->_matrixQ->get_size()[0], this->_matrixQ->get_size()[1], this->_matrixQ->get_values(), inputVertices->get_values(), polyParams));
    this->_deviceExecutor->run(make_polynomial_fill_operation(this->_matrixV->get_size()[0], this->_matrixV->get_size()[1], this->_matrixV->get_values(), outputVertices->get_values(), polyParams));

    this->_deviceExecutor->synchronize();
  }

  // Launch RBF fill kernel on device
  precice::utils::Event fillEvent("map.rbf.ginkgo.fillMatrices", precice::syncMode);
  this->_deviceExecutor->run(make_rbf_fill_operation(this->_rbfSystemMatrix->get_size()[0], this->_rbfSystemMatrix->get_size()[1], this->_rbfSystemMatrix->get_values(), inputVertices->get_values(), inputVertices->get_values(), basisFunction.getFunctor(), basisFunction.getFunctionParameters(), polynomial == Polynomial::ON, polyparams));
  this->_deviceExecutor->run(make_rbf_fill_operation(this->_matrixA->get_size()[0], this->_matrixA->get_size()[1], this->_matrixA->get_values(), inputVertices->get_values(), outputVertices->get_values(), basisFunction.getFunctor(), basisFunction.getFunctionParameters(), polynomial == Polynomial::ON, polyparams));

  // Wait for the kernels to finish
  this->_deviceExecutor->synchronize();
  fillEvent.stop();

  // TODO: Add Polynomial == SEPARATE case

  if (this->_solverType == SolverType::MG) {

  } else if (this->_solverType == SolverType::CG) {

    auto solverFactoryWithPreconditioner = [preconditionerType = this->_preconditionerType, executor = this->_deviceExecutor]() {
      if (preconditionerType == PreconditionerType::Jacobi) {
        return cg::build().with_preconditioner(jacobi::build().on(executor));
      } else if (preconditionerType == PreconditionerType::Cholesky) {
        return cg::build().with_preconditioner(cholesky::build().on(executor));
      } else {
        return cg::build().with_preconditioner(ilu::build().on(executor));
      }
    }();

    auto solverFactory = solverFactoryWithPreconditioner
                             .with_criteria(gko::stop::Iteration::build()
                                                .with_max_iters(static_cast<std::size_t>(1e6))
                                                .on(this->_deviceExecutor),
                                            gko::stop::ResidualNormReduction<>::build()
                                                .with_reduction_factor(1e-4)
                                                .on(this->_deviceExecutor))
                             .on(this->_deviceExecutor);

    this->_cgSolver = gko::share(solverFactory->generate(this->_rbfSystemMatrix));

  } else if (this->_solverType == SolverType::GMRES) {

    auto solverFactoryWithPreconditioner = [preconditionerType = this->_preconditionerType, executor = _deviceExecutor]() {
      if (preconditionerType == PreconditionerType::Jacobi) {
        return gmres::build().with_preconditioner(jacobi::build().on(executor));
      } else if (preconditionerType == PreconditionerType::Cholesky) {
        return gmres::build().with_preconditioner(cholesky::build().on(executor));
      } else {
        return gmres::build().with_preconditioner(ilu::build().on(executor));
      }
    }();

    auto solverFactory = solverFactoryWithPreconditioner
                             .with_criteria(gko::stop::Iteration::build()
                                                .with_max_iters(static_cast<std::size_t>(1e6))
                                                .on(this->_deviceExecutor),
                                            gko::stop::ResidualNormReduction<>::build()
                                                .with_reduction_factor(1e-4)
                                                .on(this->_deviceExecutor))
                             .on(this->_deviceExecutor);

    this->_gmresSolver = gko::share(solverFactory->generate(this->_rbfSystemMatrix));
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::_solveRBFSystem(const std::shared_ptr<GinkgoVector> &rhs) const
{
  precice::utils::Event e("map.rbf.ginkgo.solveSystemMatrix", precice::syncMode);
  if (this->_solverType == SolverType::CG) {
    this->_cgSolver->apply(gko::lend(rhs), gko::lend(this->_rbfCoefficients));
  } else if (this->_solverType == SolverType::GMRES) {
    this->_gmresSolver->apply(gko::lend(rhs), gko::lend(this->_rbfCoefficients));
  } else if (this->_solverType == SolverType::MG) {
    // TODO:
  }
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
std::shared_ptr<gko::Executor> GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>::_hostExecutor = ginkgoExecutorLookup.at("ginkgo-reference-executor")();

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
