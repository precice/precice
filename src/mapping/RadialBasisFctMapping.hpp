#pragma once

/**
 * @file RadialBasisFctMapping.hpp
 * @brief Concrete RBF mapping class that combines a basis function kernel with a linear solver.
 *
 * ## Design Overview
 * This class sits above `RadialBasisFctBaseMapping` and introduces the actual
 * **solver** as a second template parameter (`SOLVER_T`). The solver encapsulates:
 *   - How the interpolation matrix is assembled (CPU dense, GPU batched, PETSc parallel)
 *   - The matrix decomposition used (Cholesky for SPD kernels, QR otherwise)
 *   - How the linear system is solved at mapping time
 *
 * ## Parallel Strategy (Serial / MPI)
 * On a single process (serial), the mapping is straightforward:
 *   1. Build matrix C from all input vertices
 *   2. Decompose C
 *   3. Build evaluation matrix A from all output vertices
 *   4. At map time: solve C*lambda = data, then output = A*lambda
 *
 * On multiple ranks (MPI parallel), preCICE uses a **gather-on-primary** strategy:
 *   - Each **secondary rank** filters its owned vertices and sends them + data to rank 0.
 *   - **Primary rank** assembles the global mesh and builds the solver.
 *   - At mapping time, primary solves globally, then **scatters** results back.
 *
 * This keeps the RBF system global (no partitioned sub-problems), at the cost of
 * communication. For partition-of-unity see `PartitionOfUnityMapping`.
 *
 * ## Variadic Template Args
 * Optional `Args...` are forwarded to the `SOLVER_T` constructor. This is used to
 * pass GPU executor configuration (Ginkgo parameter structs) when a GPU solver is
 * selected without changing the interface for CPU-only builds.
 */

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "com/Communication.hpp"
#include "com/Extra.hpp"
#include "config/MappingConfiguration.hpp"
#include "mapping/RadialBasisFctBaseMapping.hpp"
#include "mesh/Filter.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "utils/IntraComm.hpp"

namespace precice::mapping {

/**
 * @brief Concrete RBF mapping: selects a solver and a basis function at compile time.
 *
 * @tparam SOLVER_T  The linear-system solver type. Must expose:
 *   - `BASIS_FUNCTION_T` alias (the RBF kernel type)
 *   - Constructor `(basisFn, inMesh, inIDs, outMesh, outIDs, deadAxis, polynomial [, extraArg])`
 *   - `solveConsistent(inputData, polynomial)` → mapped output matrix
 *   - `solveConservative(inputData, polynomial)` → mapped output matrix
 *   - `getInputSize()`, `getOutputSize()`
 *
 * @tparam Args  Optional extra constructor arguments for SOLVER_T
 *   (e.g., a `GinkgoParameter` struct for GPU executor selection).
 */
template <typename SOLVER_T, typename... Args>
class RadialBasisFctMapping : public RadialBasisFctBaseMapping<typename SOLVER_T::BASIS_FUNCTION_T> {
public:
  using RADIAL_BASIS_FUNCTION_T = typename SOLVER_T::BASIS_FUNCTION_T;
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   * @param[in] function Radial basis function used for mapping.
   * @param[in] xDead, yDead, zDead Deactivates mapping along an axis
   */
  RadialBasisFctMapping(
      Mapping::Constraint     constraint,
      int                     dimensions,
      RADIAL_BASIS_FUNCTION_T function,
      std::array<bool, 3>     deadAxis,
      Polynomial              polynomial,
      Args... args);

  /// Computes the mapping coefficients from the in- and output mesh.
  void computeMapping() final override;

  /// Removes a computed mapping.
  void clear() final override;

  /// name of the rbf mapping
  std::string getName() const final override;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctMapping"};

  // The actual solver
  std::unique_ptr<SOLVER_T> _rbfSolver;

  /// @copydoc RadialBasisFctBaseMapping::mapConservative
  void mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) final override;

  /// @copydoc RadialBasisFctBaseMapping::mapConsistent
  void mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) final override;

  /// Treatment of the polynomial
  Polynomial _polynomial;

  /// Optional constructor arguments for the solver class
  std::tuple<Args...> optionalArgs;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename SOLVER_T, typename... Args>
RadialBasisFctMapping<SOLVER_T, Args...>::RadialBasisFctMapping(
    Mapping::Constraint     constraint,
    int                     dimensions,
    RADIAL_BASIS_FUNCTION_T function,
    std::array<bool, 3>     deadAxis,
    Polynomial              polynomial,
    Args... args)
    : RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>(constraint, dimensions, function, deadAxis, Mapping::InitialGuessRequirement::None),
      _polynomial(polynomial),
      optionalArgs(std::make_tuple(std::forward<Args>(args)...))
{
  PRECICE_CHECK(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && polynomial == Polynomial::ON), "The integrated polynomial (polynomial=\"on\") is not supported for the selected radial-basis function. Please select another radial-basis function or change the polynomial configuration.");
}

template <typename SOLVER_T, typename... Args>
void RadialBasisFctMapping<SOLVER_T, Args...>::computeMapping()
{
  PRECICE_TRACE();

  precice::profiling::Event e("map.rbf.computeMapping.From" + this->input()->getName() + "To" + this->output()->getName(), profiling::Synchronize);

  PRECICE_ASSERT(this->input()->getDimensions() == this->output()->getDimensions(),
                 this->input()->getDimensions(), this->output()->getDimensions());
  PRECICE_ASSERT(this->getDimensions() == this->output()->getDimensions(),
                 this->getDimensions(), this->output()->getDimensions());

  // ROLE ASSIGNMENT FOR CONSISTENT vs CONSERVATIVE
  //
  // The RBF solver always treats its first mesh as the "source" (kernel evaluated here)
  // and second mesh as the "target" (evaluation points). For conservative mappings the
  // mathematical roles of input/output are reversed relative to the physical roles:
  //   consistent:   inMesh = source (known values), outMesh = target (mapped values)
  //   conservative: inMesh = output() = target of data flow, outMesh = input() = source
  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;

  if (this->hasConstraint(Mapping::CONSERVATIVE)) {
    inMesh  = this->output();
    outMesh = this->input();
  } else { // Consistent or scaled consistent
    inMesh  = this->input();
    outMesh = this->output();
  }

  if (utils::IntraComm::isSecondary()) {
    // SECONDARY RANK: Filter local owned vertices and send to primary for global assembly.
    //
    // The input mesh can have halo/ghost overlaps; only owned vertices contribute
    // to the global problem to avoid double-counting.
    mesh::Mesh filteredInMesh("filteredInMesh", inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
    mesh::filterMesh(filteredInMesh, *inMesh, [&](const mesh::Vertex &v) { return v.isOwner(); });

    // Send both the filtered input mesh and the local output mesh to primary (rank 0).
    com::sendMesh(*utils::IntraComm::getCommunication(), 0, filteredInMesh);
    com::sendMesh(*utils::IntraComm::getCommunication(), 0, *outMesh);

  } else { // Parallel Primary rank or Serial

    // PRIMARY RANK (or serial): Assemble the global mesh from local + all secondary contributions.
    mesh::Mesh globalInMesh(inMesh->getName(), inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
    mesh::Mesh globalOutMesh(outMesh->getName(), outMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);

    if (utils::IntraComm::isPrimary()) {
      {
        // Add primary's own owned vertices first.
        mesh::Mesh filteredInMesh("filteredInMesh", inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
        mesh::filterMesh(filteredInMesh, *inMesh, [&](const mesh::Vertex &v) { return v.isOwner(); });
        globalInMesh.addMesh(filteredInMesh);
        globalOutMesh.addMesh(*outMesh);
      }

      // Receive the filtered mesh pieces from each secondary rank and merge them.
      for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
        mesh::Mesh secondaryInMesh(inMesh->getName(), inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
        com::receiveMesh(*utils::IntraComm::getCommunication(), secondaryRank, secondaryInMesh);
        globalInMesh.addMesh(secondaryInMesh);

        mesh::Mesh secondaryOutMesh(outMesh->getName(), outMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
        com::receiveMesh(*utils::IntraComm::getCommunication(), secondaryRank, secondaryOutMesh);
        globalOutMesh.addMesh(secondaryOutMesh);
      }

    } else { // Serial (single-process) case: global == local mesh.
      globalInMesh.addMesh(*inMesh);
      globalOutMesh.addMesh(*outMesh);
    }

    // Build the RBF solver using the globally assembled meshes.
    // The solver internally assembles matrix C (input×input), decomposes it,
    // and builds the evaluation matrix A (output×input) for later mapping calls.
    //
    // The `if constexpr` branch handles optional extra arguments (e.g., GPU config)
    // passed through variadic template Args.
    if constexpr (sizeof...(Args) > 0) {
      _rbfSolver = std::make_unique<SOLVER_T>(this->_basisFunction, globalInMesh, boost::irange<Eigen::Index>(0, globalInMesh.nVertices()),
                                              globalOutMesh, boost::irange<Eigen::Index>(0, globalOutMesh.nVertices()), this->_deadAxis, _polynomial, std::get<0>(optionalArgs));
    } else {
      _rbfSolver = std::make_unique<SOLVER_T>(this->_basisFunction, globalInMesh, boost::irange<Eigen::Index>(0, globalInMesh.nVertices()),
                                              globalOutMesh, boost::irange<Eigen::Index>(0, globalOutMesh.nVertices()), this->_deadAxis, _polynomial);
    }
  }
  this->_hasComputedMapping = true;
  PRECICE_DEBUG("Compute Mapping is Completed.");
}

template <typename SOLVER_T, typename... Args>
void RadialBasisFctMapping<SOLVER_T, Args...>::clear()
{
  PRECICE_TRACE();
  _rbfSolver.reset();
  this->_hasComputedMapping = false;
}

template <typename... Args>
static std::string getNameWithArgs(const std::tuple<Args...> &optionalArgs)
{
  if constexpr (sizeof...(Args) > 0) {
    auto        param = std::get<0>(optionalArgs);
    std::string exec  = param.executor;
    if (param.solver == "qr-solver") {
      return "global-direct RBF (" + exec + ")";
    } else {
      return "global-iterative RBF (" + exec + ")";
    }
  } else {
    return "global-direct RBF (cpu-executor)";
  }
}

template <typename SOLVER_T, typename... Args>
std::string RadialBasisFctMapping<SOLVER_T, Args...>::getName() const
{
  return getNameWithArgs(optionalArgs);
}

template <typename SOLVER_T, typename... Args>
void RadialBasisFctMapping<SOLVER_T, Args...>::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.rbf.mapData.From" + this->input()->getName() + "To" + this->output()->getName(), profiling::Synchronize);

  PRECICE_DEBUG("Map conservative using {}", getName());

  /**
   * CONSERVATIVE MAPPING ALGORITHM
   *
   * In a conservative mapping, the global sum (or integral) of the field is preserved.
   * Mathematically, if s is the source data and t is the target:
   *   sum(t) == sum(s)  (weighted by integration weights)
   *
   * Implementation:
   *   The solver computes t = A^T * (C^-T * s) which is the transpose of the consistent solve.
   *   Geometrically, this distributes "loads" from output vertices onto input vertices.
   *
   * Parallel strategy (gather-on-primary):
   *   1. Each rank sends its chunk of input data to rank 0.
   *   2. Primary assembles global input, calls _rbfSolver->solveConservative(...).
   *   3. Primary filters owned-vertex results back and scatters to each secondary rank.
   *   4. Secondary ranks receive their portion of the output.
   */

  if (utils::IntraComm::isSecondary()) {
    // SECONDARY: Send all local input data and the expected output size to primary.
    // Output size tells primary how many values to scatter back (per-rank).
    const auto &localInData = inData.values;

    int localOutputSize = 0;
    for (const auto &vertex : this->output()->vertices()) {
      if (vertex.isOwner()) {
        ++localOutputSize;
      }
    }

    localOutputSize *= inData.dataDims;

    utils::IntraComm::getCommunication()->sendRange(localInData, 0);
    utils::IntraComm::getCommunication()->send(localOutputSize, 0);

  } else { // Parallel Primary rank or Serial case

    // PHASE 1: Assemble global input data vector by collecting from all ranks.
    std::vector<double> globalInValues;
    std::vector<double> outputValueSizes; // per-rank output sizes, for scatter step
    {
      const auto &localInData = inData.values;
      globalInValues.insert(globalInValues.begin(), localInData.data(), localInData.data() + localInData.size());

      int localOutputSize = 0;
      for (const auto &vertex : this->output()->vertices()) {
        if (vertex.isOwner()) {
          ++localOutputSize;
        }
      }

      localOutputSize *= inData.dataDims;
      outputValueSizes.push_back(localOutputSize); // Primary's own output size
    }

    {
      int secondaryOutputValueSize;
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        // Receive input data from each secondary and append to global buffer.
        std::vector<double> secondaryBuffer = utils::IntraComm::getCommunication()->receiveRange(rank, com::asVector<double>);
        globalInValues.insert(globalInValues.end(), secondaryBuffer.begin(), secondaryBuffer.end());

        // Also receive its output size for scatter step.
        utils::IntraComm::getCommunication()->receive(secondaryOutputValueSize, rank);
        outputValueSizes.push_back(secondaryOutputValueSize);
      }
    }

    const int valueDim = inData.dataDims;

    // Map from flat std::vector to Eigen row-major matrix for the solver.
    // Shape: (globalOutputVertices x valueDim)
    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MapToMatrix = Eigen::Map<RowMatrixXd, Eigen::Unaligned, Eigen::OuterStride<Eigen::Dynamic>>;

    // PHASE 2: Solve the global conservative system.
    // in matrix: rows = global output vertex count, cols = dataDims
    Eigen::MatrixXd in = MapToMatrix(globalInValues.data(), _rbfSolver->getOutputSize(), valueDim, Eigen::OuterStride(valueDim));

    // solveConservative returns the mapped values at all output vertices.
    // We only take the first getGlobalNumberOfVertices() rows (rest are polynomial entries).
    RowMatrixXd out = _rbfSolver->solveConservative(in, _polynomial).block(0, 0, this->output()->getGlobalNumberOfVertices(), valueDim);

    // Flatten the 2D result into a 1D vector for scattering.
    Eigen::Map<Eigen::VectorXd> outputValues(out.data(), (this->output()->getGlobalNumberOfVertices()) * valueDim);

    // PHASE 3: Scatter results back to each rank.
    if (utils::IntraComm::isPrimary()) {

      // Copy primary's own owned-vertex results from the global output vector.
      int outputCounter = 0;
      for (int i = 0; i < static_cast<int>(this->output()->nVertices()); ++i) {
        if (this->output()->vertex(i).isOwner()) {
          for (int dim = 0; dim < valueDim; ++dim) {
            outData[i * valueDim + dim] = outputValues(outputCounter);
            ++outputCounter;
          }
        }
      }

      // Send contiguous chunks from outputValues to each secondary rank.
      int beginPoint = outputValueSizes.at(0); // Offset after primary's own data
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        precice::span<const double> toSend{outputValues.data() + beginPoint, static_cast<size_t>(outputValueSizes.at(rank))};
        utils::IntraComm::getCommunication()->sendRange(toSend, rank);
        beginPoint += outputValueSizes.at(rank);
      }
    } else { // Serial (single process)
      // All output values are for the single rank — assign directly.
      outData = outputValues;
    }
  }

  // SECONDARY: Receive scattered output values from primary and fill owned-vertex slots.
  if (utils::IntraComm::isSecondary()) {
    std::vector<double> receivedValues = utils::IntraComm::getCommunication()->receiveRange(0, com::asVector<double>);

    const int valueDim = inData.dataDims;

    // Copy received values in order into the output array at owned-vertex positions.
    int outputCounter = 0;
    for (int i = 0; i < static_cast<int>(this->output()->nVertices()); ++i) {
      if (this->output()->vertex(i).isOwner()) {
        for (int dim = 0; dim < valueDim; ++dim) {
          outData[i * valueDim + dim] = receivedValues.at(outputCounter);
          ++outputCounter;
        }
      }
    }
  }
}

template <typename SOLVER_T, typename... Args>
void RadialBasisFctMapping<SOLVER_T, Args...>::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.rbf.mapData.From" + this->input()->getName() + "To" + this->output()->getName(), profiling::Synchronize);

  PRECICE_DEBUG("Map {} using {}", (this->hasConstraint(Mapping::CONSISTENT) ? "consistent" : "scaled-consistent"), getName());

  /**
   * CONSISTENT MAPPING ALGORITHM
   *
   * In a consistent mapping, values at output vertices are interpolated from input vertices.
   * Mathematically: output(y_j) = sum_i lambda_i * phi(||y_j - x_i||) + p(y_j)
   *
   * The solve has two phases:
   *   1. Solve: C * lambda = data  (find RBF coefficients from input data)
   *   2. Evaluate: output = A * lambda  (interpolate at output vertex locations)
   *
   * Only OWNED vertices of each rank contribute their data to the global system
   * (ghost/halo vertices would cause double-counting).
   *
   * Parallel strategy (gather-on-primary):
   *   1. Each secondary rank filters owned input data and sends it + output size to primary.
   *   2. Primary assembles the full global input vector, calls solveConsistent(...).
   *   3. Primary scatters the relevant output portions back to each secondary.
   *   4. Secondaries copy received values into their local output array.
   */

  if (utils::IntraComm::isSecondary()) {
    // SECONDARY: Extract only the data for owned vertices (filter out ghost/halo data),
    // then send it to primary along with how much output data to expect back.
    auto localInDataFiltered = this->input()->getOwnedVertexData(inData.values);
    int  localOutputSize     = outData.size();

    utils::IntraComm::getCommunication()->sendRange(localInDataFiltered, 0);
    utils::IntraComm::getCommunication()->send(localOutputSize, 0);

  } else { // Primary rank or Serial case

    const int valueDim = inData.dataDims;

    // Pre-allocate the global input vector (zero-initialized).
    // Polynomial rows in the solver's system are always zero in the input.
    std::vector<double> globalInValues(static_cast<std::size_t>(this->input()->getGlobalNumberOfVertices()) * valueDim, 0.0);
    std::vector<int>    outValuesSize; // per-rank output sizes for scatter

    if (utils::IntraComm::isPrimary()) {
      // PHASE 1a: Copy primary's own owned-vertex data into the global buffer.
      const auto &localInData = this->input()->getOwnedVertexData(inData.values);
      std::copy(localInData.data(), localInData.data() + localInData.size(), globalInValues.begin());
      outValuesSize.push_back(outData.size());

      int inputSizeCounter = localInData.size();
      int secondaryOutDataSize{0};

      // PHASE 1b: Receive each secondary's owned-vertex data and append to global buffer.
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        std::vector<double> secondaryBuffer = utils::IntraComm::getCommunication()->receiveRange(rank, com::asVector<double>);
        std::copy(secondaryBuffer.begin(), secondaryBuffer.end(), globalInValues.begin() + inputSizeCounter);
        inputSizeCounter += secondaryBuffer.size();

        utils::IntraComm::getCommunication()->receive(secondaryOutDataSize, rank);
        outValuesSize.push_back(secondaryOutDataSize);
      }

    } else { // Serial (single process): the local data is already global.
      const auto &localInData = inData.values;
      std::copy(localInData.data(), localInData.data() + localInData.size(), globalInValues.begin());
      outValuesSize.push_back(outData.size());
    }

    // PHASE 2: Build the solver's input matrix.
    // Shape: (inputSize x valueDim) where inputSize includes polynomial placeholder rows.
    // The first globalNumberOfVertices rows hold the actual data.
    // Polynomial rows at the bottom remain zero as required by the augmented system.
    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using MapToMatrix = Eigen::Map<RowMatrixXd, Eigen::Unaligned, Eigen::OuterStride<Eigen::Dynamic>>;

    Eigen::MatrixXd in = Eigen::MatrixXd::Zero(_rbfSolver->getInputSize(), valueDim);
    in.block(0, 0, this->input()->getGlobalNumberOfVertices(), valueDim) =
        MapToMatrix(globalInValues.data(), this->input()->getGlobalNumberOfVertices(), valueDim, Eigen::OuterStride(valueDim));

    // PHASE 2: Solve: coefficients = C^{-1} * in, output = A * coefficients.
    RowMatrixXd out = _rbfSolver->solveConsistent(in, _polynomial);

    // Primary extracts its own output slice.
    outData = Eigen::Map<Eigen::VectorXd>(out.data(), outValuesSize.at(0));

    // PHASE 3: Scatter output slices to secondary ranks.
    int beginPoint = outValuesSize.at(0);
    if (utils::IntraComm::isPrimary()) {
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        precice::span<const double> toSend{out.data() + beginPoint, static_cast<size_t>(outValuesSize.at(rank))};
        utils::IntraComm::getCommunication()->sendRange(toSend, rank);
        beginPoint += outValuesSize.at(rank);
      }
    }
  }

  // SECONDARY: Receive the scattered output and store it in the local output array.
  if (utils::IntraComm::isSecondary()) {
    std::vector<double> receivedValues = utils::IntraComm::getCommunication()->receiveRange(0, com::asVector<double>);
    outData                            = Eigen::Map<Eigen::VectorXd>(receivedValues.data(), receivedValues.size());
  }
}
} // namespace precice::mapping
