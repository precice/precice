#pragma once

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

namespace precice {
namespace mapping {

/**
 * @brief Mapping with radial basis functions.
 *
 * With help of the input data points and values an interpolant is constructed.
 * The interpolant is formed by a weighted sum of conditionally positive radial
 * basis functions and a (low order) polynomial, and evaluated at the output
 * data points.
 *
 * The radial basis function type has to be given as template parameter, and has
 * to be one of the defined types in this file.
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
  void computeMapping() final;

  /// Removes a computed mapping.
  void clear() final;

  /// name of the rbf mapping
  std::string getName() const final;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctMapping"};

  // The actual solver
  std::unique_ptr<SOLVER_T> _rbfSolver;

  /// @copydoc RadialBasisFctBaseMapping::mapConservative
  void mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) final;

  /// @copydoc RadialBasisFctBaseMapping::mapConsistent
  void mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) final;

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

    // Input mesh may have overlaps
    mesh::Mesh filteredInMesh("filteredInMesh", inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
    mesh::filterMesh(filteredInMesh, *inMesh, [&](const mesh::Vertex &v) { return v.isOwner(); });

    // Send the mesh
    com::sendMesh(*utils::IntraComm::getCommunication(), 0, filteredInMesh);
    com::sendMesh(*utils::IntraComm::getCommunication(), 0, *outMesh);

  } else { // Parallel Primary rank or Serial

    mesh::Mesh globalInMesh(inMesh->getName(), inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
    mesh::Mesh globalOutMesh(outMesh->getName(), outMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);

    if (utils::IntraComm::isPrimary()) {
      {
        // Input mesh may have overlaps
        mesh::Mesh filteredInMesh("filteredInMesh", inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
        mesh::filterMesh(filteredInMesh, *inMesh, [&](const mesh::Vertex &v) { return v.isOwner(); });
        globalInMesh.addMesh(filteredInMesh);
        globalOutMesh.addMesh(*outMesh);
      }

      // Receive mesh
      for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
        mesh::Mesh secondaryInMesh(inMesh->getName(), inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
        com::receiveMesh(*utils::IntraComm::getCommunication(), secondaryRank, secondaryInMesh);
        globalInMesh.addMesh(secondaryInMesh);

        mesh::Mesh secondaryOutMesh(outMesh->getName(), outMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
        com::receiveMesh(*utils::IntraComm::getCommunication(), secondaryRank, secondaryOutMesh);
        globalOutMesh.addMesh(secondaryOutMesh);
      }

    } else { // Serial
      globalInMesh.addMesh(*inMesh);
      globalOutMesh.addMesh(*outMesh);
    }

    // Forwarding the tuples here requires some template magic I don't want to implement
    if constexpr (std::tuple_size_v<std::tuple<Args...>>> 0) {
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

template <typename SOLVER_T, typename... Args>
std::string RadialBasisFctMapping<SOLVER_T, Args...>::getName() const
{
  if constexpr (std::tuple_size_v<std::tuple<Args...>>> 0) {
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
void RadialBasisFctMapping<SOLVER_T, Args...>::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.rbf.mapData.From" + this->input()->getName() + "To" + this->output()->getName(), profiling::Synchronize);

  PRECICE_DEBUG("Map conservative using {}", getName());

  // Gather input data
  if (utils::IntraComm::isSecondary()) {

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

    std::vector<double> globalInValues;
    std::vector<double> outputValueSizes;
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

      outputValueSizes.push_back(localOutputSize);
    }

    {
      int secondaryOutputValueSize;
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        std::vector<double> secondaryBuffer = utils::IntraComm::getCommunication()->receiveRange(rank, com::asVector<double>);
        globalInValues.insert(globalInValues.end(), secondaryBuffer.begin(), secondaryBuffer.end());

        utils::IntraComm::getCommunication()->receive(secondaryOutputValueSize, rank);
        outputValueSizes.push_back(secondaryOutputValueSize);
      }
    }

    const int valueDim = inData.dataDims;

    // Construct Eigen vectors
    Eigen::Map<Eigen::VectorXd> inputValues(globalInValues.data(), globalInValues.size());
    Eigen::VectorXd             outputValues((this->output()->getGlobalNumberOfVertices()) * valueDim);

    Eigen::VectorXd in;                     // rows == outputSize
    in.resize(_rbfSolver->getOutputSize()); // rows == outputSize

    outputValues.setZero();

    for (int dim = 0; dim < valueDim; dim++) {
      for (int i = 0; i < in.size(); i++) { // Fill input data values
        in[i] = inputValues(i * valueDim + dim);
      }

      Eigen::VectorXd out;
      out = _rbfSolver->solveConservative(in, _polynomial);

      // Copy mapped data to output data values
      for (int i = 0; i < this->output()->getGlobalNumberOfVertices(); i++) {
        outputValues[i * valueDim + dim] = out[i];
      }
    }

    // Data scattering to secondary ranks
    if (utils::IntraComm::isPrimary()) {

      // Filter data
      int outputCounter = 0;
      for (int i = 0; i < static_cast<int>(this->output()->nVertices()); ++i) {
        if (this->output()->vertex(i).isOwner()) {
          for (int dim = 0; dim < valueDim; ++dim) {
            outData[i * valueDim + dim] = outputValues(outputCounter);
            ++outputCounter;
          }
        }
      }

      // Data scattering to secondary ranks
      int beginPoint = outputValueSizes.at(0);
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        precice::span<const double> toSend{outputValues.data() + beginPoint, static_cast<size_t>(outputValueSizes.at(rank))};
        utils::IntraComm::getCommunication()->sendRange(toSend, rank);
        beginPoint += outputValueSizes.at(rank);
      }
    } else { // Serial
      outData = outputValues;
    }
  }
  if (utils::IntraComm::isSecondary()) {
    std::vector<double> receivedValues = utils::IntraComm::getCommunication()->receiveRange(0, com::asVector<double>);

    const int valueDim = inData.dataDims;

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

  // Gather input data
  if (utils::IntraComm::isSecondary()) {
    // Input data is filtered
    auto localInDataFiltered = this->input()->getOwnedVertexData(inData.values);
    int  localOutputSize     = outData.size();

    // Send data and output size
    utils::IntraComm::getCommunication()->sendRange(localInDataFiltered, 0);
    utils::IntraComm::getCommunication()->send(localOutputSize, 0);

  } else { // Primary rank or Serial case

    const int valueDim = inData.dataDims;

    std::vector<double> globalInValues(static_cast<std::size_t>(this->input()->getGlobalNumberOfVertices()) * valueDim, 0.0);
    std::vector<int>    outValuesSize;

    if (utils::IntraComm::isPrimary()) { // Parallel case

      // Filter input data
      const auto &localInData = this->input()->getOwnedVertexData(inData.values);
      std::copy(localInData.data(), localInData.data() + localInData.size(), globalInValues.begin());
      outValuesSize.push_back(outData.size());

      int inputSizeCounter = localInData.size();
      int secondaryOutDataSize{0};

      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        std::vector<double> secondaryBuffer = utils::IntraComm::getCommunication()->receiveRange(rank, com::asVector<double>);
        std::copy(secondaryBuffer.begin(), secondaryBuffer.end(), globalInValues.begin() + inputSizeCounter);
        inputSizeCounter += secondaryBuffer.size();

        utils::IntraComm::getCommunication()->receive(secondaryOutDataSize, rank);
        outValuesSize.push_back(secondaryOutDataSize);
      }

    } else { // Serial case
      const auto &localInData = inData.values;
      std::copy(localInData.data(), localInData.data() + localInData.size(), globalInValues.begin());
      outValuesSize.push_back(outData.size());
    }

    Eigen::VectorXd in;

    in.resize(_rbfSolver->getInputSize()); // rows == n
    in.setZero();

    // Construct Eigen vectors
    Eigen::Map<Eigen::VectorXd> inputValues(globalInValues.data(), globalInValues.size());

    Eigen::VectorXd outputValues;
    outputValues.resize((_rbfSolver->getOutputSize()) * valueDim);

    Eigen::VectorXd out;
    outputValues.setZero();

    // For every data dimension, perform mapping
    for (int dim = 0; dim < valueDim; dim++) {
      // Fill input from input data values (last polyparams entries remain zero)
      for (int i = 0; i < this->input()->getGlobalNumberOfVertices(); i++) {
        in[i] = inputValues[i * valueDim + dim];
      }

      out = _rbfSolver->solveConsistent(in, _polynomial);

      // Copy mapped data to output data values
      for (int i = 0; i < out.size(); i++) {
        outputValues[i * valueDim + dim] = out[i];
      }
    }

    outData = Eigen::Map<Eigen::VectorXd>(outputValues.data(), outValuesSize.at(0));

    // Data scattering to secondary ranks
    int beginPoint = outValuesSize.at(0);

    if (utils::IntraComm::isPrimary()) {
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        precice::span<const double> toSend{outputValues.data() + beginPoint, static_cast<size_t>(outValuesSize.at(rank))};
        utils::IntraComm::getCommunication()->sendRange(toSend, rank);
        beginPoint += outValuesSize.at(rank);
      }
    }
  }
  if (utils::IntraComm::isSecondary()) {
    std::vector<double> receivedValues = utils::IntraComm::getCommunication()->receiveRange(0, com::asVector<double>);
    outData                            = Eigen::Map<Eigen::VectorXd>(receivedValues.data(), receivedValues.size());
  }
}
} // namespace mapping
} // namespace precice
