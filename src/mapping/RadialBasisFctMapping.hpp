#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "com/Communication.hpp"
#include "com/Extra.hpp"
#include "config/MappingConfiguration.hpp"
#include "mapping/GinkgoRadialBasisFctSolver.hpp"
#include "mapping/RadialBasisFctBaseMapping.hpp"
#include "mapping/RadialBasisFctSolver.hpp"
#include "mesh/Filter.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"
#include "utils/IntraComm.hpp"

namespace precice {
extern bool syncMode;

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
template <typename RADIAL_BASIS_FUNCTION_T>
class RadialBasisFctMapping : public RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T> {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   * @param[in] function Radial basis function used for mapping.
   * @param[in] xDead, yDead, zDead Deactivates mapping along an axis
   */
  RadialBasisFctMapping(
      Mapping::Constraint                          constraint,
      int                                          dimensions,
      RADIAL_BASIS_FUNCTION_T                      function,
      std::array<bool, 3>                          deadAxis,
      Polynomial                                   polynomial,
      const bool                                   useEigen        = true,
      const MappingConfiguration::GinkgoParameter &ginkgoParameter = MappingConfiguration::GinkgoParameter());

  /// Computes the mapping coefficients from the in- and output mesh.
  void computeMapping() final override;

  /// Removes a computed mapping.
  void clear() final override;

  /// name of the rbf mapping
  std::string getName() const final override;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctMapping"};

#ifndef PRECICE_NO_GINKGO
  GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T> _ginkgoRbfSolver;
#endif

  RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T> _rbfSolver;

  /// @copydoc RadialBasisFctBaseMapping::mapConservative
  void mapConservative(DataID inputDataID, DataID outputDataID) final override;

  /// @copydoc RadialBasisFctBaseMapping::mapConsistent
  void mapConsistent(DataID inputDataID, DataID outputDataID) final override;

  /// Treatment of the polynomial
  Polynomial _polynomial;

  /// Determines which backend to use (Eigen or Ginkgo)
  bool _useEigen;

  /// Ginkgo Configuration
  MappingConfiguration::GinkgoParameter _ginkgoParameter;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctMapping(
    Mapping::Constraint                          constraint,
    int                                          dimensions,
    RADIAL_BASIS_FUNCTION_T                      function,
    std::array<bool, 3>                          deadAxis,
    Polynomial                                   polynomial,
    const bool                                   useEigen,
    const MappingConfiguration::GinkgoParameter &ginkgoParameter)
    : RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>(constraint, dimensions, function, deadAxis),
      _polynomial(polynomial),
      _useEigen(useEigen),
      _ginkgoParameter(ginkgoParameter)
{
  PRECICE_CHECK(!(RADIAL_BASIS_FUNCTION_T::isStrictlyPositiveDefinite() && polynomial == Polynomial::ON), "The integrated polynomial (polynomial=\"on\") is not supported for the selected radial-basis function. Please select another radial-basis function or change the polynomial configuration.");
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  PRECICE_TRACE();

  precice::utils::Event e("map.rbf.computeMapping.From" + this->input()->getName() + "To" + this->output()->getName(), precice::syncMode);

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

    if (_useEigen) {
      _rbfSolver = RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>{this->_basisFunction, globalInMesh, boost::irange<Eigen::Index>(0, globalInMesh.vertices().size()),
                                                                 globalOutMesh, boost::irange<Eigen::Index>(0, globalOutMesh.vertices().size()), this->_deadAxis, _polynomial};
    }
#ifndef PRECICE_NO_GINKGO
    else {
      _ginkgoRbfSolver.~GinkgoRadialBasisFctSolver();
      new (&_ginkgoRbfSolver) GinkgoRadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>{this->_basisFunction, globalInMesh, boost::irange<Eigen::Index>(0, globalInMesh.vertices().size()),
                                                                                  globalOutMesh, boost::irange<Eigen::Index>(0, globalOutMesh.vertices().size()), this->_deadAxis, _polynomial, _ginkgoParameter};
    }
#endif
  }
  this->_hasComputedMapping = true;
  PRECICE_DEBUG("Compute Mapping is Completed.");
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  _rbfSolver.clear();
  this->_hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::string RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::getName() const
{
  return "global-direct RBF";
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::mapConservative(DataID inputDataID, DataID outputDataID)
{
  precice::utils::Event e("map.rbf.mapData.From" + this->input()->getName() + "To" + this->output()->getName(), precice::syncMode);
  PRECICE_TRACE(inputDataID, outputDataID);
  using precice::com::AsVectorTag;

  // Gather input data
  if (utils::IntraComm::isSecondary()) {

    const auto &localInData = this->input()->data(inputDataID)->values();

    int localOutputSize = 0;
    for (const auto &vertex : this->output()->vertices()) {
      if (vertex.isOwner()) {
        ++localOutputSize;
      }
    }

    localOutputSize *= this->output()->data(outputDataID)->getDimensions();

    utils::IntraComm::getCommunication()->sendRange(localInData, 0);
    utils::IntraComm::getCommunication()->send(localOutputSize, 0);

  } else { // Parallel Primary rank or Serial case

    std::vector<double> globalInValues;
    std::vector<double> outputValueSizes;
    {
      const auto &localInData = this->input()->data(inputDataID)->values();
      globalInValues.insert(globalInValues.begin(), localInData.data(), localInData.data() + localInData.size());

      int localOutputSize = 0;
      for (const auto &vertex : this->output()->vertices()) {
        if (vertex.isOwner()) {
          ++localOutputSize;
        }
      }

      localOutputSize *= this->output()->data(outputDataID)->getDimensions();

      outputValueSizes.push_back(localOutputSize);
    }

    {
      int secondaryOutputValueSize;
      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        std::vector<double> secondaryBuffer = utils::IntraComm::getCommunication()->receiveRange(rank, AsVectorTag<double>{});
        globalInValues.insert(globalInValues.end(), secondaryBuffer.begin(), secondaryBuffer.end());

        utils::IntraComm::getCommunication()->receive(secondaryOutputValueSize, rank);
        outputValueSizes.push_back(secondaryOutputValueSize);
      }
    }

    int valueDim = this->output()->data(outputDataID)->getDimensions();

    // Construct Eigen vectors
    Eigen::Map<Eigen::VectorXd> inputValues(globalInValues.data(), globalInValues.size());
    Eigen::VectorXd             outputValues((this->output()->getGlobalNumberOfVertices()) * valueDim);

    Eigen::VectorXd in; // rows == outputSize

    if (_useEigen) {
      in.resize(_rbfSolver.getEvaluationMatrix().rows()); // rows == outputSize
    }
#ifndef PRECICE_NO_GINKGO
    else {
      in.resize(_ginkgoRbfSolver.getEvaluationMatrix()->get_size()[0]); // rows == outputSize
    };
#endif

    outputValues.setZero();

    for (int dim = 0; dim < valueDim; dim++) {
      for (int i = 0; i < in.size(); i++) { // Fill input data values
        in[i] = inputValues(i * valueDim + dim);
      }

      Eigen::VectorXd out;

      if (_useEigen) {
        out = _rbfSolver.solveConservative(in, _polynomial);
      }
#ifndef PRECICE_NO_GINKGO
      else {
        out = _ginkgoRbfSolver.solveConservative(in, _polynomial);
      }
#endif

      // Copy mapped data to output data values
      for (int i = 0; i < this->output()->getGlobalNumberOfVertices(); i++) {
        outputValues[i * valueDim + dim] = out[i];
      }
    }

    // Data scattering to secondary ranks
    if (utils::IntraComm::isPrimary()) {

      // Filter data
      int outputCounter = 0;
      for (int i = 0; i < static_cast<int>(this->output()->vertices().size()); ++i) {
        if (this->output()->vertices()[i].isOwner()) {
          for (int dim = 0; dim < valueDim; ++dim) {
            this->output()->data(outputDataID)->values()[i * valueDim + dim] = outputValues(outputCounter);
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
      this->output()->data(outputDataID)->values() = outputValues;
    }
  }
  if (utils::IntraComm::isSecondary()) {
    std::vector<double> receivedValues = utils::IntraComm::getCommunication()->receiveRange(0, AsVectorTag<double>{});

    int valueDim = this->output()->data(outputDataID)->getDimensions();

    int outputCounter = 0;
    for (int i = 0; i < static_cast<int>(this->output()->vertices().size()); ++i) {
      if (this->output()->vertices()[i].isOwner()) {
        for (int dim = 0; dim < valueDim; ++dim) {
          this->output()->data(outputDataID)->values()[i * valueDim + dim] = receivedValues.at(outputCounter);
          ++outputCounter;
        }
      }
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  precice::utils::Event e("map.rbf.mapData.From" + this->input()->getName() + "To" + this->output()->getName(), precice::syncMode);
  PRECICE_TRACE(inputDataID, outputDataID);
  using precice::com::AsVectorTag;

  // Gather input data
  if (utils::IntraComm::isSecondary()) {
    // Input data is filtered
    auto localInDataFiltered = this->input()->getOwnedVertexData(inputDataID);
    int  localOutputSize     = this->output()->data(outputDataID)->values().size();

    // Send data and output size
    utils::IntraComm::getCommunication()->sendRange(localInDataFiltered, 0);
    utils::IntraComm::getCommunication()->send(localOutputSize, 0);

  } else { // Primary rank or Serial case

    int valueDim = this->output()->data(outputDataID)->getDimensions();

    std::vector<double> globalInValues(static_cast<std::size_t>(this->input()->getGlobalNumberOfVertices()) * valueDim, 0.0);
    std::vector<int>    outValuesSize;

    if (utils::IntraComm::isPrimary()) { // Parallel case

      // Filter input data
      const auto &localInData = this->input()->getOwnedVertexData(inputDataID);
      std::copy(localInData.data(), localInData.data() + localInData.size(), globalInValues.begin());
      outValuesSize.push_back(this->output()->data(outputDataID)->values().size());

      int inputSizeCounter = localInData.size();
      int secondaryOutDataSize{0};

      for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
        std::vector<double> secondaryBuffer = utils::IntraComm::getCommunication()->receiveRange(rank, AsVectorTag<double>{});
        std::copy(secondaryBuffer.begin(), secondaryBuffer.end(), globalInValues.begin() + inputSizeCounter);
        inputSizeCounter += secondaryBuffer.size();

        utils::IntraComm::getCommunication()->receive(secondaryOutDataSize, rank);
        outValuesSize.push_back(secondaryOutDataSize);
      }

    } else { // Serial case
      const auto &localInData = this->input()->data(inputDataID)->values();
      std::copy(localInData.data(), localInData.data() + localInData.size(), globalInValues.begin());
      outValuesSize.push_back(this->output()->data(outputDataID)->values().size());
    }

    Eigen::VectorXd in;

    if (_useEigen) {
      in.resize(_rbfSolver.getEvaluationMatrix().cols()); // rows == n
    }
#ifndef PRECICE_NO_GINKGO
    else {
      in.resize(_ginkgoRbfSolver.getEvaluationMatrix()->get_size()[1]); // rows == n
    }
#endif

    in.setZero();

    // Construct Eigen vectors
    Eigen::Map<Eigen::VectorXd> inputValues(globalInValues.data(), globalInValues.size());

    Eigen::VectorXd outputValues;

    if (_useEigen) {
      outputValues.resize((_rbfSolver.getEvaluationMatrix().rows()) * valueDim);
    }
#ifndef PRECICE_NO_GINKGO
    else {
      outputValues.resize((_ginkgoRbfSolver.getEvaluationMatrix()->get_size()[0]) * valueDim);
    }
#endif

    Eigen::VectorXd out;
    outputValues.setZero();

    // For every data dimension, perform mapping
    for (int dim = 0; dim < valueDim; dim++) {
      // Fill input from input data values (last polyparams entries remain zero)
      for (int i = 0; i < this->input()->getGlobalNumberOfVertices(); i++) {
        in[i] = inputValues[i * valueDim + dim];
      }

      if (_useEigen) {
        out = _rbfSolver.solveConsistent(in, _polynomial);
      }
#ifndef PRECICE_NO_GINKGO
      else {
        out = _ginkgoRbfSolver.solveConsistent(in, _polynomial);
      }
#endif

      // Copy mapped data to output data values
      for (int i = 0; i < out.size(); i++) {
        outputValues[i * valueDim + dim] = out[i];
      }
    }

    this->output()->data(outputDataID)->values() = Eigen::Map<Eigen::VectorXd>(outputValues.data(), outValuesSize.at(0));

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
    std::vector<double> receivedValues           = utils::IntraComm::getCommunication()->receiveRange(0, AsVectorTag<double>{});
    this->output()->data(outputDataID)->values() = Eigen::Map<Eigen::VectorXd>(receivedValues.data(), receivedValues.size());
  }
}
} // namespace mapping
} // namespace precice
