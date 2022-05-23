#pragma once

#include <Eigen/Core>
#include <Eigen/QR>

#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "impl/BasisFunctions.hpp"
#include "mapping/RadialBasisFctBaseMapping.hpp"
#include "mesh/Filter.hpp"
#include "precice/types.hpp"
#include "utils/EigenHelperFunctions.hpp"
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
      Mapping::Constraint     constraint,
      int                     dimensions,
      RADIAL_BASIS_FUNCTION_T function,
      bool                    xDead,
      bool                    yDead,
      bool                    zDead);

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Removes a computed mapping.
  virtual void clear() override;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctMapping"};

  Eigen::MatrixXd _matrixA;

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qr;

  /// @copydoc RadialBasisFctBaseMapping::mapConservative
  virtual void mapConservative(DataID inputDataID, DataID outputDataID) override;

  /// @copydoc RadialBasisFctBaseMapping::mapConsistent
  virtual void mapConsistent(DataID inputDataID, DataID outputDataID) override;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctMapping(
    Mapping::Constraint     constraint,
    int                     dimensions,
    RADIAL_BASIS_FUNCTION_T function,
    bool                    xDead,
    bool                    yDead,
    bool                    zDead)
    : RadialBasisFctBaseMapping<RADIAL_BASIS_FUNCTION_T>(constraint, dimensions, function, xDead, yDead, zDead)
{
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
    com::CommunicateMesh(utils::IntraComm::getCommunication()).sendMesh(filteredInMesh, 0);
    com::CommunicateMesh(utils::IntraComm::getCommunication()).sendMesh(*outMesh, 0);

  } else { // Parallel Primary rank or Serial

    mesh::Mesh globalInMesh("globalInMesh", inMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
    mesh::Mesh globalOutMesh("globalOutMesh", outMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);

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
        com::CommunicateMesh(utils::IntraComm::getCommunication()).receiveMesh(secondaryInMesh, secondaryRank);
        globalInMesh.addMesh(secondaryInMesh);

        mesh::Mesh secondaryOutMesh(outMesh->getName(), outMesh->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
        com::CommunicateMesh(utils::IntraComm::getCommunication()).receiveMesh(secondaryOutMesh, secondaryRank);
        globalOutMesh.addMesh(secondaryOutMesh);
      }

    } else { // Serial
      globalInMesh.addMesh(*inMesh);
      globalOutMesh.addMesh(*outMesh);
    }

    _matrixA = buildMatrixA(this->_basisFunction, globalInMesh, globalOutMesh, this->_deadAxis);
    _qr      = buildMatrixCLU(this->_basisFunction, globalInMesh, this->_deadAxis).colPivHouseholderQr();

    PRECICE_CHECK(_qr.isInvertible(),
                  "The interpolation matrix of the RBF mapping from mesh {} to mesh {} is not invertable. "
                  "This means that the mapping problem is not well-posed. "
                  "Please check if your coupling meshes are correct. Maybe you need to fix axis-aligned mapping setups "
                  "by marking perpendicular axes as dead?",
                  this->input()->getName(), this->output()->getName());
  }
  this->_hasComputedMapping = true;
  PRECICE_DEBUG("Compute Mapping is Completed.");
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  _matrixA                  = Eigen::MatrixXd();
  _qr                       = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>();
  this->_hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::mapConservative(DataID inputDataID, DataID outputDataID)
{
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
    Eigen::VectorXd             outputValues((_matrixA.cols() - this->getPolynomialParameters()) * valueDim);
    outputValues.setZero();

    Eigen::VectorXd Au(_matrixA.cols());  // rows == n
    Eigen::VectorXd in(_matrixA.rows());  // rows == outputSize
    Eigen::VectorXd out(_matrixA.cols()); // rows == n

    for (int dim = 0; dim < valueDim; dim++) {
      for (int i = 0; i < in.size(); i++) { // Fill input data values
        in[i] = inputValues(i * valueDim + dim);
      }

      Au  = _matrixA.transpose() * in;
      out = _qr.solve(Au);

      // Copy mapped data to output data values
      for (int i = 0; i < out.size() - this->getPolynomialParameters(); i++) {
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

    std::vector<double> globalInValues((_matrixA.cols() - this->getPolynomialParameters()) * valueDim, 0.0);
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

    Eigen::VectorXd p(_matrixA.cols());   // rows == n
    Eigen::VectorXd in(_matrixA.cols());  // rows == n
    Eigen::VectorXd out(_matrixA.rows()); // rows == outputSize
    in.setZero();

    // Construct Eigen vectors
    Eigen::Map<Eigen::VectorXd> inputValues(globalInValues.data(), globalInValues.size());

    Eigen::VectorXd outputValues((_matrixA.rows()) * valueDim);
    outputValues.setZero();

    // For every data dimension, perform mapping
    for (int dim = 0; dim < valueDim; dim++) {
      // Fill input from input data values (last polyparams entries remain zero)
      for (int i = 0; i < in.size() - this->getPolynomialParameters(); i++) {
        in[i] = inputValues[i * valueDim + dim];
      }

      p   = _qr.solve(in);
      out = _matrixA * p;

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

// ------- Non-Member Functions ---------

template <typename RADIAL_BASIS_FUNCTION_T>
static Eigen::MatrixXd buildMatrixCLU(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, std::vector<bool> deadAxis)
{
  int inputSize  = inputMesh.vertices().size();
  int dimensions = inputMesh.getDimensions();

  int deadDimensions = 0;
  for (int d = 0; d < dimensions; d++) {
    if (deadAxis[d])
      deadDimensions += 1;
  }

  int polyparams = 1 + dimensions - deadDimensions;
  PRECICE_ASSERT(inputSize >= 1 + polyparams, inputSize);
  int n = inputSize + polyparams; // Add linear polynom degrees

  Eigen::MatrixXd matrixCLU(n, n);
  matrixCLU.setZero();

  for (int i = 0; i < inputSize; ++i) {
    for (int j = i; j < inputSize; ++j) {
      const auto &u   = inputMesh.vertices()[i].getCoords();
      const auto &v   = inputMesh.vertices()[j].getCoords();
      matrixCLU(i, j) = basisFunction.evaluate(utils::reduceVector((u - v), deadAxis).norm());
    }

    const auto reduced = utils::reduceVector(inputMesh.vertices()[i].getCoords(), deadAxis);

    for (int dim = 0; dim < dimensions - deadDimensions; dim++) {
      matrixCLU(i, inputSize + 1 + dim) = reduced[dim];
    }
    matrixCLU(i, inputSize) = 1.0;
  }

  matrixCLU.triangularView<Eigen::Lower>() = matrixCLU.transpose();

  return matrixCLU;
}

template <typename RADIAL_BASIS_FUNCTION_T>
static Eigen::MatrixXd buildMatrixA(RADIAL_BASIS_FUNCTION_T basisFunction, const mesh::Mesh &inputMesh, const mesh::Mesh &outputMesh, std::vector<bool> deadAxis)
{
  int inputSize  = inputMesh.vertices().size();
  int outputSize = outputMesh.vertices().size();
  int dimensions = inputMesh.getDimensions();

  int deadDimensions = 0;
  for (int d = 0; d < dimensions; d++) {
    if (deadAxis[d])
      deadDimensions += 1;
  }

  int polyparams = 1 + dimensions - deadDimensions;
  PRECICE_ASSERT(inputSize >= 1 + polyparams, inputSize);
  int n = inputSize + polyparams; // Add linear polynom degrees

  Eigen::MatrixXd matrixA(outputSize, n);
  matrixA.setZero();

  // Fill _matrixA with values
  for (int i = 0; i < outputSize; ++i) {
    for (int j = 0; j < inputSize; ++j) {
      const auto &u = outputMesh.vertices()[i].getCoords();
      const auto &v = inputMesh.vertices()[j].getCoords();
      matrixA(i, j) = basisFunction.evaluate(utils::reduceVector((u - v), deadAxis).norm());
    }

    const auto reduced = utils::reduceVector(outputMesh.vertices()[i].getCoords(), deadAxis);

    for (int dim = 0; dim < dimensions - deadDimensions; dim++) {
      matrixA(i, inputSize + 1 + dim) = reduced[dim];
    }
    matrixA(i, inputSize) = 1.0;
  }
  return matrixA;
}

} // namespace mapping
} // namespace precice
