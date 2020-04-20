#pragma once

#include "Mapping.hpp"
#include "impl/BasisFunctions.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"

#include <Eigen/Core>
#include <Eigen/QR>

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
class RadialBasisFctMapping : public Mapping {
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
      Constraint              constraint,
      int                     dimensions,
      RADIAL_BASIS_FUNCTION_T function,
      bool                    xDead,
      bool                    yDead,
      bool                    zDead);

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map(int inputDataID, int outputDataID) override;

  virtual void tagMeshFirstRound() override;

  virtual void tagMeshSecondRound() override;

private:
  precice::logging::Logger _log{"mapping::RadialBasisFctMapping"};

  bool _hasComputedMapping = false;

  /// Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  Eigen::MatrixXd _matrixA;

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> _qr;

  /// true if the mapping along some axis should be ignored
  std::vector<bool> _deadAxis;

  /// Deletes all dead directions from fullVector and returns a vector of reduced dimensionality.
  Eigen::VectorXd reduceVector(const Eigen::VectorXd &fullVector);

  void setDeadAxis(bool xDead, bool yDead, bool zDead)
  {
    _deadAxis.resize(getDimensions());
    if (getDimensions() == 2) {
      _deadAxis[0] = xDead;
      _deadAxis[1] = yDead;
      PRECICE_CHECK(not(xDead && yDead), "You cannot choose all axes to be dead for a RBF mapping");
      if (zDead)
        PRECICE_WARN("Setting the z-axis to dead on a 2-dimensional problem has no effect.");
    } else if (getDimensions() == 3) {
      _deadAxis[0] = xDead;
      _deadAxis[1] = yDead;
      _deadAxis[2] = zDead;
      PRECICE_CHECK(not(xDead && yDead && zDead), "You cannot choose all axes to be dead for a RBF mapping");
    } else {
      PRECICE_ASSERT(false);
    }
  }
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::RadialBasisFctMapping(
    Constraint              constraint,
    int                     dimensions,
    RADIAL_BASIS_FUNCTION_T function,
    bool                    xDead,
    bool                    yDead,
    bool                    zDead)
    : Mapping(constraint, dimensions),
      _basisFunction(function)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  setDeadAxis(xDead, yDead, zDead);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  PRECICE_TRACE();

  precice::utils::Event e("map.rbf.computeMapping.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PRECICE_ASSERT(input()->getDimensions() == output()->getDimensions(),
                 input()->getDimensions(), output()->getDimensions());
  PRECICE_ASSERT(getDimensions() == output()->getDimensions(),
                 getDimensions(), output()->getDimensions());
  int           dimensions = getDimensions();

  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;
  if (getConstraint() == CONSERVATIVE) {
    inMesh  = output();
    outMesh = input();
  } else {
    inMesh  = input();
    outMesh = output();
  }

  int inputSize      = (int) inMesh->vertices().size();
  int outputSize     = (int) outMesh->vertices().size();
  int deadDimensions = 0;
  for (int d = 0; d < dimensions; d++) {
    if (_deadAxis[d])
      deadDimensions += 1;
  }
  int polyparams = 1 + dimensions - deadDimensions;
  PRECICE_ASSERT(inputSize >= 1 + polyparams, inputSize);
  int             n = inputSize + polyparams; // Add linear polynom degrees
  
  std::vector<int> inputIndices(inputSize, -1);

  if(utils::MasterSlave::isSlave()){
    int numberOfVertices = -1;
    numberOfVertices = (int) inMesh->vertices().size();
    utils::MasterSlave::_communication::send(numberOfVertices, 0);
    if(numberOfVertices != 0){
      std::vector<int> vertexIDs(numberOfVertices, -1);
      for(int i = 0; i < numberOfVertices; ++i){
        vertexIDs[i] = inMesh->vertices()[i].getGlobalIndex();
      }
      utils::MasterSlave::_communication.send(vertexIDs, 0);
    }
  } else {
    int numberOfVertices = (int) inMesh->vertices().size();
    for(int i = 0; i < numberOfVertices; ++i){
      inputIndices.push_back(inMesh->vertices()[i].getGlobalIndex());
    }
    for(int slaveRank = 1; slaveRank < utils::MasterSlave::getSize(); ++slaveRank){
      int slaveNumberOfVertices = -1;
      utils::MasterSlave::_communication::receive(slaveNumberOfVertices, slaveRank);
      if(slaveNumberOfVertices != 0{
        std::vector<int> slaveIndices(slaveNumberOfVertices, -1);
        utils::MasterSlave::_communication::receive(slaveIndices, slaveRank);
        for(int i = 0; i < slaveNumberOfVertices; ++i){
          inputIndices.push_back(slaveIndices[i]);
        }
      }
    }
  }
  
  Eigen::MatrixXd matrixCLU(n, n);
  matrixCLU.setZero();
  _matrixA = Eigen::MatrixXd(outputSize, n);
  _matrixA.setZero();

  // Fill upper right part (due to symmetry) of _matrixCLU with values
  int             i = 0;
  Eigen::VectorXd difference(dimensions);
  for (const mesh::Vertex &iVertex : inMesh->vertices()) {
    for (int j = iVertex.getID(); j < inputSize; j++) {
      difference = iVertex.getCoords();
      difference -= inMesh->vertices()[j].getCoords();
      matrixCLU(i, j) = _basisFunction.evaluate(reduceVector(difference).norm());
    }
    matrixCLU(i, inputSize) = 1.0;
    for (int dim = 0; dim < dimensions - deadDimensions; dim++) {
      matrixCLU(i, inputSize + 1 + dim) = reduceVector(iVertex.getCoords())[dim];
    }
    i++;
  }
  // Copy values of upper right part of C to lower left part
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      matrixCLU(j, i) = matrixCLU(i, j);
    }
  }

  // Fill _matrixA with values
  i = 0;
  for (const mesh::Vertex &iVertex : outMesh->vertices()) {
    int j = 0;
    for (const mesh::Vertex &jVertex : inMesh->vertices()) {
      difference = iVertex.getCoords();
      difference -= jVertex.getCoords();
      _matrixA(i, j) = _basisFunction.evaluate(reduceVector(difference).norm());
      j++;
    }
    _matrixA(i, inputSize) = 1.0;
    for (int dim = 0; dim < dimensions - deadDimensions; dim++) {
      _matrixA(i, inputSize + 1 + dim) = reduceVector(iVertex.getCoords())[dim];
    }
    i++;
  }

  _qr = matrixCLU.colPivHouseholderQr();
  if (not _qr.isInvertible()) {
    PRECICE_ERROR("RBF interpolation matrix is not invertible! "
                  "Try to fix axis-aligned mapping setups by marking perpendicular axes as dead.");
  }

  _hasComputedMapping = true;
}

template <typename RADIAL_BASIS_FUNCTION_T>
bool RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::hasComputedMapping() const
{
  return _hasComputedMapping;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  _matrixA            = Eigen::MatrixXd();
  _qr                 = Eigen::ColPivHouseholderQR<Eigen::MatrixXd>();
  _hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::map(
    int inputDataID,
    int outputDataID)
{
  PRECICE_TRACE(inputDataID, outputDataID);

  precice::utils::Event e("map.rbf.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(input()->getDimensions() == output()->getDimensions(),
                 input()->getDimensions(), output()->getDimensions());
  PRECICE_ASSERT(getDimensions() == output()->getDimensions(),
                 getDimensions(), output()->getDimensions());

  Eigen::VectorXd &inValues  = input()->data(inputDataID)->values();
  Eigen::VectorXd &outValues = output()->data(outputDataID)->values();
  int              valueDim  = input()->data(inputDataID)->getDimensions();
  PRECICE_ASSERT(valueDim == output()->data(outputDataID)->getDimensions(),
                 valueDim, output()->data(outputDataID)->getDimensions());
  int deadDimensions = 0;
  for (int d = 0; d < getDimensions(); d++) {
    if (_deadAxis[d])
      deadDimensions += 1;
  }
  int polyparams = 1 + getDimensions() - deadDimensions;

  if (getConstraint() == CONSERVATIVE) {
    PRECICE_DEBUG("Map conservative");
    static int      mappingIndex = 0;
    Eigen::VectorXd Au(_matrixA.cols());  // rows == n
    Eigen::VectorXd in(_matrixA.rows());  // rows == outputSize
    Eigen::VectorXd out(_matrixA.cols()); // rows == n

    // PRECICE_DEBUG("C rows=" << _matrixCLU.rows() << " cols=" << _matrixCLU.cols());
    PRECICE_DEBUG("A rows=" << _matrixA.rows() << " cols=" << _matrixA.cols());
    PRECICE_DEBUG("in size=" << in.size() << ", out size=" << out.size());

    for (int dim = 0; dim < valueDim; dim++) {
      for (int i = 0; i < in.size(); i++) { // Fill input data values
        in[i] = inValues(i * valueDim + dim);
      }

      Au  = _matrixA.transpose() * in;
      out = _qr.solve(Au);

      // Copy mapped data to output data values
      for (int i = 0; i < out.size() - polyparams; i++) {
        outValues(i * valueDim + dim) = out[i];
      }
    }
    mappingIndex++;
  } else { // Map consistent
    PRECICE_DEBUG("Map consistent");
    Eigen::VectorXd p(_matrixA.cols());   // rows == n
    Eigen::VectorXd in(_matrixA.cols());  // rows == n
    Eigen::VectorXd out(_matrixA.rows()); // rows == outputSize
    in.setZero();

    // For every data dimension, perform mapping
    for (int dim = 0; dim < valueDim; dim++) {
      // Fill input from input data values (last polyparams entries remain zero)
      for (int i = 0; i < in.size() - polyparams; i++) {
        in[i] = inValues(i * valueDim + dim);
      }

      p   = _qr.solve(in);
      out = _matrixA * p;

      // Copy mapped data to ouptut data values
      for (int i = 0; i < out.size(); i++) {
        outValues(i * valueDim + dim) = out[i];
      }
    }
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
Eigen::VectorXd RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::reduceVector(
    const Eigen::VectorXd &fullVector)
{
  int deadDimensions = 0;
  for (int d = 0; d < getDimensions(); d++) {
    if (_deadAxis[d])
      deadDimensions += 1;
  }
  PRECICE_ASSERT(getDimensions() > deadDimensions, getDimensions(), deadDimensions);
  Eigen::VectorXd reducedVector(getDimensions() - deadDimensions);
  int             k = 0;
  for (int d = 0; d < getDimensions(); d++) {
    if (not _deadAxis[d]) {
      reducedVector[k] = fullVector[d];
      k++;
    }
  }
  return reducedVector;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshFirstRound()
{
  PRECICE_CHECK(not utils::MasterSlave::isSlave() && not utils::MasterSlave::isMaster(),
                "RBF mapping is not supported for a participant in master mode, use petrbf instead");
}

template <typename RADIAL_BASIS_FUNCTION_T>
void RadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  PRECICE_CHECK(not utils::MasterSlave::isSlave() && not utils::MasterSlave::isMaster(),
                "RBF mapping is not supported for a participant in master mode, use petrbf instead");
}

} // namespace mapping
} // namespace precice
