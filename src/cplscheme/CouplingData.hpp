#pragma once

#include <Eigen/Core>
#include <vector>
#include "cplscheme/CouplingScheme.hpp"
#include "mesh/SharedPointer.hpp"
#include "time/Storage.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

class CouplingData {
public:
  CouplingData(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization,
      int           extrapolationOrder);

  int getDimensions() const;

  int getSize() const;

  /// Returns a reference to the data values.
  Eigen::VectorXd &values();

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

  /// Returns a reference to the gradient data values.
  Eigen::MatrixXd &gradients();

  /// Returns a const reference to the gradient data values.
  const Eigen::MatrixXd &gradients() const;

  /// Returns a reference to the gradient data Sample.
  time::Sample &sample();

  /// Returns a const reference to the data Sample.
  const time::Sample &sample() const;

  /// Returns a reference to the time step storage of the data.
  time::Storage &timeStepsStorage();

  /// Returns a const reference to the time step storage of the data.
  const time::Storage &timeStepsStorage() const;

  /// Returns the stamples in _timeStepsStorage.
  auto stamples() const
  {
    return timeStepsStorage().stamples();
  }

  /// Add sample at given time to _timeStepsStorage.
  void setSampleAtTime(double time, time::Sample sample);

  /// Returns if the data contains gradient data
  bool hasGradient() const;

  /// Returns the dimensions of the current mesh (2D or 3D)
  int meshDimensions() const;

  /// store _data->values() in read-only variable _previousIteration for convergence checks etc.
  void storeIteration();

  /// returns data value from previous iteration
  const Eigen::VectorXd previousIteration() const;

  /// returns gradient data from previous iteration
  const Eigen::MatrixXd &previousIterationGradients() const;

  /// returns size of previous iteration
  int getPreviousIterationSize() const;

  /// get ID of this CouplingData's mesh. See Mesh::getID().
  int getMeshID();

  /// get ID of this CouplingData's data. See Data::getID().
  int getDataID();

  /// get name of this CouplingData's data. See Data::getName().
  std::string getDataName();

  /// get vertex offsets of this CouplingData's mesh. See Mesh::getVertexOffsets().
  std::vector<int> getVertexOffsets();

  ///  True, if the data values of this CouplingData require to be initialized by this participant.
  const bool requiresInitialization;

  /// returns keys in _timeStepsStorage in ascending order.
  Eigen::VectorXd getStoredTimesAscending();

  /// move to next window and initialize data via extrapolation
  void moveToNextWindow();

  /**
   * @brief Returns the values of all time steps stored in this coupling data in a serialized fashion
   *
   * Serialization of the data is performed per mesh node. The dimension of one node is then data dimension (scalar or vector) * number of time steps.
   *
   * @return Eigen::VectorXd a vector containing all data for all time steps in serialized fashion.
   */
  Eigen::VectorXd getSerializedValues();

  /**
   * @brief Returns the gradients of all time steps stored in this coupling data in a serialized fashion
   *
   * Serialization of the data is performed per mesh node. The dimension of one node is then data dimension (scalar or vector) * number of time steps.
   *
   * @return Eigen::VectorXd a vector containing all gradient data for all time steps in serialized fashion.
   */
  Eigen::VectorXd getSerializedGradients();

  /**
   * @brief accepts serialized values and stores them in this coupling data
   *
   * @param timesAscending times associated with data
   * @param serializedValues values in serialized form
   */
  void storeFromSerialized(Eigen::VectorXd timesAscending, Eigen::VectorXd serializedValues);

  /**
   * @brief accepts serialized values and gradients and stores them in this coupling data
   *
   * @param timesAscending times associated with data
   * @param serializedValues values in serialized form
   * @param serializedGradients gradients in serialized form
   */
  void storeFromSerialized(Eigen::VectorXd timesAscending, Eigen::VectorXd serializedValues, Eigen::MatrixXd serializedGradients);

private:
  logging::Logger _log{"cplscheme::CouplingData"};

  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  CouplingData()
      : requiresInitialization(false)
  {
    PRECICE_ASSERT(false);
  }

  /// Sample values of previous iteration (end of time window).
  time::Sample _previousIteration;

  /// Data associated with this CouplingData
  mesh::PtrData _data;

  /// Mesh associated with this CouplingData
  mesh::PtrMesh _mesh;
};

} // namespace cplscheme
} // namespace precice
