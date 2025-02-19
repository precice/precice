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
  enum struct Direction : bool { Send,
                                 Receive };

  CouplingData(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization,
      bool          exchangeSubsteps,
      Direction     direction);

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

  /// returns previous data interpolated to the relativeDt time
  Eigen::VectorXd getPreviousValuesAtTime(double relativeDt);

  Eigen::MatrixXd getPreviousGradientsAtTime(double relativeDt);

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

  /// Reshape the past iterations and initial sample during remeshing
  void reinitialize();

  /// store _data->values() in read-only variable _previousIteration for convergence checks etc.
  void storeIteration();

  /// returns data value from previous iteration
  const Eigen::VectorXd &previousIteration() const;

  /// returns gradient data from previous iteration
  const Eigen::MatrixXd &previousIterationGradients() const;

  /// returns size of previous iteration
  int getPreviousIterationSize() const;

  /// get ID of this CouplingData's mesh. See Mesh::getID().
  int getMeshID();

  /// get ID of this CouplingData's data. See Data::getID().
  int getDataID();

  /// get name of this CouplingData's data. See Data::getName().
  std::string getDataName() const;

  /// get name of this CouplingData's mesh. See Mesh::getName().
  std::string getMeshName() const;

  /// get vertex offsets of this CouplingData's mesh. See Mesh::getVertexOffsets().
  std::vector<int> getVertexOffsets();

  /// get direction of this coupling data
  Direction getDirection() const;

  ///  True, if the data values of this CouplingData require to be initialized by this participant.
  const bool requiresInitialization;

  /// move to next window and initialize data via extrapolation
  void moveToNextWindow();

  bool exchangeSubsteps() const;

private:
  logging::Logger _log{"cplscheme::CouplingData"};

  /// Mesh associated with this CouplingData
  mesh::PtrMesh _mesh;

  /// Data associated with this CouplingData
  mesh::PtrData _data;

  /// Sample values of previous iteration (end of time window).
  time::Storage _previousTimeStepsStorage;

  /// If true, all substeps will be sent / received for this coupling data
  bool _exchangeSubsteps;

  Direction _direction;
};

} // namespace cplscheme
} // namespace precice
