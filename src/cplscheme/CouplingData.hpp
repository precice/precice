#pragma once

#include <Eigen/Core>
#include <vector>
#include "cplscheme/CouplingScheme.hpp"
#include "mesh/SharedPointer.hpp"
#include "time/Waveform.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme {

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

  int nVertices() const;

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

  /// Returns a const reference to the gradient data values.
  const Eigen::MatrixXd &gradients() const;

  /// Returns number of rows of the stored gradients.
  int gradientsRows() const;

  /// Returns number of columns of the stored gradients.
  int gradientsCols() const;

  /// Returns a const reference to the data Sample.
  const time::Sample &sample() const;

  /// Returns a reference to the time step storage of the data.
  time::Waveform &waveform();

  /// returns previous data interpolated to the relativeDt time
  time::SampleResult getPreviousValuesAtTime(double relativeDt);

  Eigen::MatrixXd getPreviousGradientsAtTime(double relativeDt);

  /// Returns a const reference to the time step storage of the data.
  const time::Waveform &waveform() const;

  /// Returns the stamples in the Waveform
  auto stamples() const
  {
    return waveform().stamples();
  }

  /// Add sample at given time to the Waveform
  void setSampleAtTime(double time, time::Sample sample);

  /// Set _data::_sample
  void setGlobalSample(const time::Sample &sample); // @todo try to remove this function

  /// Add sample with zero values at given time to the Waveform
  void initializeWithZeroAtTime(double time);

  /// Creates an empty sample at given time
  void emplaceSampleAtTime(double time);

  /// Creates a sample at given time with given values
  void emplaceSampleAtTime(double time, std::initializer_list<double> values);

  /// Creates a sample at given time with given values and gradients
  void emplaceSampleAtTime(double time, std::initializer_list<double> values, std::initializer_list<double> gradients);

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

  std::vector<double> getLowerBound() const;

  std::vector<double> getUpperBound() const;

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
  time::Waveform _previousTimeStepsStorage;

  /// If true, all substeps will be sent / received for this coupling data
  bool _exchangeSubsteps;

  Direction _direction;
};

} // namespace precice::cplscheme
