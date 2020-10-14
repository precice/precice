#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace cplscheme {

struct Waveform {
  void addNewWindowData(Eigen::VectorXd data)
  {
    // For extrapolation, treat the initial value as old time windows value
    utils::shiftSetFirst(this->lastTimeWindows, data);
  }

  Eigen::VectorXd extrapolateData(int order, int timeWindows)
  {
    Eigen::VectorXd extrapolatedValue;
    if ((order == 1) || (timeWindows == 2 && order == 2)) { //timesteps is increased before extrapolate is called
      // PRECICE_INFO("Performing first order extrapolation");
      PRECICE_ASSERT(this->lastTimeWindows.cols() > 1);
      extrapolatedValue = this->lastTimeWindows.col(0) * 2.0;          // = 2*x^t
      extrapolatedValue -= this->lastTimeWindows.col(1); // = 2*x^t - x^(t-1)
    } else if (order == 2) {
      // PRECICE_INFO("Performing second order extrapolation");
      PRECICE_ASSERT(this->lastTimeWindows.cols() > 2);
      extrapolatedValue = this->lastTimeWindows.col(0) * 2.5;  // = 2.5*x^t
      extrapolatedValue -= this->lastTimeWindows.col(1) * 2.0; // = 2.5*x^t - 2*x^(t-1)
      extrapolatedValue += this->lastTimeWindows.col(2) * 0.5; // = 2.5*x^t - 2*x^(t-1) + 0.5*x^(t-2)
    } else {
      PRECICE_ASSERT(false, "Extrapolation order is invalid.");
    }
    return extrapolatedValue;
  }

  /// Data values of previous time windows.
  Eigen::MatrixXd lastTimeWindows;
};

struct CouplingData {  // @todo: should be a class from a design standpoint. See https://github.com/precice/precice/pull/865#discussion_r495825098
  /// Returns a reference to the data values.
  Eigen::VectorXd &values()
  {
    PRECICE_ASSERT(data != nullptr);
    return data->values();
  }

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const
  {
    PRECICE_ASSERT(data != nullptr);
    return data->values();
  }

  void storeIteration()
  {
    lastIteration = this->values();
  }

  void extrapolateData(int order, int timeWindows)
  {
    waveform.addNewWindowData(this->values());
    this->values() = waveform.extrapolateData(order, timeWindows);
  }

  void initializeWaveform(int extrapolationOrder)
  {
    waveform.lastTimeWindows = Eigen::MatrixXd::Zero(this->values().size(), extrapolationOrder + 1);
  }

  /// Stores data of this and previous time windows and allows to extrapolate.
  Waveform waveform;

  /// Data values of previous iteration.
  Eigen::VectorXd lastIteration;

  mesh::PtrData data;

  mesh::PtrMesh mesh;

  ///  True, if the data values if this CouplingData requires to be initialized by a participant.
  bool requiresInitialization;

  int getDimensions()
  {
    PRECICE_ASSERT(data != nullptr);
    return data->getDimensions();
  }

  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  CouplingData()
  {
    PRECICE_ASSERT(false);
  }

  CouplingData(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization)
      : data(data),
        mesh(mesh),
        requiresInitialization(requiresInitialization),
        waveform()
  {
    PRECICE_ASSERT(data != nullptr);
    PRECICE_ASSERT(mesh != nullptr);
    PRECICE_ASSERT(mesh.use_count() > 0);
  }
};

} // namespace cplscheme
} // namespace precice
