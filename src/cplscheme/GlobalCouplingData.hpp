#pragma once

#include <Eigen/Core>
#include <vector>
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/impl/Extrapolation.hpp"
#include "mesh/SharedPointer.hpp"
#include "time/Storage.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

class GlobalCouplingData {
public:
  GlobalCouplingData(
      mesh::PtrData data,
      bool          requiresInitialization,
      int           extrapolationOrder = CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER);

  int getDimensions() const;

  int getSize() const;

  /// Returns a reference to the data values.
  Eigen::VectorXd &values();

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

  /// Returns a reference to the data Sample.
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

  /// store _globalData->values() in read-only variable _previousIteration for convergence checks etc.
  void storeIteration();

  /// returns data value from previous iteration
  const Eigen::VectorXd previousIteration() const;

  /// returns size of previous iteration
  int getPreviousIterationSize() const;

  /// get ID of this GlobalCouplingData's data. See GlobalData::getID().
  int getDataID();

  /// get name of this GlobalCouplingData's data. See GlobalData::getName().
  std::string getDataName();

  ///  True, if the data values of this CouplingData require to be initialized by this participant.
  const bool requiresInitialization;

  /// initialize _extrapolation
  void initializeExtrapolation();

  /// move to next window and initialize data via extrapolation
  void moveToNextWindow();

  /// store current value in _extrapolation
  void storeExtrapolationData();

private:
  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  GlobalCouplingData()
      : requiresInitialization(false),
        _extrapolation(CouplingScheme::UNDEFINED_EXTRAPOLATION_ORDER)
  {
    PRECICE_ASSERT(false);
  }

  /// Sample values of previous iteration (end of time window).
  time::Sample _previousIteration;

  /// Data associated with this CouplingData
  mesh::PtrData _data;

  /// Extrapolation associated with this CouplingData
  cplscheme::impl::Extrapolation _extrapolation;
};

} // namespace cplscheme
} // namespace precice
