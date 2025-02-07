#pragma once

#include <Eigen/Core>
#include <numeric>
#include <vector>

#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace time {

/**
 * @brief Interface for storing the time grids in the Quasi-Newton and Aitken methods.
 * A time grid is a ordered vector containing the time stamps from the samples in the waveform of the coupled data.
 *
 *
 */
class TimeGrids {
public:
  /// Map from data ID to data values.
  using DataMap = std::map<int, cplscheme::PtrCouplingData>;

  /**
   * @brief saves the current time grid of the couplig data
   * @param cplData
   * @param dataIDs of which time grids we want to save
   * @param reduced if true then only the last timestep is saved in the timeGrid otherwise all timesteps are saved in the time grid
   * This is used in the Quasi-Newton method to toggle between the full and reduced variant detailed here
   */
  TimeGrids(const DataMap &cplData, std::vector<int> dataIDs, bool reduced);

  Eigen::VectorXd getTimeGrid(int dataID) const;

  Eigen::VectorXd getTimeGridAfter(int dataID, double time) const;

  //  Linearly transforms the time grid from the old time window [t_{N-1}, t_N] to the new time window [t_N, t_{N+1}]
  // This is done to allow the QN methods to sample from the new time window while keeping the structure and vector dimensions inside the QN method
  void moveTimeGridToNewWindow(const DataMap &cplData);

private:
  logging::Logger _log{"acceleration::WaveformTimeGrids"};

  /// @brief List of the time grid to which all the data will be interpolated to
  /// Stored in a map, since each data entry has its own time grid
  std::map<int, Eigen::VectorXd> _timeGrids;
};

} // namespace time
} // namespace precice
