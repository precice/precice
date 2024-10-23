#pragma once

#include <Eigen/Core>
#include <numeric>
#include <vector>

#include "cplscheme/CouplingData.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace acceleration {
namespace impl {

/**
 * @brief Interface for storing the time grids in the Quasi-Newton and Aitken methods
 *
 *
 * initialize saves the time grids of the waveforms
 * moveTimeGridToNewWindow Moves the time grid to the new time window. This is done by translating and scaling the QN time grid such that its startpoint and end point is the same as the waveforms
 *
 */
class WaveformTimeGrids {
public:
  /// Map from data ID to data values.
  using DataMap = std::map<int, cplscheme::PtrCouplingData>;

  /**
   * @brief sets the time grid
   * @param cplData
   * @param dataIDs
   * @param reduced if true then only the last timestep is saved in the timeGrid otherwise all timesteps are saved in the time grid
   * This is used in the Quasi-Newton method to toggle between the full and reduced variant detailed here
   */
  void setTimeGrid(const DataMap &cplData, std::vector<int> dataIDs, bool reduced);

  Eigen::VectorXd getTimeGrid(int dataID);

  void moveTimeGridToNewWindow(const DataMap &cplData);

private:
  logging::Logger _log{"acceleration::WaveformTimeGrids"};

  /// @brief List of the time grid to which all the data will be interpolated to
  /// Stored in a map, since each data entry has its own time grid
  std::map<int, Eigen::VectorXd> _timeGrids;

  //These two are needed for the reduced variants, since then the _primaryTimeGrids only contain one timestep
  std::map<int, double> _timeGridsStartTime;
  std::map<int, double> _timeGridsEndTime;
};

} // namespace impl
} // namespace acceleration
} // namespace precice
