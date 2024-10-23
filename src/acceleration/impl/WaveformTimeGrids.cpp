#pragma once

#include <Eigen/Core>
#include <numeric>
#include <vector>

#include "WaveformTimeGrids.hpp"
#include "utils/assertion.hpp"

namespace precice::acceleration::impl {

void WaveformTimeGrids::setTimeGrid(const DataMap &cplData, std::vector<int> dataIDs, bool reduced)
{
  for (int dataID : dataIDs) {
    Eigen::VectorXd timeGrid = cplData.at(dataID)->timeStepsStorage().getTimes();
    if (reduced) {
      _timeGrids.insert(std::pair<int, Eigen::VectorXd>(dataID, timeGrid.tail<1>()));
    } else {
      _timeGrids.insert(std::pair<int, Eigen::VectorXd>(dataID, timeGrid));
    }
    _timeGridsStartTime.insert(std::pair<int, double>(dataID, timeGrid(0)));
    _timeGridsEndTime.insert(std::pair<int, double>(dataID, timeGrid(timeGrid.size() - 1)));
  }
}

Eigen::VectorXd WaveformTimeGrids::getTimeGrid(int dataID)
{

  if (_timeGrids.count(dataID)) {
    PRECICE_ASSERT("there does not exists a stored time grid corresponding to this dataID");
  }

  return _timeGrids.at(dataID);
}

void WaveformTimeGrids::moveTimeGridToNewWindow(const DataMap &cplData)
{
  for (auto &pair : _timeGrids) {
    int dataID = pair.first;
    // Only way to access the first time stamp is through the whole vector
    Eigen::VectorXd newtimeGrid = cplData.at(dataID)->timeStepsStorage().getTimes();
    double          newTimesMin = newtimeGrid(0);
    double          newTimesMax = newtimeGrid(newtimeGrid.size() - 1);

    Eigen::VectorXd timeGrid    = pair.second;
    double          oldTimesMin = _timeGridsStartTime.at(dataID);
    double          oldTimesMax = _timeGridsEndTime.at(dataID);

    // transform the time to the new time grid
    auto transformNewTime = [oldTimesMin = oldTimesMin, oldTimesMax = oldTimesMax, newTimesMin = newTimesMin, newTimesMax = newTimesMax](double t) -> double { return (t - oldTimesMin) / (oldTimesMax - oldTimesMin) * (newTimesMax - newTimesMin) + newTimesMin; };
    timeGrid              = timeGrid.unaryExpr(transformNewTime);
    _timeGrids.at(dataID) = timeGrid;

    _timeGridsStartTime.at(dataID) = newTimesMin;
    _timeGridsEndTime.at(dataID)   = newTimesMax;
  }
}
} // namespace precice::acceleration::impl
