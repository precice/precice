#include <Eigen/Core>
#include <numeric>
#include <vector>

#include "TimeGrids.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

TimeGrids::TimeGrids(const DataMap &cplData, std::vector<int> dataIDs, bool reducedTimeGrid)
{
  for (int dataID : dataIDs) {
    Eigen::VectorXd timeGrid = cplData.at(dataID)->waveform().getTimes();
    if (reducedTimeGrid) {
      _timeGrids.insert(std::pair<int, Eigen::VectorXd>(dataID, timeGrid.tail<1>()));
    } else {
      _timeGrids.insert(std::pair<int, Eigen::VectorXd>(dataID, timeGrid));
    }
  }
}

Eigen::VectorXd TimeGrids::getTimeGridAfter(int dataID, double time) const
{
  PRECICE_ASSERT(_timeGrids.count(dataID), "there does not exists a stored time grid corresponding to this dataID");

  std::vector<double> reduced;
  for (double d : _timeGrids.at(dataID)) {
    if (math::greater(d, time)) {
      reduced.push_back(d);
    }
  }

  return Eigen::Map<Eigen::VectorXd>(reduced.data(), reduced.size());
}

void TimeGrids::moveTimeGridToNewWindow(const DataMap &cplData)
{
  for (auto &pair : _timeGrids) {
    if (pair.second.size() == 1) {
      _timeGrids.at(pair.first) = cplData.at(pair.first)->waveform().getTimes().tail<1>();
    } else {
      int dataID = pair.first;
      // Only way to access the first time stamp is through the whole vector
      Eigen::VectorXd newtimeGrid = cplData.at(dataID)->waveform().getTimes();
      double          newTimesMin = newtimeGrid(0);
      double          newTimesMax = newtimeGrid(newtimeGrid.size() - 1);

      Eigen::VectorXd timeGrid    = pair.second;
      double          oldTimesMin = _timeGrids.at(dataID)(0);
      double          oldTimesMax = _timeGrids.at(dataID)(timeGrid.size() - 1);

      //  Linearly transforms the time grid from the old time window [t_{N-1}, t_N] to the new time window [t_N, t_{N+1}] by translating and scaling the time grid
      auto transformNewTime = [oldTimesMin, oldTimesMax, newTimesMin, newTimesMax](double t) -> double { return (t - oldTimesMin) / (oldTimesMax - oldTimesMin) * (newTimesMax - newTimesMin) + newTimesMin; };
      timeGrid              = timeGrid.unaryExpr(transformNewTime);
      _timeGrids.at(dataID) = timeGrid;
    }
  }
}
} // namespace precice::time
