#pragma once

#include <string>
#include <vector>
#include "Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace action {

/// Action that adds multiple source data into target data
class SummationAction : public Action {
public:
  /**
	 * @brief Constructor
	 * 
	 * @param[in] Timing When to apply the action
	 * @param[in] sourceDataIDs Data indexes which are to be added
	 * @param[in] targetDataID Data in which the action will be applied
	 * 
	 */
  SummationAction(
      Timing               timing,
      std::vector<int>     sourceDataIDs,
      int                  targetDataID,
      const mesh::PtrMesh &mesh);

  virtual ~SummationAction() {}

  /// Adding data and applying them to target
  virtual void performAction(
      double time,
      double dt,
      double computedPartFullDt,
      double fullDt);

private:
  logging::Logger _log{"action::SummationAction"};

  mesh::PtrData              _targetData;
  std::vector<mesh::PtrData> _sourceDataVector;
};

} // namespace action
} // namespace precice