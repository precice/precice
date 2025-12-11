#pragma once

#include <string>
#include <vector>
#include "Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::action {

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
      Timing                  timing,
      const std::vector<int> &sourceDataIDs,
      int                     targetDataID,
      const mesh::PtrMesh    &mesh);

  /// Adding data and applying them to target
  void performAction() final override;

private:
  logging::Logger _log{"action::SummationAction"};

  mesh::PtrData              _targetData;
  std::vector<mesh::PtrData> _sourceDataVector;
};

} // namespace precice::action
