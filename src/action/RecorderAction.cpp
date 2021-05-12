#include "action/RecorderAction.hpp"
#include <vector>
#include "action/Action.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace action {

RecorderAction::RecorderAction(
    Timing               timing,
    const mesh::PtrMesh &mesh)
    : Action(timing, mesh) {}

void RecorderAction::performAction(
    double time,
    double timeStepSize,
    double computedTimeWindowPart,
    double timeWindowSize)
{
  records.push_back(Record{
      getTiming(), time, timeStepSize, computedTimeWindowPart, timeWindowSize});
}

std::vector<RecorderAction::Record> RecorderAction::records{};

void RecorderAction::reset()
{
  records.clear();
}

} // namespace action
} // namespace precice
