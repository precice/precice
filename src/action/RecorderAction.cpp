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
    double dt,
    double computedPartFullDt,
    double fullDt)
{
  records.push_back(Record{
      getTiming(), time, dt, computedPartFullDt, fullDt});
}

std::vector<RecorderAction::Record> RecorderAction::records{};

void RecorderAction::reset()
{
  records.clear();
}

} // namespace action
} // namespace precice
