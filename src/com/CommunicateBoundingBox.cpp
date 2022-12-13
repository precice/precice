#include <cstddef>
#include <iterator>
#include <memory>
#include <numeric>
#include <utility>

#include "com/CommunicateBoundingBox.hpp"
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice::com {
CommunicateBoundingBox::CommunicateBoundingBox(
    com::PtrCommunication communication)
    : _communication(std::move(communication))
{
}

void CommunicateBoundingBox::sendBoundingBox(
    const mesh::BoundingBox &bb,
    int                      rankReceiver)
{
  PRECICE_TRACE(rankReceiver);
  _communication->sendRange(bb.dataVector(), rankReceiver);
}

void CommunicateBoundingBox::receiveBoundingBox(
    mesh::BoundingBox &bb,
    int                rankSender)
{
  PRECICE_TRACE(rankSender);
  auto              receivedData = _communication->receiveRange(rankSender, AsVectorTag<double>{});
  mesh::BoundingBox tempBB(receivedData);
  bb = std::move(tempBB);
}

} // namespace precice::com
