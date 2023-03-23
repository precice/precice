#pragma once

#include <memory>

namespace precice {
namespace impl {

class Participant;
class Coupling;
class WatchPoint;
class WatchIntegral;
struct MeshContext;

using PtrParticipant   = std::shared_ptr<Participant>;
using PtrCoupling      = std::shared_ptr<Coupling>;
using PtrWatchPoint    = std::shared_ptr<WatchPoint>;
using PtrWatchIntegral = std::shared_ptr<WatchIntegral>;
using PtrMeshContext   = std::shared_ptr<MeshContext>;

} // namespace impl
} // namespace precice
