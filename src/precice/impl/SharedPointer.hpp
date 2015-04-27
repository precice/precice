#pragma once

#include <memory>

namespace precice {
namespace impl {

class Participant;
class Coupling;
class WatchPoint;
class AbstractDataAction;
struct MeshContext;

using PtrParticipant        = std::shared_ptr<Participant>;
using PtrCoupling           = std::shared_ptr<Coupling>;
using PtrWatchPoint         = std::shared_ptr<WatchPoint>;
using PtrAbstractDataAction = std::shared_ptr<AbstractDataAction>;
using PtrMeshContext        = std::shared_ptr<MeshContext>;

}} // namespace precice, impl
