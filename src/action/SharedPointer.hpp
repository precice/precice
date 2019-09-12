#pragma once

#include <memory>

namespace precice {
namespace action {

class Action;
class ActionConfiguration;

using PtrAction              = std::shared_ptr<Action>;
using PtrActionConfiguration = std::shared_ptr<ActionConfiguration>;

} // namespace action
} // namespace precice
