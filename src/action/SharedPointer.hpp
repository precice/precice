#pragma once

#include <memory>

namespace precice::action {

class Action;
class ActionConfiguration;

using PtrAction              = std::unique_ptr<Action>;
using PtrActionConfiguration = std::shared_ptr<ActionConfiguration>;

} // namespace precice::action
