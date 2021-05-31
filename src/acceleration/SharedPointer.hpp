#pragma once

#include <memory>

namespace precice {
namespace acceleration {
class Acceleration;
class AccelerationConfiguration;

using PtrAcceleration              = std::shared_ptr<Acceleration>;
using PtrAccelerationConfiguration = std::shared_ptr<AccelerationConfiguration>;
} // namespace acceleration
} // namespace precice
