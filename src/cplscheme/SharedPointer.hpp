#pragma once

#include <memory>

namespace precice {
namespace time {

class Waveform;
using PtrWaveform = std::shared_ptr<Waveform>;

} // namespace time
} // namespace precice

namespace precice {
namespace cplscheme {

class CouplingScheme;
class CouplingSchemeConfiguration;
class CouplingData;

using PtrCouplingScheme              = std::shared_ptr<CouplingScheme>;
using PtrCouplingSchemeConfiguration = std::shared_ptr<CouplingSchemeConfiguration>;
using PtrCouplingData                = std::shared_ptr<CouplingData>;

} // namespace cplscheme
} // namespace precice
