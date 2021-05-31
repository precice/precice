#pragma once

#include <memory>

namespace precice {
namespace cplscheme {

class CouplingScheme;
class CouplingSchemeConfiguration;
struct CouplingData;

using PtrCouplingScheme              = std::shared_ptr<CouplingScheme>;
using PtrCouplingSchemeConfiguration = std::shared_ptr<CouplingSchemeConfiguration>;
using PtrCouplingData                = std::shared_ptr<CouplingData>;

} // namespace cplscheme
} // namespace precice
