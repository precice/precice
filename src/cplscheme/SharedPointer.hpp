#pragma once

#include <map>
#include <memory>

namespace precice {
namespace cplscheme {

class CouplingScheme;
class CouplingSchemeConfiguration;
class CouplingData;
class GlobalCouplingData;

using PtrCouplingScheme              = std::shared_ptr<CouplingScheme>;
using PtrCouplingSchemeConfiguration = std::shared_ptr<CouplingSchemeConfiguration>;
using PtrCouplingData                = std::shared_ptr<CouplingData>;
using PtrGlobalCouplingData          = std::shared_ptr<GlobalCouplingData>;
using DataMap                        = std::map<int, PtrCouplingData>; /// Map that links DataID to CouplingData

} // namespace cplscheme
} // namespace precice
