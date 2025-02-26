#pragma once

#include <map>
#include <memory>

namespace precice::cplscheme {

class CouplingScheme;
class CouplingSchemeConfiguration;
class CouplingData;

using PtrCouplingScheme              = std::shared_ptr<CouplingScheme>;
using PtrCouplingSchemeConfiguration = std::shared_ptr<CouplingSchemeConfiguration>;
using PtrCouplingData                = std::shared_ptr<CouplingData>;
using DataMap                        = std::map<int, PtrCouplingData>; /// Map that links DataID to CouplingData

} // namespace precice::cplscheme
