#pragma once

#include <memory>

namespace precice {
namespace cplscheme {

class CouplingScheme;
class CouplingSchemeConfiguration;
class PostProcessingConfiguration;
struct CouplingData;

using PtrCouplingScheme              = std::shared_ptr<CouplingScheme>;
using PtrCouplingSchemeConfiguration = std::shared_ptr<CouplingSchemeConfiguration>;
using PtrPostProcessingConfiguration = std::shared_ptr<PostProcessingConfiguration>;
using PtrCouplingData                = std::shared_ptr<CouplingData>;

}} // namespace precice, cplscheme
