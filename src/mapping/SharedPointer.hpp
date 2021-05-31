#pragma once

#include <memory>

namespace precice {
namespace mapping {

class Mapping;
class MappingConfiguration;

using PtrMapping              = std::shared_ptr<Mapping>;
using PtrMappingConfiguration = std::shared_ptr<MappingConfiguration>;

} // namespace mapping
} // namespace precice
