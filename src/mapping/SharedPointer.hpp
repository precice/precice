#pragma once

#include <memory>

namespace precice::mapping {

class Mapping;
class MappingConfiguration;

using PtrMapping              = std::shared_ptr<Mapping>;
using PtrMappingConfiguration = std::shared_ptr<MappingConfiguration>;

} // namespace precice::mapping
