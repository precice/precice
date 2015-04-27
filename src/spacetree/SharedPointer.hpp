#pragma once

#include <memory>

namespace precice {
namespace spacetree {

class Spacetree;
class SpacetreeConfiguration;

using PtrSpacetree              = std::shared_ptr<Spacetree>;
using PtrSpacetreeConfiguration = std::shared_ptr<SpacetreeConfiguration>;

}} // namespace precice, spacetree
