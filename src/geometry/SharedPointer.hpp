#pragma once

#include <memory>

namespace precice {
namespace geometry {

class Geometry;
class CustomGeometry;
class GeometryConfiguration;

using PtrGeometry              = std::shared_ptr<Geometry>;
using PtrCustomGeometry        = std::shared_ptr<CustomGeometry>;
using PtrGeometryConfiguration = std::shared_ptr<GeometryConfiguration>;

}} // namespace precice, geometry
