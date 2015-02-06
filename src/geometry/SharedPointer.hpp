// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_SHAREDPOINTER_HPP_
#define PRECICE_GEOMETRY_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace geometry {

class Geometry;
class CustomGeometry;
class GeometryConfiguration;

typedef boost::shared_ptr<Geometry>      PtrGeometry;
typedef boost::shared_ptr<CustomGeometry>        PtrCustomGeometry;
typedef boost::shared_ptr<GeometryConfiguration> PtrGeometryConfiguration;


}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_SHAREDPOINTER_HPP_ */
