// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_DISTRIBUTEDGEOMETRY_HPP_
#define PRECICE_GEOMETRY_DISTRIBUTEDGEOMETRY_HPP_

#include "CommunicatedGeometry.hpp"
#include "com/SharedPointer.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include <string>
#include <map>

namespace precice {
namespace geometry {

/**
 * @brief Creates the mesh by copying another remote mesh via communication. In contrast to CommunicatedGeometry
 * the mesh is distributed amongst several processors.
 */
class DistributedGeometry : public CommunicatedGeometry
{
public:

  DistributedGeometry (
    const utils::DynVector& offset,
    const std::string&      accessor,
    const std::string&      provider,
    com::PtrCommunication   masterSlaveCom,
    const int               rank,
    const int               size);

  virtual ~DistributedGeometry() {}


protected:

  /**
   * @brief Is called from Geometry  and transmits the mesh.
   */
  virtual void specializedCreate ( mesh::Mesh& seed );

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  com::PtrCommunication _masterSlaveCommunication;

  int _rank;

  int _size;
};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_DISTRIBUTEDGEOMETRY_HPP_ */
