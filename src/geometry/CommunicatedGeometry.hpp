// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_COMMUNICATEDGEOMETRY_HPP_
#define PRECICE_GEOMETRY_COMMUNICATEDGEOMETRY_HPP_

#include "Geometry.hpp"
#include "com/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "utils/Dimensions.hpp"
#include "utils/MasterSlave.hpp"
#include "tarch/logging/Log.h"
#include <string>
#include <map>

namespace precice {
namespace geometry {

/**
 * @brief Creates the mesh by copying another remote mesh via communication.
 */
class CommunicatedGeometry : public Geometry
{
public:

  CommunicatedGeometry (
    const utils::DynVector& offset,
    const std::string&      accessor,
    const std::string&      provider,
    int                     dimensions);

  virtual ~CommunicatedGeometry() {}

  void addReceiver (
    const std::string&     receiver,
    com::PtrCommunication com );

  void setBoundingFromMapping(mapping::PtrMapping mapping);

  void setBoundingToMapping(mapping::PtrMapping mapping);

protected:

  /**
   * @brief Is called from Geometry  and transmits the mesh.
   */
  virtual void specializedCreate ( mesh::Mesh& seed );

private:

  void gatherMesh(
    mesh::Mesh& seed);

  void scatterMesh(
    mesh::Mesh& seed);

  void computeBoundingMappings();

  void clearBoundingMappings();

  bool doesVertexContribute(int vertexID);

  // @brief Logging device.
  static tarch::logging::Log _log;

  std::string _accessorName;

  std::string _providerName;

  std::map<std::string,com::PtrCommunication> _receivers;

  int _dimensions;

  mapping::PtrMapping _boundingFromMapping;

  mapping::PtrMapping _boundingToMapping;
};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_COMMUNICATEDGEOMETRY_HPP_ */
