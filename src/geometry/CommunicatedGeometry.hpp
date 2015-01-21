// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_COMMUNICATEDGEOMETRY_HPP_
#define PRECICE_GEOMETRY_COMMUNICATEDGEOMETRY_HPP_

#include "Geometry.hpp"
#include "m2n/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "utils/Dimensions.hpp"
#include "utils/MasterSlave.hpp"
#include "tarch/logging/Log.h"
#include <string>
#include <map>

namespace precice {
namespace geometry {

/**
 * @brief Creates the mesh by copying another remote mesh via communication. In case of parallel
 * solvers using a master-slave approach the meshes are gathered or scattered. In case of scattering,
 * globals meshes are sent first to each slave and filtered then depending on both mappings defined.
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
    m2n::PtrM2N m2n );

  void setBoundingFromMapping(mapping::PtrMapping mapping);

  void setBoundingToMapping(mapping::PtrMapping mapping);

protected:

  /**
   * @brief Is called from Geometry  and transmits the mesh.
   */
  virtual void specializedCreate ( mesh::Mesh& seed );

private:

  /**
   * @brief The received mesh is scattered amongst the slaves.
   */
  void scatterMesh(
    mesh::Mesh& seed);

  void sendMesh(
    mesh::Mesh& seed);

  void receiveMesh(
    mesh::Mesh& seed);

  /**
   * @brief Compute the preliminary mappings between the global mesh and the slave's own mesh.
   */
  void computeBoundingMappings();

  void clearBoundingMappings();

  void mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb);

  std::vector<int> filterMesh(mesh::Mesh& seed, mesh::Mesh& filteredMesh, bool filterByMapping);

  /**
   * @brief Returns true if a vertex contributes to one of the 2 mappings. If false, the vertex can be erased.
   */
  bool doesVertexContribute(const mesh::Vertex& vertex, bool filterByMapping);

  // @brief Logging device.
  static tarch::logging::Log _log;

  std::string _accessorName;

  std::string _providerName;

  std::map<std::string,m2n::PtrM2N> _receivers;

  int _dimensions;

  mapping::PtrMapping _boundingFromMapping;

  mapping::PtrMapping _boundingToMapping;

  mesh::Mesh::BoundingBox _bb;

  double _safetyGap;
};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_COMMUNICATEDGEOMETRY_HPP_ */
