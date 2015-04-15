#pragma once

#include "Geometry.hpp"
#include "m2n/M2N.hpp"
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
    const std::string& receiver,
    m2n::M2N::SharedPointer m2n );

  void setBoundingFromMapping(mapping::PtrMapping mapping);

  void setBoundingToMapping(mapping::PtrMapping mapping);

  void setSafetyFactor(double safetyFactor);

protected:

  /**
   * Is called from Geometry and sends the mesh if the accessor is provider,
   * receives if the accessor is contained in receivers, fails otherwise.
   */
  void specializedCreate ( mesh::Mesh& seed );

private:

  /// The received mesh is scattered amongst the slaves.
  void scatterMesh(
    mesh::Mesh& seed);

  void sendMesh(
    mesh::Mesh& seed);

  void receiveMesh(
    mesh::Mesh& seed);

  /// Compute the preliminary mappings between the global mesh and the slave's own mesh.
  void computeBoundingMappings();

  void clearBoundingMappings();

  void mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb);

  std::vector<int> filterMesh(mesh::Mesh& seed, mesh::Mesh& filteredMesh, bool filterByMapping);

  /// Returns true if a vertex contributes to one of the 2 mappings. If false, the vertex can be erased.
 bool doesVertexContribute(const mesh::Vertex& vertex, bool filterByMapping);

  /// Logging device.
  static tarch::logging::Log _log;

  std::string _accessorName;

  std::string _providerName;

  std::map<std::string,m2n::M2N::SharedPointer> _receivers;

  int _dimensions;

  mapping::PtrMapping _boundingFromMapping;

  mapping::PtrMapping _boundingToMapping;

  mesh::Mesh::BoundingBox _bb;

  double _safetyGap;

  double _safetyFactor;
};

}} // namespace precice, geometry
