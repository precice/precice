#pragma once

#include <vector>
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace partition {
/**
 * @brief this class is supposed to:
 * 1- receive a set of bounding boxes from other solver
 * 2- each rank compares its bounding box with received set of boundingboxs
 * 3- if there is intersection, local rank is added to the connection map
 * 4- connection map sent to the other participant  
 * @todo add documentation
 */
class ReceivedBoundingBox : public Partition {
public:
  /// Constructor
  ReceivedBoundingBox(mesh::PtrMesh mesh, double safetyFactor);
  virtual ~ReceivedBoundingBox() {}

  /// bounding box map is received from other participant
  virtual void communicateBoundingBox();

  /// bounding boxes are compared and feedback sent to master of other participant
  virtual void computeBoundingBox();

  /// These functions will be implemented in 3rd package
  virtual void communicate();
  virtual void compute();

private:
  logging::Logger _log{"partition::ReceivedBoundingBox"};

  /// compares to bounding box and if they have intersection, returns true, otherwise flase!
  static bool overlapping(mesh::Mesh::BoundingBox const &currentBB, mesh::Mesh::BoundingBox const &receivedBB);

  /// Sets _bb to the union with the mesh from fromMapping resp. toMapping, also enlage by _safetyFactor
  void prepareBoundingBox();

  /// will be implemented in 3rd work package
  virtual void createOwnerInformation();

  /// number of other particpant ranks
  int _remoteParComSize = 0;

  /// bounding box map of other participant
  mesh::Mesh::BoundingBoxMap _remoteBBM;

  mesh::Mesh::BoundingBox _bb;

  int _dimensions;

  double _safetyFactor;
};

} // namespace partition
} // namespace precice
