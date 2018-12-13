#pragma once

#include "Partition.hpp"
#include "logging/Logger.hpp"
#include <vector>
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"


namespace precice {
namespace partition {
/**
 * @brief A partition that is computed from a mesh received from another participant.
 *
 * @todo add documentation
 */
class ReceivedBoundingBox : public Partition
{
public:

   /// Constructor
  ReceivedBoundingBox (mesh::PtrMesh mesh, double safetyFactor);
  virtual ~ReceivedBoundingBox() {}

  /// bounding box map is received from other participant
  virtual void communicateBoundingBox();

  /// bounding boxes are compared and feedback sent to master of other participant
  virtual void computeBoundingBox();
  
  /// These functions will be implemented in 3rd package
  virtual void communicate ();
  virtual void compute (); 
    
private:

  logging::Logger _log{"partition::ReceivedBoundingBox"};
  
  bool overlapping(mesh::Mesh::BoundingBox currentBB, mesh::Mesh::BoundingBox receivedBB);

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

}} // namespace precice, partition
