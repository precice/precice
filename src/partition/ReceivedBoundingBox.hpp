#pragma once

#include "Partition.hpp"
#include "logging/Logger.hpp"
#include <vector>
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"


// Forward delcration to friend the boost test struct

namespace precice {
namespace partition {
/**
 * @brief A partition that is computed from a mesh received from another participant.
 *
 * A mesh is received by the master rank and re-partitioned among all slave ranks.
 * Afterwards necessary distribution data structures are set up.
 */
class ReceivedBoundingBox : public Partition
{
public:
  
  /// Constructor
  ReceivedBoundingBox (mesh::PtrMesh mesh, double safetyFactor, GeometricFilter geometricFilter);
  virtual ~ReceivedBoundingBox() {}
  /// The mesh is received from another participant.
//  virtual void communicate ();
  /// The mesh is re-partitioned and all distribution data structures are set up.
//  virtual void compute ();
  virtual void communicateBoundingBox();
  virtual void computeBoundingBox();  
  
private:

  logging::Logger _log{"partition::ReceivedBoundingBox"};
  
  bool compareBoundingBox(mesh::Mesh::BoundingBox currentBB, mesh::Mesh::BoundingBox receivedBB);

  // void filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB);

  void prepareBoundingBox();

  // bool isVertexInBB(const mesh::Vertex& vertex);

  // virtual void createOwnerInformation();

  // /// Helper function for 'createOwnerFunction' to set local owner information
  // void setOwnerInformation(const std::vector<int> &ownerVec);

  // number of other particpant ranks
  int _remoteParComSize = 0;
  
  mesh::Mesh::BoundingBoxMap _globalBB;
  
  mesh::Mesh::BoundingBox _bb;

  // feedback from each rank which contains list of connected ranks of the other solver
  std::vector<int> _feedback; 

  // int : each rank, vect: connected ranks
  mesh::Mesh::FeedbackMap  _feedbackMap; 

  // list of connected ranks -> list of vertices
  mesh::Mesh::FeedbackMap _localCommunicationMap; 

  int _dimensions;

  double _safetyFactor;

  int _numberOfVertices = 0;  
  
};

}} // namespace precice, partition
