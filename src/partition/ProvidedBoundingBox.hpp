#pragma once
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"



namespace PartitionTests {
namespace ProvidedBoundingBoxTests {
struct TestProvidedBoundingBox2D;
struct TestProvidedBoundingBox3D;
struct TestInitialCommunicationMap;
struct TestM2NMeshExchange;
}
}

namespace PartitionTests {
namespace ReceivedBoundingBoxTests {
struct TestReceivedBoundingBox2D;
struct TestReceivedBoundingBox3D;
}
}


namespace precice {
namespace partition {


/**
 * @brief this class is supposed to:
 * 1- creat bounding boxes around each ranks mesh partition
 * 2- gather these bounding boxes in the master
 * 3- send them to the other master
 */
class ProvidedBoundingBox :  public Partition 
{
public:

   /// Constructor
  ProvidedBoundingBox(mesh::PtrMesh mesh, bool hasToSend, double safetyFactor);

  virtual ~ProvidedBoundingBox() {}

   /// The boundingbox is gathered and sent to another participant (if required)
  virtual void communicate();
  virtual void compute();
  // virtual void communicateBoundingBox();
  // virtual void computeBoundingBox();
  
  friend struct PartitionTests::ProvidedBoundingBoxTests::TestProvidedBoundingBox2D;
  friend struct PartitionTests::ProvidedBoundingBoxTests::TestProvidedBoundingBox3D;
  friend struct PartitionTests::ProvidedBoundingBoxTests::TestInitialCommunicationMap;
  friend struct PartitionTests::ProvidedBoundingBoxTests::TestM2NMeshExchange;
  friend struct PartitionTests::ReceivedBoundingBoxTests::TestReceivedBoundingBox2D;
  friend struct PartitionTests::ReceivedBoundingBoxTests::TestReceivedBoundingBox3D;
    
private:
  
//  virtual void createOwnerInformation();
  
  static logging::Logger _log;
  
  bool _hasToSend;
   
  int _dimensions;
  
  std::vector<int> _connectedRanks;
  
  double _safetyFactor;
  
  std::vector<int> _vertexCounters;
  
  int _remoteParComSize = 0 ;

  mesh::Mesh::BoundingBox _bb;
  
  mesh::Mesh::BoundingBoxMap _globalBB;

  mesh::Mesh::FeedbackMap _receivedFeedbackMap;
  
};

}} // namespace precice, partition
