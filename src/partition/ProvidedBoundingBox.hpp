#pragma once
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace partition {


/**
 * @brief this class is supposed to:
 * 1- gather bounding boxes around the mesh partition of each rank in the master
 * 2- send them to the other master
 * 3- receive the feedback map from the other master (i.e. who is connected to whom from the other participant's perspective)
 * 4- create own initial connection map (i.e. who is connected to whom from this participant's perspective)
 */
class ProvidedBoundingBox :  public Partition 
{
public:

   /// Constructor
  ProvidedBoundingBox(mesh::PtrMesh mesh, bool hasToSend, double safetyFactor);

  virtual ~ProvidedBoundingBox() {}
  
  /// The bounding box is gathered and sent to another participant (if required)
  virtual void communicateBoundingBox();

  /// The feedback from the other participant is received and the initial connection map is build
  virtual void computeBoundingBox();

  /// sends mesh partition to remote connected ranks
  virtual void communicate() override;

  /// receive final communication map from remote connected ranks
  virtual void compute() override;

  virtual void createOwnerInformation();

private:

  logging::Logger _log{"partition::ProvidedBoundingBox"};
  
  bool _hasToSend;
   
  int _dimensions;
  
  double _safetyFactor;
  
};

}} // namespace precice, partition
