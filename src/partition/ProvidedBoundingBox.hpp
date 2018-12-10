#pragma once
#include "Partition.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace partition {


/**
 * @brief this class is supposed to:
 * 1- create bounding boxes around the mesh partition of each rank.
 * 2- gather these bounding boxes in the master
 * 3- send them to the other master
 */
class ProvidedBoundingBox :  public Partition 
{
public:

   /// Constructor
  ProvidedBoundingBox(mesh::PtrMesh mesh, bool hasToSend, double safetyFactor);

  virtual ~ProvidedBoundingBox() {}

  // These functions will be iplemented in 3rd package
  virtual void communicate();
  virtual void compute();
  virtual void createOwnerInformation();
  
  /// The boundingbox is gathered and sent to another participant (if required)
  virtual void communicateBoundingBox();

  /// The feedback from other participants received here and initial communication map is build 
  virtual void computeBoundingBox();

private:

  logging::Logger _log{"partition::ProvidedBoundingBox"};
  
//  virtual void createOwnerInformation();
    
  bool _hasToSend;
   
  int _dimensions;
  
  double _safetyFactor;
  
};

}} // namespace precice, partition
