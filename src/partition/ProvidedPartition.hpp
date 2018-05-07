#pragma once

#include "Partition.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace partition {


/**
 * @brief A partition that is provided by the participant.
 *
 * The participant already provides a partition by calling setMeshVertices etc.
 * If required the mesh needs to be sent to another participant.
 * Furthermore, distribution data structures need to be set up.
 */
class ProvidedPartition : public Partition
{
public:

   /// Constructor
   ProvidedPartition (mesh::PtrMesh mesh, bool hasToSend);

   virtual ~ProvidedPartition() {}

   /// The mesh is gathered and sent to another participant (if required)
   virtual void communicate ();

   /// All distribution data structures are set up.
   virtual void compute ();

private:

   virtual void createOwnerInformation();

   logging::Logger _log{"partition::ProvidedPartition"};

   bool _hasToSend;

};

}} // namespace precice, partition
