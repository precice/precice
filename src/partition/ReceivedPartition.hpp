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
 * A mesh is received by the master rank and re-partitioned among all slave ranks.
 * Afterwards necessary distribution data structures are set up.
 */
class ReceivedPartition : public Partition
{
public:

   /// Constructor
   ReceivedPartition (bool filterFirst, int dimensions, double safetyFactor);

   virtual ~ReceivedPartition() {}

   /// The mesh is received from another participant.
   virtual void communicate ();

   /// The mesh is re-partitioned and all distribution data structures are set up.
   virtual void compute ();

private:

   void filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB);

   void prepareBoundingBox();

   bool isVertexInBB(const mesh::Vertex& vertex);

   virtual void createOwnerInformation();

   /// Helper function for 'createOwnerFunction' to set local owner information
   void setOwnerInformation(const std::vector<int> &ownerVec);

   bool _filterFirst;

   mesh::Mesh::BoundingBox _bb;

   int _dimensions;

   double _safetyFactor;

   static logging::Logger _log;

};

}} // namespace precice, partition
