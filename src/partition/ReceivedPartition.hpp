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

  /// Defines the typ of geometric filter used
  enum GeometricFilter {
    // @brief undefined
    UNDEFINED,
    // @brief No geometric filter used (e.g. for RBF mappings)
    NO_FILTER,
    // @brief Filter at master and communicate only filtered mesh.
    FILTER_FIRST,
    // @brief Broadcast first and filter then
    BROADCAST_FILTER
  };

   /// Constructor
   ReceivedPartition (mesh::PtrMesh mesh, GeometricFilter geometricFilter, double safetyFactor);

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

   GeometricFilter _geometricFilter;

   mesh::Mesh::BoundingBox _bb;

   int _dimensions;

   double _safetyFactor;

   logging::Logger _log{"partition::ReceivedPartition"};

};

}} // namespace precice, partition
