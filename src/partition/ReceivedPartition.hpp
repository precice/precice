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

   /**
    * @brief Constructor
    */
  ReceivedPartition (bool filterFirst, int dimensions, double safetyFactor);

   virtual ~ReceivedPartition() {}

   /**
    * @brief The mesh is received from another participant.
    */
   virtual void communicate ();

   /**
    * @brief The mesh is re-partitioned and all distribution data structures are set up.
    */
   virtual void compute ();

private:

   std::vector<int> filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB);

   void prepareBoundingBox();

   bool isVertexInBB(const mesh::Vertex& vertex);

   bool _filterFirst;

   mesh::Mesh::BoundingBox _bb;

   int _dimensions;

   double _safetyFactor;

   static logging::Logger _log;

};

}} // namespace precice, partition
