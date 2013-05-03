// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_EXPORTVTKVOXELQUERIES_HPP_
#define PRECICE_QUERY_EXPORTVTKVOXELQUERIES_HPP_

#include "utils/Dimensions.hpp"
#include <vector>
#include <string>

namespace precice {
namespace mesh {
   class Group;
   class IVisitable;
}
namespace spacetree {
   class Spacetree;
}
namespace query {

/**
 * @brief Visualizes results of FindVoxelContainedVisitor searches
 *
 * The results are visualized as voxels.
 *
 * @author Bernhard Gatzhammer
 */
class ExportVTKVoxelQueries
{
public:

   void addQuery (
      const utils::DynVector& voxelCenter,
      const utils::DynVector& voxelHalflengths,
      int                     containedElementsCount );

   /**
    * @brief Writes performed queries' results into vtk file
    */
   void exportQueries ( std::string filename );

//   void exportVoxels ( std::string filename );

   /**
    * @brief Resets query results
    */
   void resetQueries();

private:

   // @brief Stores voxel centers queried
   std::vector<utils::DynVector> _voxelCenters;

   // @brief Stores voxel halflengths queried
   std::vector<utils::DynVector> _voxelHalflengths;

   // @brief Stores in queried voxels contained visitables
   std::vector<int> _containedElementCount;
};

}}

#endif /* PRECICE_QUERY_EXPORTVTKVOXELQUERIES_HPP_ */
