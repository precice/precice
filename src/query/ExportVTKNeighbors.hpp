// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_EXPORTVTKNEIGHBORS_HPP_
#define PRECICE_QUERY_EXPORTVTKNEIGHBORS_HPP_

#include "FindClosest.hpp"
#include "utils/Dimensions.hpp"
#include <string>

namespace precice {
namespace query {

/*
 * @brief Visualizes results of FindClosestOnMeshVisitor searches
 *
 * The visualization consists of lines connecting search points with
 * closest points on found visitables. The lines have assigned a value
 * distance, that has the value of the distance to the found visitables.
 *
 * @author Bernhard Gatzhammer
 */
class ExportVTKNeighbors
{
public:

   /**
    * @brief Adds a neighborhood relation already computed.
    *
    * @param queryPoint [IN] Point whose next neighbor was determined.
    * @param neighborPoint [IN] ClosestNeighbor of queryPoint.
    */
   void addNeighbors (
     const utils::DynVector& queryPoint,
     const ClosestElement&   closestNeighbor );

   /**
    * @brief Visualizes results of performed searches into vtk file
    */
   void exportNeighbors ( const std::string& filename );

   /**
    * @brief Resets results of performed searches
    */
   void resetElements() { _neighbors.clear(); };

private:

   // @brief Results of performed searches
   std::vector<std::pair<utils::DynVector,query::ClosestElement> >   _neighbors;
};

}} // namespace precice, query

#endif /* PRECICE_QUERY_EXPORTVTKNEIGHBORS_HPP_ */
