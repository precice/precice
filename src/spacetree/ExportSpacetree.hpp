// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_EXPORTSPACETREE_HPP_
#define PRECICE_NEWSPACETREE_EXPORTSPACETREE_HPP_

#include "spacetree/Spacetree.hpp"
#include "boost/array.hpp"
#include "utils/Dimensions.hpp"
#include <string>
#include <list>
#include <cmath>

namespace precice {
namespace spacetree {

/**
 * @brief Writes the leave cells of a spacetree as voxels into a vtk file
 */
class ExportSpacetree : public Spacetree::Visitor
{
public:

   /**
    * @brief Constructor, defining vtk filename.
    */
   ExportSpacetree ( const std::string& filename );

   virtual ~ExportSpacetree() {}

   void doExport ( Spacetree& toExport );

   void nodeCallback (
     const utils::DynVector& center,
     const utils::DynVector& halflengths,
     int                     position );

   void leafCallback (
     const utils::DynVector& center,
     const utils::DynVector& halflengths,
     int                     position,
     const mesh::Group&      content );

//   /**
//    * @brief Exports the cells of the spacetree into a vtk file.
//    */
//   void doExport ( const Spacetree& spacetree );

private:

   std::string _filename;

   size_t _vertexCounter;

   std::list<utils::DynVector> _vertices;

   std::list<std::vector<int> > _cells;

   std::list<int> _cellPositions;

   std::list<int> _cellContents;
};

}} // namespace precice, spacetree

#endif // PRECICE_NEWSPACETREE_EXPORTSPACETREE_HPP_
