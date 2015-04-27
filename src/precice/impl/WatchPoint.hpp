// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_WATCHPOINT_HPP_
#define PRECICE_IMPL_WATCHPOINT_HPP_

#include "SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include <string>
#include <vector>

namespace precice {
  namespace mesh {
    class Vertex;
  }
}

// ----------------------------------------------------------- CLASS_DEFINITION

namespace precice {
namespace impl {

/**
 * @brief Observes and exports coordinates of a point on the geometry.
 */
class WatchPoint
{
public:

   /**
    * @brief Constructor.
    *
    * @param meshToWatch [IN] Mesh to be watched, can be empty on construction.
    */
   WatchPoint (
     const utils::DynVector& pointCoords,
     const mesh::PtrMesh&    meshToWatch,
     const std::string&      exportFilename );

   const mesh::PtrMesh& mesh() const;

   const std::string& filename() const;

   /**
    * @brief Initializes the watch point for exporting point data.
    */
   void initialize();

   /**
    * @brief Writes one line with data of the watchpoint into the output file.
    */
   void exportPointData(double time);

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   utils::DynVector _point;

   mesh::PtrMesh _mesh;

   io::TXTTableWriter _txtWriter;

   double _shortestDistance;

   std::vector<double> _weights;

   std::vector<mesh::Vertex*> _vertices;

   std::vector<mesh::PtrData> _dataToExport;

   // @bried Holds the information if this processor is the closest
   bool _isClosest;

   void getValue (
     utils::DynVector& value,
     mesh::PtrData&    data );

   void getValue (
     double&        value,
     mesh::PtrData& data );
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_WATCHPOINT_HPP_ */
