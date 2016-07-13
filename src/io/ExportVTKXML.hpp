// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_EXPORTVTKXML_HPP_
#define PRECICE_IO_EXPORTVTKXML_HPP_

#include "Export.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include "utils/Dimensions.hpp"
#include <string>

namespace precice {
   namespace mesh {
      class Mesh;
      class Quad;
      class Edge;
      class Triangle;
   }
}

namespace precice {
   namespace utils {
      class Parallel;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief Writes polygonal, triangle, or quadrangle meshes to xml-vtk files. Works in serial and parallel.
 */
class ExportVTKXML : public Export
{
public:

  /**
   * @brief Standard constructor
   *
   * @param exportNormals  [IN] boolean: write normals to file?
   * @param parallelWrite  [IN] boolean: write as parallel?
   */
  ExportVTKXML ( bool writeNormals, bool parallelWrite );

  /**
   * @brief Returns the VTK type ID.
   */
  virtual int getType() const;


  /**
   * @brief Perform writing to vtk file
   */
  virtual void doExport (
    const std::string& filename,
    mesh::Mesh&        mesh );

  static void writeVertex (
    const utils::DynVector& position,
    std::ofstream&           outFile );

  static void writeLine (
    mesh::Edge&    edge,
    std::ofstream& outFile );

  static void writeTriangle (
    mesh::Triangle& triangle,
    std::ofstream&  outFile );

  static void writeQuadrangle (
    mesh::Quad&    quad,
    std::ofstream& outFile );

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief By default set true: plot vertex normals, false: no normals plotting
   bool _writeNormals;

   // @brief true: write as parallel file, false: write as serial file
   bool _parallelWrite;

   // @brief true: data names and dimensions have been processed
   bool _isDataNamesAndDimensionsProcessed;

   // @brief true: mesh contains cells
   bool _isCellPresent;

   // @ brief dimensions of mesh
   int _meshDimensions;

   // @brief List of names of all scalar data on mesh
   std::vector<std::string> _scalarDataNames;

   // @brief List of names of all vector data on mesh
   std::vector<std::string> _vectorDataNames;

   // @brief List of vector data dimensions
   std::vector<int> _vectorDataDimensions;
   /**
    * @brief Stores scalar and vector data names and dimensions in string vectors
    * Needed for writing master file and sub files
    */
   void processDataNamesAndDimensions
   (
     mesh::Mesh& mesh);

   /**
    * @brief Writes the master file (called only by rank 0)
    */
   void writeMasterFile
   (
     const std::string& filename,
     mesh::Mesh&        mesh);

   /**
    * @brief Writes the sub file for each process
    */
   void writeSubFile
   (
     const std::string& filename,
 	 mesh::Mesh&        mesh);

   void openFile (
    std::ofstream&     outFile,
    const std::string& filename ) const;

   void exportGeometry (
     std::ofstream& outFile,
     mesh::Mesh&    mesh );

   void exportData (
     std::ofstream& outFile,
     mesh::Mesh&    mesh );
};

}} // namespace precice, io

#endif /* PRECICE_IO_EXPORTVTKXML_HPP_ */
