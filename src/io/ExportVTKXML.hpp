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
  //virtual int getType() const;

  /**
   * @brief Writes the master file
   */
  void ExportVTKXML::writeMasterFile
  (
    const std::string& filename,
    mesh::Mesh&        mesh);

  /**
   * @brief Perform writing to vtk file
   */
  virtual void doExport (
    const std::string& filename,
    mesh::Mesh&        mesh );

  static void initializeWriting (
    const std::string& filename,
    std::ofstream&     filestream );

  static void writeHeader ( std::ostream& outFile );

  static void writeVertex (
    const utils::DynVector& position,
    std::ostream&           outFile );

  static void writeLine (
    int           vertexIndices[2],
    std::ostream& outFile );

  static void writeTriangle (
    int           vertexIndices[3],
    std::ostream& outFile );

  static void writeQuadrangle (
    int           vertexIndices[4],
    std::ostream& outFile );

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief By default set true: plot vertex normals, false: no normals plotting
   bool _writeNormals;

   // @brief true: write as parallel file, false: write as serial file
   bool _parallelWrite;

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
