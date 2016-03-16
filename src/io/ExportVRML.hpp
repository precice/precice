// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_EXPORTVRML_HPP_
#define PRECICE_IO_EXPORTVRML_HPP_

#include "Export.hpp"
#include "logging/Logger.hpp"
#include <string>
#include <map>

namespace precice {
  namespace mesh {
     class Vertex;
     class Edge;
     class Triangle;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief Writes polygonal or triangulated meshes to vrml 1.0 files.
 */
class ExportVRML : public Export
{
public:

  /**
   * @brief Standard constructor
   *
   * @param container [IN] Container holding geometry to be visualized
   */
  ExportVRML ( bool plotNormals );

  /**
   * @brief Returns export type ID, i.e., VRML.
   */
  virtual int getType() const;

  /**
   * @brief Perform writing to VRML file
   */
  virtual void doExport (
    const std::string& filname,
    mesh::Mesh&        mesh );

  /**
   * @brief Perform writing to VRML file of a full geometry checkpoint.
   *
   * A full checkpoint also  includes vertex data and property containers.
   */
  void doExportCheckpoint (
    const std::string& filename,
    mesh::Mesh&        mesh );

private:

  /// @brief Logging device.
  static logging::Logger _log;

  void openFile (
    std::ofstream&     outFile,
    const std::string& filename ) const;

  void writeHeader ( std::ofstream& outFile ) const;

  void writeGeometry (
    std::ofstream& outFile,
    mesh::Mesh&    mesh ) const;

  void writeVertexData (
    std::ofstream& outFile,
    mesh::Mesh&    mesh ) const;

  void writePropertyContainer (
    std::ofstream& outFile,
    mesh::Mesh&    mesh ) const;
};

}} // namespace precice, io

#endif /* PRECICE_IO_EXPORTVRML_HPP_ */
